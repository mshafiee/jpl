package jpl

import (
	"encoding/binary"
	"errors"
	"fmt"
	"math"
	"os"
	"strings"
	"unsafe"
)

// CelestialBody represents various celestial bodies as integer constants.
type CelestialBody int

const (
	Mercury CelestialBody = iota
	Venus
	Earth
	Mars
	Jupiter
	Saturn
	Uranus
	Neptune
	Pluto
	Moon
	Sun
	SolarSystemBarycenter
	EarthMoonBarycenter
	Nutations
	Librations
)

// NotAvailable is used to indicate unavailable or undefined data.
const NotAvailable = -1

// JPL represents the local globals for the JPL ephemerides reader
type JPL struct {
	JplFilePath string   // Path to the JPL data file
	JplFile     *os.File // File pointer to the JPL data file
	DoReorder   bool     // Flag to determine if byte reordering is necessary due to different endianess
	Constants   struct { // Struct holding constants from the JPL ephemeris files
		Cval  [400]float64 // Holds various constants from the file
		SS    [3]float64   // Start, stop, and step times for the ephemeris
		AU    float64      // Astronomical unit in kilometers
		Emrat float64      // Earth/Moon mass ratio
		Denum int32        // Ephemeris number (version)
		Ncon  int32        // Number of constants
		IPT   [39]int32    // Array holding pointers into the coefficient arrays
	}
	ChCnam    [6 * 400]byte // Names of constants
	PV        [78]float64   // Position and velocity array for planets
	PVSun     [6]float64    // Position and velocity of the Sun
	Buf       [1500]float64 // Buffer to hold coefficients from file
	Chebyshev struct {      // Struct for Chebyshev coefficients
		PC [18]float64 // Position coefficients
		VC [18]float64 // Velocity coefficients
		AC [18]float64 // Acceleration coefficients
		JC [18]float64 // Jerk coefficients
	}
	DoKm      bool // Determines if output is in km (true) or AU (false)
	DebugMode bool // Flag to enable debug mode which may print additional info
}

// NewJPL creates a new JPL structure with default values
func NewJPL(filepath string) *JPL {
	// Automatically determine byte order to set DoReorder
	var x uint32 = 0x01020304
	doReorder := (*[4]byte)(unsafe.Pointer(&x))[0] != 0x04

	return &JPL{
		JplFilePath: filepath,
		DoReorder:   doReorder,
		DoKm:        false,
	}
}

// handleEarthMoonInteraction handles specific interactions between Earth and Moon in the pv array.
func (js *JPL) handleEarthMoonInteraction(list []int32, pv []float64) {
	if list[Earth] == 2 {
		for i := 0; i < 6; i++ {
			pv[i+6*int(Earth)] -= pv[i+6*int(Moon)] / (js.Constants.Emrat + 1.0)
		}
	}
	if list[Moon] == 2 {
		for i := 0; i < 6; i++ {
			pv[i+6*int(Moon)] += pv[i+6*int(Earth)]
		}
	}
}

// adjustPositionsForSunEMBBary adjusts positions in the pv array for the Sun, EarthMoonBarycenter, and barycenter.
func adjustPositionsForSunEMBBary(ntarg, ncent CelestialBody, pv, pvsun []float64) {
	if ntarg == Sun || ncent == Sun {
		for i := 0; i < 6; i++ {
			pv[i+6*int(Sun)] = pvsun[i]
		}
	}
	if ntarg == SolarSystemBarycenter || ncent == SolarSystemBarycenter {
		for i := 0; i < 6; i++ {
			pv[i+6*int(SolarSystemBarycenter)] = 0.0
		}
	}
	if ntarg == EarthMoonBarycenter || ncent == EarthMoonBarycenter {
		for i := 0; i < 6; i++ {
			pv[i+6*int(EarthMoonBarycenter)] = pv[i+6*int(Earth)]
		}
	}
}

// setupListForState sets up the list array based on the target and center.
func setupListForState(ntarg, ncent CelestialBody) []int32 {
	list := make([]int32, 12)

	if ntarg < Sun {
		list[ntarg] = 2
	}
	if ntarg == Moon || ntarg == Earth || ntarg == EarthMoonBarycenter {
		list[Earth] = 2
	}
	if ntarg == Earth {
		list[Moon] = 2
	}
	if ncent < Sun {
		list[ncent] = 2
	}
	if ncent == Moon || ncent == Earth || ncent == EarthMoonBarycenter {
		list[Earth] = 2
	}
	if ncent == Earth {
		list[Moon] = 2
	}

	return list
}

// handleLibration handles the case when ntarg is Librations (libration).
// Returns 0 if successful, otherwise returns an error code.
func (js *JPL) handleLibration(et float64, pv, pvsun []float64) ([]float64, error) {
	rrd := make([]float64, 6)
	if js.Constants.IPT[37] > 0 {
		list := make([]int32, 12)
		list[11] = 2
		if err := js.state(et, list, false, pv, pvsun, rrd); err != nil {
			return nil, err
		}
		copy(rrd, pv[60:66])
		return rrd, nil
	} else {
		return nil, fmt.Errorf("no librations on the ephemeris file")
	}
}

// handleNutation handles the case when ntarg is Nutations (nutation).
// Returns 0 if successful, otherwise returns an error code.
func (js *JPL) handleNutation(et float64, pv, pvsun []float64) ([]float64, error) {
	rrd := make([]float64, 6)
	if js.Constants.IPT[34] > 0 {
		list := make([]int32, 12)
		list[10] = 2
		err := js.state(et, list, false, pv, pvsun, rrd)
		return rrd, err
	} else {
		return nil, fmt.Errorf("no nutations on the JPL ephemeris file")
	}
}

// calculateKsize calculates the record size (ksize) based on the planetary interpolation parameters.
func calculateKsize(ipt []int32, numde int32) int32 {
	kmx, khi := int32(0), int32(0)
	for i := 0; i < 13; i++ {
		if ipt[i*3] > kmx {
			kmx = ipt[i*3]
			khi = int32(i + 1)
		}
	}
	nd := int32(3)
	if khi == 12 {
		nd = 2
	}
	ksize := (ipt[khi*3-3] + nd*ipt[khi*3-2]*ipt[khi*3-1] - 1) * 2

	if ksize == 1546 {
		ksize = 1652
	}

	return ksize
}

// readConstantsAndParameters reads the constants and parameters from the JPL file.
func (js *JPL) readConstantsAndParameters(au *float64, emrat *float64, lpt []int32, numde *int32, ncon *int32) error {
	var buf [8]byte

	if _, err := js.JplFile.Read(buf[:4]); err != nil {
		return err
	}
	*ncon = int32(binary.LittleEndian.Uint32(buf[:4]))
	if js.DoReorder {
		Reorder(buf[:4])
		*ncon = int32(binary.BigEndian.Uint32(buf[:4]))
	}

	if _, err := js.JplFile.Read(buf[:8]); err != nil {
		return err
	}
	*au = float64(binary.LittleEndian.Uint64(buf[:8]))
	if js.DoReorder {
		Reorder(buf[:8])
		*au = float64(binary.BigEndian.Uint64(buf[:8]))
	}

	if _, err := js.JplFile.Read(buf[:8]); err != nil {
		return err
	}
	*emrat = float64(binary.LittleEndian.Uint64(buf[:8]))
	if js.DoReorder {
		Reorder(buf[:8])
		*emrat = float64(binary.BigEndian.Uint64(buf[:8]))
	}

	for i := 0; i < 36; i++ {
		if _, err := js.JplFile.Read(buf[:4]); err != nil {
			return err
		}
		js.Constants.IPT[i] = int32(binary.LittleEndian.Uint32(buf[:4]))
		if js.DoReorder {
			Reorder(buf[:4])
			js.Constants.IPT[i] = int32(binary.BigEndian.Uint32(buf[:4]))
		}
	}

	if _, err := js.JplFile.Read(buf[:4]); err != nil {
		return err
	}
	*numde = int32(binary.LittleEndian.Uint32(buf[:4]))
	if js.DoReorder {
		Reorder(buf[:4])
		*numde = int32(binary.BigEndian.Uint32(buf[:4]))
	}

	for i := 0; i < 3; i++ {
		if _, err := js.JplFile.Read(buf[:4]); err != nil {
			return err
		}
		lpt[i] = int32(binary.LittleEndian.Uint32(buf[:4]))
		if js.DoReorder {
			Reorder(buf[:4])
			lpt[i] = int32(binary.BigEndian.Uint32(buf[:4]))
		}
		js.Constants.IPT[i+36] = lpt[i]
	}

	if _, err := js.JplFile.Seek(0, 0); err != nil {
		return err
	}

	return nil
}

// readAndValidateHeader reads and validates the header information of the JPL file.
// Returns an error if unsuccessful.
func (js *JPL) readAndValidateHeader() error {
	var ttl [6 * 14 * 3]byte
	var ss [3]float64

	// Read the TTL array
	nrd, err := js.JplFile.Read(ttl[:])
	if err != nil || nrd != 252 {
		return errors.New("failed to read TTL array")
	}

	// Read the CNAM array
	nrd, err = js.JplFile.Read(js.ChCnam[:])
	if err != nil || nrd != 6*400 {
		return errors.New("failed to read CNAM array")
	}

	// Read the SS array
	for i := 0; i < 3; i++ {
		err = binary.Read(js.JplFile, binary.LittleEndian, &ss[i])
		if err != nil {
			return errors.New("failed to read SS array")
		}
	}

	js.DoReorder = ss[2] < 1 || ss[2] > 200

	// Store the values in SS, applying reordering if necessary
	for i := 0; i < 3; i++ {
		js.Constants.SS[i] = ss[i]
	}

	if js.DoReorder {
		Reorder(js.Constants.SS[:])
	}

	// Validate the SS values
	if js.Constants.SS[0] < -5583942 || js.Constants.SS[1] > 9025909 || js.Constants.SS[2] < 1 || js.Constants.SS[2] > 200 {
		return fmt.Errorf("alleged ephemeris file (%s) has an invalid format", js.JplFilePath)
	}

	return nil
}

// openJPLFile opens the JPL file and checks for errors.
// Returns 0 if successful, NotAvailable (-1) if an error occurs.
func (js *JPL) openJPLFile() error {
	// Attempt to open the file
	file, err := os.Open(js.JplFilePath)
	if err != nil {
		return errors.New("JPL file not available")
	}

	// Store the file pointer in the struct
	js.JplFile = file
	return nil
}

// fsizer opens the JPL file, reads the first record, and uses the info to compute ksize.
// Returns the calculated ksize or NotAvailable in case of an error.
func (js *JPL) fsizer() (int32, error) {
	var au, emrat float64
	var lpt [3]int32
	var numde, ncon int32

	err := js.openJPLFile()
	if err != nil {
		return NotAvailable, err
	}

	err = js.readAndValidateHeader()
	if err != nil {
		return NotAvailable, err
	}

	err = js.readConstantsAndParameters(&au, &emrat, lpt[:], &numde, &ncon)
	if err != nil {
		return NotAvailable, err
	}

	ksize := calculateKsize(js.Constants.IPT[:], numde)
	if ksize < 1000 || ksize > 5000 {
		return NotAvailable, fmt.Errorf("JPL ephemeris file does not provide valid ksize (%d)", ksize)
	}

	return ksize, nil
}

// EphemerisLookup reads the JPL planetary ephemeris data and computes the position and velocity
// of a specified celestial body ('ntarg') relative to another celestial body ('ncent') at a given
// ephemeris time ('et').
//
// Parameters:
//   - et (float64): The Julian ephemeris date for which the position and velocity are requested.
//   - ntarg (CelestialBody): The target celestial body whose position and velocity are desired.
//   - ncent (CelestialBody): The reference celestial body relative to which the position and velocity
//     of 'ntarg' are computed.
//
// Returns:
//   - []float64: A slice containing six elements: the x, y, z position components and the
//     dx, dy, dz velocity components of 'ntarg' relative to 'ncent'.
//   - error: Returns an error if the ephemeris lookup fails.
//
// Notes:
//   - If 'ntarg' and 'ncent' are the same, the function returns an empty slice as they are at the
//     same position.
//   - The function handles special cases for nutation (J_NUT) and libration (J_LIB) using
//     `handleNutation` and `handleLibration`, respectively.
//   - For other celestial bodies, the function sets up the state list for interpolation and
//     computes the positions and velocities using `state`.
//   - Adjustments are made for the Sun-Earth-Moon barycenter (EMB) as needed.
//   - In specific cases, such as the interaction between Earth and Moon, additional handling
//     is performed using `handleEarthMoonInteraction`.
func (js *JPL) EphemerisLookup(et float64, ntarg, ncent CelestialBody) ([]float64, error) {
	if ntarg == ncent {
		return []float64{}, nil
	}

	switch ntarg {
	case Nutations:
		return js.handleNutation(et, js.PV[:], js.PVSun[:])
	case Librations:
		return js.handleLibration(et, js.PV[:], js.PVSun[:])
	}

	rrd := make([]float64, 6)
	list := setupListForState(ntarg, ncent)
	if err := js.state(et, list, true, js.PV[:], js.PVSun[:], rrd); err != nil {
		return []float64{}, err
	}

	adjustPositionsForSunEMBBary(ntarg, ncent, js.PV[:], js.PVSun[:])

	if !(ntarg == Earth && ncent == Moon) && !(ntarg == Moon && ncent == Earth) {
		js.handleEarthMoonInteraction(list, js.PV[:])
	}

	for i := int32(0); i < 6; i++ {
		rrd[i] = js.PV[i+int32(ntarg)*6] - js.PV[i+int32(ncent)*6]
	}

	return rrd, nil
}

// computeSubInterval computes the sub-interval number and normalized Chebyshev time.
//
// Parameters:
//   - t: fractional time in the interval covered by coefficients.
//   - na: number of sets of coefficients in the full array.
//   - tc: pointer to store the normalized Chebyshev time (-1 <= tc <= 1).
//   - ni: pointer to store the sub-interval number.
//
// Returns:
//   - 0 on success.
func computeSubInterval(t float64, na int32, tc *float64, ni *int32) int32 {
	var dt1 float64
	if t >= 0 {
		dt1 = math.Floor(t)
	} else {
		dt1 = -math.Floor(-t)
	}
	temp := float64(na) * t
	*ni = int32(temp - dt1)
	*tc = (math.Mod(temp, 1.0)+dt1)*2.0 - 1.0
	return 0
}

// evaluatePolynomials evaluates Chebyshev polynomials at the given time.
//
// Parameters:
//   - pc: slice to store the evaluated polynomial values.
//   - tc: normalized Chebyshev time.
//   - ncf: number of coefficients per component.
//   - twot: pointer to store the double of tc.
func evaluatePolynomials(pc []float64, tc float64, ncf int32, twot *float64) {
	*twot = tc + tc
	pc[1] = tc
	for i := int32(2); i < ncf; i++ {
		pc[i] = (*twot)*pc[i-1] - pc[i-2]
	}
}

// interpolatePosition interpolates the position for each component using Chebyshev coefficients.
//
// Parameters:
//   - pv: slice to store the interpolated position values.
//   - pc: slice of evaluated polynomial values.
//   - buf: slice of Chebyshev coefficients of position.
//   - ncf: number of coefficients per component.
//   - ncm: number of components per set of coefficients.
//   - ni: sub-interval number for this set of coefficients.
func interpolatePosition(pv []float64, pc []float64, buf []float64, ncf, ncm, ni int32) {
	for i := int32(0); i < ncm; i++ {
		pv[i] = 0.0
		for j := ncf - 1; j >= 0; j-- {
			pv[i] += pc[j] * buf[j+(i+ni*ncm)*ncf]
		}
	}
}

// interpolateVelocity interpolates the velocity for each component using Chebyshev coefficients.
//
// Parameters:
//   - pv: slice to store the interpolated velocity values.
//   - vc: slice of evaluated velocity polynomial values.
//   - pc: slice of evaluated polynomial values.
//   - buf: slice of Chebyshev coefficients of position.
//   - ncf: number of coefficients per component.
//   - ncm: number of components per set of coefficients.
//   - ni: sub-interval number for this set of coefficients.
//   - bma: factor for converting Chebyshev time to real time.
func interpolateVelocity(pv []float64, vc, pc [18]float64, buf []float64, ncf, ncm, ni int32, bma float64) {
	for i := int32(0); i < ncm; i++ {
		pv[i+ncm] = 0.0
		for j := int32(ncf) - 1; j >= 1; j-- { // Start from the second coefficient (index 1)
			pv[i+ncm] += vc[j] * buf[j+(i+ni*ncm)*ncf]
		}
		pv[i+ncm] *= bma // Apply the time conversion factor
	}
}

// interpolateAcceleration interpolates the acceleration for each component using Chebyshev coefficients.
//
// Parameters:
//   - pv: slice to store the interpolated acceleration values.
//   - ac: slice of evaluated acceleration polynomial values.
//   - buf: slice of Chebyshev coefficients of position.
//   - ncf: number of coefficients per component.
//   - ncm: number of components per set of coefficients.
//   - ni: sub-interval number for this set of coefficients.
//   - bma2: factor for converting Chebyshev time to real time, squared.
func interpolateAcceleration(pv []float64, ac [18]float64, buf []float64, ncf, ncm, ni int32, bma2 float64) {
	for i := int32(0); i < ncm; i++ {
		pv[i+int32(ncm)*2] = 0.0
		for j := int32(ncf) - 1; j >= 2; j-- { // Start from the third coefficient (index 2)
			pv[i+int32(ncm)*2] += ac[j] * buf[j+(i+ni*ncm)*ncf]
		}
		pv[i+int32(ncm)*2] *= bma2 // Apply the squared time conversion factor
	}
}

// interpolateJerk interpolates the jerk for each component using Chebyshev coefficients.
//
// Parameters:
//   - pv: slice to store the interpolated jerk values.
//   - jc: slice of evaluated jerk polynomial values.
//   - buf: slice of Chebyshev coefficients of position.
//   - ncf: number of coefficients per component.
//   - ncm: number of components per set of coefficients.
//   - ni: sub-interval number for this set of coefficients.
//   - bma3: factor for converting Chebyshev time to real time, cubed.
func interpolateJerk(pv []float64, jc [18]float64, buf []float64, ncf, ncm, ni int32, bma3 float64) {
	for i := int32(0); i < ncm; i++ {
		pv[i+int32(ncm)*3] = 0.0
		for j := int32(ncf) - 1; j >= 3; j-- {
			pv[i+int32(ncm)*3] += jc[j] * buf[j+(i+ni*ncm)*ncf]
		}
		pv[i+int32(ncm)*3] *= bma3
	}
}

// interpolation performs the main interpolation routine.
//
// Parameters:
//   - buf: array of Chebyshev coefficients for position.
//   - t: fractional time in the interval covered by coefficients.
//   - intv: length of the whole interval in input time units.
//   - ncfin: number of coefficients per component.
//   - ncmin: number of components per set of coefficients.
//   - nain: number of sets of coefficients in the full array.
//   - ifl: flag indicating which interpolated quantities are required:
//     1 = position only, 2 = position and velocity,
//     3 = position, velocity, and acceleration,
//     4 = position, velocity, acceleration, and jerk.
//   - pv: array to store the interpolated quantities.
//
// Returns 0 on success.
func (js *JPL) interpolation(buf []float64, t, intv float64, ncfin, ncmin, nain, ifl int32, pv []float64) error {
	// Static-like variables
	var nv, nac, njk int32
	var twot float64

	pc := js.Chebyshev.PC
	vc := js.Chebyshev.VC
	ac := js.Chebyshev.AC
	jc := js.Chebyshev.JC
	ncf := ncfin
	ncm := ncmin
	na := nain

	var tc float64
	var ni int32

	// Calculate the sub-interval number and normalized Chebyshev time
	computeSubInterval(t, na, &tc, &ni)

	// Evaluate Chebyshev polynomials if time has changed
	if tc != pc[1] {
		nv = 3
		nac = 4
		njk = 5
		evaluatePolynomials(pc[:], tc, ncf, &twot)
	}

	// Interpolate position
	interpolatePosition(pv, pc[:], buf, ncf, ncm, ni)

	// If velocity is needed
	if ifl > 1 {
		bma := (float64(na) + float64(na)) / intv
		vc[2] = twot + twot
		if nv < ncf {
			for i := nv; i < ncf; i++ {
				vc[i] = twot*vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2]
			}
			nv = ncf
		}
		interpolateVelocity(pv, vc, pc, buf, ncf, ncm, ni, bma)
	}

	// If acceleration is needed
	if ifl > 2 {
		bma2 := (float64(na) + float64(na)) / intv
		ac[3] = pc[1] * 24.0
		if nac < ncf {
			nac = ncf
			for i := nac; i < ncf; i++ {
				ac[i] = twot*ac[i-1] + vc[i-1]*4.0 - ac[i-2]
			}
		}
		interpolateAcceleration(pv, ac, buf, ncf, ncm, ni, bma2)
	}

	// If jerk is needed
	if ifl > 3 {
		bma3 := (float64(na) + float64(na)) / intv
		jc[4] = pc[1] * 192.0
		if njk < ncf {
			njk = ncf
			for i := njk; i < ncf; i++ {
				jc[i] = twot*jc[i-1] + ac[i-1]*6.0 - jc[i-2]
			}
		}
		interpolateJerk(pv, jc, buf, ncf, ncm, ni, bma3)
	}

	return nil
}

// interpolateLibrations interpolates lunar librations based on the input parameters.
func (js *JPL) interpolateLibrations(list []int32, t, intv float64, pv []float64) {
	if list[11] > 0 && js.Constants.IPT[37] > 0 {
		js.interpolation(js.Buf[js.Constants.IPT[36]-1:], t, intv, js.Constants.IPT[37], 3, js.Constants.IPT[38], list[11], pv[60:])
	}
}

// interpolateNutations interpolates nutations in longitude and obliquity based on the input parameters.
func (js *JPL) interpolateNutations(list []int32, t, intv float64, nut []float64) {
	if list[10] > 0 && js.Constants.IPT[34] > 0 {
		js.interpolation(js.Buf[js.Constants.IPT[33]-1:], t, intv, js.Constants.IPT[34], 2, js.Constants.IPT[35], list[10], nut)
	}
}

// interpolateBodies interpolates positions and velocities of celestial bodies based on the input parameters.
func (js *JPL) interpolateBodies(list []int32, doBary bool, aufac, t, intv float64, pv, pvsun []float64) {
	for i := 0; i < 10; i++ {
		if list[i] > 0 {
			js.interpolation(js.Buf[js.Constants.IPT[i*3]-1:], t, intv, js.Constants.IPT[i*3+1], 3, js.Constants.IPT[i*3+2], list[i], pv[i*6:])
			for j := 0; j < 6; j++ {
				if i < 9 && doBary == false {
					pv[j+i*6] = pv[j+i*6]*aufac - pvsun[j]
				} else {
					pv[j+i*6] *= aufac
				}
			}
		}
	}
}

// interpolateSunPosition interpolates the position of the Sun.
func (js *JPL) interpolateSunPosition(t, intv float64, pvsun []float64) {
	js.interpolation(js.Buf[js.Constants.IPT[30]-1:], t, intv, js.Constants.IPT[31], 3, js.Constants.IPT[32], 2, pvsun)
	for i := 0; i < 6; i++ {
		pvsun[i] *= 1.0 / js.Constants.AU
	}
}

// loadJPLRecord loads a specific record from the JPL ephemeris file.
func (js *JPL) loadJPLRecord(nr int32, irecsz int32, ncoeffs int32) error {
	// Seek to the correct position in the file
	offset := int64(nr) * int64(irecsz)
	if _, err := js.JplFile.Seek(offset, 0); err != nil {
		return fmt.Errorf("Read error in JPL eph. at record %d\n", nr)
	}

	// Read the coefficients
	for k := int32(0); k < ncoeffs; k++ {
		if err := binary.Read(js.JplFile, binary.LittleEndian, &js.Buf[k]); err != nil {
			return fmt.Errorf("Read error in JPL eph. at record %d\n", nr)
		}
		if js.DoReorder {
			Reorder(&js.Buf[k]) // Assuming double is 8 bytes
		}
	}

	return nil
}

// checkJplFileIntegrity verifies the integrity of the JPL ephemeris file by checking its length.
func (js *JPL) checkJplFileIntegrity(irecsz int32) error {
	// Seek to the end of the file to determine the file length
	file, err := os.Open(js.JplFilePath)
	if err != nil {
		return fmt.Errorf("Error opening file: %v", err)
	}
	defer file.Close()

	flen, err := file.Seek(0, os.SEEK_END)
	if err != nil {
		return fmt.Errorf("Error determining file length: %v", err)
	}

	nseg := int32((js.Constants.SS[1] - js.Constants.SS[0]) / js.Constants.SS[2])
	nb := int64(0)

	for i := int32(0); i < 13; i++ {
		k := int32(3)
		if i == 11 {
			k = 2
		}
		nb += int64(js.Constants.IPT[i*3+1] * js.Constants.IPT[i*3+2] * k * nseg)
	}
	nb += int64(2 * nseg)
	nb *= 8
	nb += int64(2 * irecsz)

	if flen != nb && flen-nb != int64(irecsz) {
		return fmt.Errorf("JPL ephemeris file %s is mutilated; length = %d instead of %d.", js.JplFilePath, flen, nb)
	}

	return nil
}

// readJplHeader reads and processes the header information from the JPL ephemeris file.
func (js *JPL) readJplHeader(irecsz *int32, ncoeffs *int32) error {
	var ksize int32
	var nrecl int32
	var ch_ttl [252]byte
	var lpt [3]int32
	var nrd int

	ksize, err := js.fsizer() // Assume fsizer is implemented
	if err != nil || ksize == NotAvailable {
		return err
	}

	nrecl = 4
	*irecsz = nrecl * ksize // record size in bytes
	*ncoeffs = ksize / 2    // number of coefficients, doubles

	// Read header
	nrd, err = js.JplFile.Read(ch_ttl[:])
	if err != nil || nrd != 252 {
		return err
	}

	nrd, err = js.JplFile.Read(js.ChCnam[:])
	if err != nil || nrd != 2400 {
		return err
	}

	err = binary.Read(js.JplFile, binary.LittleEndian, &js.Constants.SS)
	if err != nil {
		return err
	}
	if js.DoReorder {
		Reorder(&js.Constants.SS)
	}

	err = binary.Read(js.JplFile, binary.LittleEndian, &js.Constants.Ncon)
	if err != nil {
		return err
	}
	if js.DoReorder {
		Reorder(&js.Constants.Ncon)
	}

	err = binary.Read(js.JplFile, binary.LittleEndian, &js.Constants.AU)
	if err != nil {
		return err
	}
	if js.DoReorder {
		Reorder(&js.Constants.AU)
	}

	err = binary.Read(js.JplFile, binary.LittleEndian, &js.Constants.Emrat)
	if err != nil {
		return err
	}
	if js.DoReorder {
		Reorder(&js.Constants.Emrat)
	}

	var ehIpt [36]int32
	err = binary.Read(js.JplFile, binary.LittleEndian, &ehIpt)
	if err != nil {
		return err
	}

	if js.DoReorder {
		Reorder(&ehIpt)
	}

	copy(js.Constants.IPT[:], ehIpt[:])

	err = binary.Read(js.JplFile, binary.LittleEndian, &js.Constants.Denum)
	if err != nil {
		return err
	}
	if js.DoReorder {
		Reorder(&js.Constants.Denum)
	}

	err = binary.Read(js.JplFile, binary.LittleEndian, &lpt)
	if err != nil {
		return err
	}
	if js.DoReorder {
		Reorder(&lpt)
	}

	// Seek and read constants
	js.JplFile.Seek(int64(1*(*irecsz)), 0)
	err = binary.Read(js.JplFile, binary.LittleEndian, &js.Constants.Cval)
	if err != nil {
		return err
	}
	if js.DoReorder {
		Reorder(&js.Constants.Cval)
	}

	for i := 0; i < 3; i++ {
		js.Constants.IPT[i+36] = lpt[i]
	}

	return nil
}

// state calculates and interpolates the state vectors for celestial bodies based on the input parameters.
/*
	Parameters:
		et (float64): Julian ephemeris epoch at which interpolation is desired.
		interpolationList ([]int32): A 12-element integer array specifying interpolation requests for each celestial body.
			- 0: No interpolation for body i.
			- 1: Position only.
			- 2: Position and velocity.
			Designations for astronomical bodies (index 'i'):
				0: Mercury
				1: Venus
				2: Earth-Moon barycenter (NOT Earth!)
				3: Mars
				4: Jupiter
				5: Saturn
				6: Uranus
				7: Neptune
				8: Pluto
				9: Geocentric Moon
				10: Nutations in longitude and obliquity
				11: Lunar librations (if available)
		doBary (bool): Determines if calculations are barycentric (true) or heliocentric (false).
		pv ([]float64): A 6x11 array that will contain the requested interpolated quantities.
			- The state for each body specified by interpolationList(i) will start at pv(1,i).
			- Order: x, y, z, dx, dy, dz.
		pvsun ([]float64): A 6-element array containing the barycentric position and velocity of the Sun.
		nut ([]float64): A 4-element array that will contain nutations and their rates based on interpolationList(10).
			- Order: d psi (nutation in longitude), d epsilon (nutation in obliquity), d psi dot, d epsilon dot.
*/
func (js *JPL) state(et float64, list []int32, doBary bool, pv, pvsun, nut []float64) error {
	var s, t, etMn, etFr, aufac, intv float64
	var irecsz, ncoeffs int32

	if js.JplFile != nil {
		if err := js.readJplHeader(&irecsz, &ncoeffs); err != nil {
			return err
		}
		if err := js.checkJplFileIntegrity(irecsz); err != nil {
			return err
		}
	}

	if list == nil {
		return nil
	}

	s = et - 0.5
	etMn = math.Floor(s)
	etFr = s - etMn
	etMn += 0.5

	if et < js.Constants.SS[0] || et > js.Constants.SS[1] {
		return fmt.Errorf("jd %f outside JPL eph. range %.2f .. %.2f", et, js.Constants.SS[0], js.Constants.SS[1])
	}

	nr := int32(((etMn - js.Constants.SS[0]) / js.Constants.SS[2]) + 2)
	if etMn == js.Constants.SS[1] {
		nr--
	}

	t = (etMn - ((float64(nr)-2)*js.Constants.SS[2] + js.Constants.SS[0]) + etFr) / js.Constants.SS[2]

	if err := js.loadJPLRecord(nr, irecsz, ncoeffs); err != nil {
		return err
	}

	if js.DoKm {
		intv = js.Constants.SS[2] * 86400.0
		aufac = 1.0
	} else {
		intv = js.Constants.SS[2]
		aufac = 1.0 / js.Constants.AU
	}

	js.interpolateSunPosition(t, intv, pvsun)
	js.interpolateBodies(list, doBary, aufac, t, intv, pv, pvsun)
	js.interpolateNutations(list, t, intv, nut)
	js.interpolateLibrations(list, t, intv, pv)

	return nil
}

func (js *JPL) readConstJpl(ss []float64) error {
	// Call state to initialize the ephemeris and read in the constants
	err := js.state(0.0, nil, false, nil, nil, nil)
	if err != nil {
		return err
	}

	// Copy the first three constants into ss
	copy(ss, js.Constants.SS[:3])

	if !js.DebugMode {
		bname := []string{
			"Mercury", "Venus", "EarthMoonBarycenter", "Mars", "Jupiter", "Saturn",
			"Uranus", "Neptune", "Pluto", "Moon", "SunBary", "Nut", "Libr",
		}

		fmt.Println("JPL EPHEMERIS TEST")
		for i := 0; i < 13; i++ {
			j := i * 3
			k := 3
			if i == 11 {
				k = 2
			}
			nb := js.Constants.IPT[j+1] * js.Constants.IPT[j+2] * int32(k)
			nc := float64(nb*36525) / js.Constants.SS[2] * 8
			fmt.Printf("%s\t%d\tipt[%d]\t%3d %2d %2d,\t",
				bname[i], i, j, js.Constants.IPT[j], js.Constants.IPT[j+1], js.Constants.IPT[j+2])
			fmt.Printf("%3d double, bytes per century = %6d\n", nb, int32(nc))
		}

		fmt.Printf("%16.2f %16.2f %16.2f\n", js.Constants.SS[0], js.Constants.SS[1], js.Constants.SS[2])

		n := int(js.Constants.Ncon)
		if n > len(js.Constants.Cval) {
			n = len(js.Constants.Cval)
		}
		for i := 0; i < n; i++ {
			fmt.Printf("%.6s\t%24.16f\n", strings.TrimSpace(string(js.ChCnam[i*6:(i+1)*6])), js.Constants.Cval[i])
		}
	}
	return nil
}

func (js *JPL) CloseJplFile() {
	if js == nil {
		return // Early return if the JPL struct is nil
	}
	js.JplFile.Close() // Attempt to close the file
	js.JplFile = nil   // Reset the file pointer after closing
}

func (js *JPL) OpenJplFile() ([]float64, error) {
	// Early return if the file is already open
	if js != nil && js.JplFile != nil {
		return []float64{}, nil
	}

	// Attempt to open the JPL file
	var err error
	js.JplFile, err = os.Open(js.JplFilePath)
	if err != nil {
		return nil, err
	}

	// Allocate space for storing constants
	ss := make([]float64, 3)

	// Attempt to read constants from the JPL file
	if err := js.readConstJpl(ss[:]); err != nil {
		js.CloseJplFile() // Clean up by closing file on error
		return nil, err
	}

	// Set initial values for interpolation coefficients
	js.Chebyshev.PC[0] = 1
	js.Chebyshev.PC[1] = 2
	js.Chebyshev.VC[1] = 1
	js.Chebyshev.AC[2] = 4
	js.Chebyshev.JC[3] = 24

	return ss, nil
}

func (js *JPL) GetDenum() int32 {
	if js == nil {
		return -1
	}
	return js.Constants.Denum
}

// Unified Reorder function that dynamically handles different data types and sizes.
func Reorder(data interface{}) {
	switch v := data.(type) {
	case []byte:
		reorderBytes(v)
	case *int32:
		reorderInt32(v)
	case *float64:
		reorderFloat64(v)
	case []int32:
		for i := range v {
			reorderInt32(&v[i])
		}
	case []float64:
		for i := range v {
			reorderFloat64(&v[i])
		}
	default:
		logUnsupportedType(v)
	}
}

// Logs unsupported data types.
func logUnsupportedType(v interface{}) {
	fmt.Printf("Unsupported type provided to Reorder: %T\n", v)
}

// Helper function to reverse byte order for byte slices.
func reorderBytes(data []byte) {
	for i, j := 0, len(data)-1; i < j; i, j = i+1, j-1 {
		data[i], data[j] = data[j], data[i]
	}
}

// Helper function to reverse byte order for int32 pointers.
func reorderInt32(data *int32) {
	bytes := make([]byte, 4)
	binary.LittleEndian.PutUint32(bytes, uint32(*data))
	*data = int32(binary.BigEndian.Uint32(bytes))
}

// Helper function to reverse byte order for float64 pointers.
func reorderFloat64(data *float64) {
	bytes := make([]byte, 8)
	binary.LittleEndian.PutUint64(bytes, math.Float64bits(*data))
	*data = math.Float64frombits(binary.BigEndian.Uint64(bytes))
}
