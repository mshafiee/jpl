package jpl

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"strings"
	"unsafe"
)

// NewJPL initializes a new JPL structure with default values.
// It takes a FileReader interface as input, which represents the JPL ephemeris file to be read.
// The function sets up the JPL structure with necessary handlers, determines the system's endianness,
// and attempts to open and read the ephemeris file. It returns a pointer to the JPL structure,
// an array of constants read from the file, and any error encountered during the process.
func NewJPL(fileReader FileReader) (*JPL, []float64, error) {
	jpl := &JPL{
		JplFile: fileReader,
		Endian:  getEndian(),
	}

	jpl.initHandlers()

	ss, err := jpl.fileManager.openJplFile()
	if err != nil {
		return nil, nil, fmt.Errorf("Error opening JPL file: %w", err)
	}

	return jpl, ss, nil
}

// initHandlers initializes all the necessary handlers within the JPL structure.
// These handlers manage various aspects of ephemeris data processing, such as
// reading constants, handling celestial bodies, and performing interpolation operations.
func (jpl *JPL) initHandlers() {
	jpl.ephemerisHandler = &ephemerisHandlerImpl{jpl}
	jpl.headerValidator = &headerValidatorImpl{jpl}
	jpl.constantsReader = &constantsReaderImpl{jpl}
	jpl.fileManager = &jplFileManagerImpl{jpl}
	jpl.celestialBodyHandler = &celestialBodyHandlerImpl{jpl}
	jpl.interpolationHelper = &interpolationHelperImpl{jpl}
}

// getEndian determines the byte order (endianness) of the system.
// It returns binary.BigEndian if the system is big-endian, and binary.LittleEndian if it is little-endian.
func getEndian() binary.ByteOrder {
	if isBigEndian() {
		return binary.BigEndian
	}
	return binary.LittleEndian
}

// isBigEndian checks if the system's byte order is big-endian.
// It does so by setting an integer value and inspecting the byte order of its memory representation.
func isBigEndian() bool {
	var i int32 = 0x01020304
	var b [4]byte
	*(*int32)(unsafe.Pointer(&b[0])) = i
	return b[0] == 0x01
}

// --- celestialBodyHandler Interface Implementation ---
// handleNutation calculates the nutation effects at a given ephemeris time (et).
// It checks if nutation data is available in the ephemeris file and computes the nutation state if so.
// The function returns the computed nutation values or an error if the computation fails.
func (jpl *JPL) handleNutation(et float64, pv, pvsun []float64) ([]float64, error) {
	rrd := make([]float64, 6)

	// Check if nutation data is available in the ephemeris
	if jpl.Constants.IPT[34] <= 0 {
		return nil, fmt.Errorf("nutation data not available in the JPL ephemeris file")
	}

	list := make([]int32, 12)
	list[10] = 2

	// Calculate state for nutation
	if err := jpl.ephemerisHandler.state(et, list, false, pv, pvsun, rrd); err != nil {
		return nil, fmt.Errorf("failed to calculate nutation state: %w", err)
	}

	return rrd, nil
}

// handleLibration calculates the libration effects at a given ephemeris time (et).
// Similar to handleNutation, it checks if libration data is available in the ephemeris file
// and computes the libration state if so. The function returns the computed libration values
// or an error if the computation fails.
func (jpl *JPL) handleLibration(et float64, pv, pvsun []float64) ([]float64, error) {
	rrd := make([]float64, 6)

	// Check if libration data is available in the ephemeris
	if jpl.Constants.IPT[37] <= 0 {
		return nil, fmt.Errorf("libration data not available in the JPL ephemeris file")
	}

	list := make([]int32, 12)
	list[11] = 2

	// Calculate state for libration
	if err := jpl.ephemerisHandler.state(et, list, false, pv, pvsun, rrd); err != nil {
		return nil, fmt.Errorf("failed to calculate libration state: %w", err)
	}

	// Ensure the pv array is correctly populated
	copy(rrd, pv[60:66])

	return rrd, nil
}

// handleEarthMoonInteraction adjusts the positions of Earth and Moon in the pv array based on
// the specified list of celestial bodies. This method applies the Earth-Moon ratio (EMRAT)
// to correct the relative positions of Earth and Moon, ensuring accurate interaction between these bodies.
func (jpl *JPL) handleEarthMoonInteraction(list []int32, pv []float64) {
	if len(list) == 0 {
		return
	}

	if list[Earth] == 2 {
		for i := 0; i < 6; i++ {
			pv[i+6*int(Earth)] -= pv[i+6*int(Moon)] / (jpl.Constants.Emrat + 1.0)
		}
	}
	if list[Moon] == 2 {
		for i := 0; i < 6; i++ {
			pv[i+6*int(Moon)] += pv[i+6*int(Earth)]
		}
	}
}

// --- interpolationHelper Interface Implementation ---

// setupListForState creates a list of celestial bodies for which state vectors (position and velocity)
// need to be computed. The list is based on the target and center bodies provided. This list is used
// for interpolation of celestial body states within the ephemeris.
func (jpl *JPL) setupListForState(ntarg, ncent CelestialBody) []int32 {
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

// adjustPositionsForSunEMBBary adjusts the positions in the pv array for cases where
// the target or center celestial body is the Sun, the Solar System Barycenter, or the Earth-Moon Barycenter.
// This adjustment ensures the positions are correctly aligned for these special cases.
func (jpl *JPL) adjustPositionsForSunEMBBary(ntarg, ncent CelestialBody, pv, pvsun []float64) {
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

// computeSubInterval computes the sub-interval within the ephemeris record for a given time (t).
// This method returns the normalized Chebyshev time (tc) and the sub-interval index (ni).
// These values are essential for evaluating the Chebyshev polynomials used in ephemeris interpolation.
func (jpl *JPL) computeSubInterval(t float64, na int32) (tc float64, ni int32) {
	var dt1 float64

	if t >= 0 {
		dt1 = math.Floor(t)
	} else {
		dt1 = -math.Floor(-t)
	}

	temp := float64(na) * t
	ni = int32(temp - dt1)
	tc = (math.Mod(temp, 1.0)+dt1)*2.0 - 1.0

	return
}

// evaluatePolynomials computes the Chebyshev polynomials based on the provided normalized Chebyshev time (tc)
// and the number of coefficients (ncf). The method fills the polynomial coefficient array (pc) and
// returns the value of two times the Chebyshev time (twot), which is used in further interpolations.
func (jpl *JPL) evaluatePolynomials(pc []float64, tc float64, ncf int32) (twot float64) {
	if ncf < 2 {
		return tc + tc
	}

	twot = tc + tc
	pc[1] = tc
	for i := int32(2); i < ncf; i++ {
		pc[i] = twot*pc[i-1] - pc[i-2]
	}
	return
}

// interpolatePosition calculates the position of a celestial body using the provided Chebyshev polynomial coefficients (pc),
// the interpolation buffer (buf), and the sub-interval parameters. The results are stored in the position vector (pv).
func (jpl *JPL) interpolatePosition(pv []float64, pc []float64, buf []float64, ncf, ncm, ni int32) {
	for i := int32(0); i < ncm; i++ {
		pv[i] = 0.0
		for j := ncf - 1; j >= 0; j-- {
			pv[i] += pc[j] * buf[j+(i+ni*ncm)*ncf]
		}
	}
}

// interpolateVelocity calculates the velocity of a celestial body based on the velocity coefficients (vc),
// polynomial coefficients (pc), and the interpolation buffer (buf). The results are stored in the velocity vector (pv),
// and the method accounts for the scaling factor (bma) related to the time interval.
func (jpl *JPL) interpolateVelocity(pv []float64, vc, pc [18]float64, buf []float64, ncf, ncm, ni int32, bma float64) {
	for i := int32(0); i < ncm; i++ {
		pv[i+ncm] = 0.0
		for j := int32(ncf) - 1; j >= 1; j-- {
			pv[i+ncm] += vc[j] * buf[j+(i+ni*ncm)*ncf]
		}
		pv[i+ncm] *= bma
	}
}

// interpolateAcceleration calculates the acceleration of a celestial body using the acceleration coefficients (ac),
// polynomial coefficients (pc), and the interpolation buffer (buf). The results are stored in the acceleration vector (pv),
// and the method uses a scaling factor (bma2) to adjust the acceleration values based on the time interval.
func (jpl *JPL) interpolateAcceleration(pv []float64, ac [18]float64, buf []float64, ncf, ncm, ni int32, bma2 float64) {
	for i := int32(0); i < ncm; i++ {
		pv[i+int32(ncm)*2] = 0.0
		for j := int32(ncf) - 1; j >= 2; j-- {
			pv[i+int32(ncm)*2] += ac[j] * buf[j+(i+ni*ncm)*ncf]
		}
		pv[i+int32(ncm)*2] *= bma2
	}
}

// interpolateJerk calculates the jerk (the derivative of acceleration) of a celestial body using the jerk coefficients (jc),
// polynomial coefficients (pc), and the interpolation buffer (buf). The results are stored in the jerk vector (pv),
// and the method uses a scaling factor (bma3) to adjust the jerk values based on the time interval.
func (jpl *JPL) interpolateJerk(pv []float64, jc [18]float64, buf []float64, ncf, ncm, ni int32, bma3 float64) {
	for i := int32(0); i < ncm; i++ {
		pv[i+int32(ncm)*3] = 0.0
		for j := int32(ncf) - 1; j >= 3; j-- {
			pv[i+int32(ncm)*3] += jc[j] * buf[j+(i+ni*ncm)*ncf]
		}
		pv[i+int32(ncm)*3] *= bma3
	}
}

// --- ephemerisHandler Interface Implementation ---

// state computes the state vectors (position, velocity, etc.) of celestial bodies for a given ephemeris time (et).
// It performs various operations, including validation of the JPL file header, checking file integrity,
// interpolation of positions, velocities, and handling of special celestial body interactions. The method returns
// any error encountered during these operations.
func (jpl *JPL) state(et float64, list []int32, doBary bool, pv, pvsun, nut []float64) error {
	var s, t, etMn, etFr, aufac, intv float64
	var irecsz, ncoeffs int32
	var err error

	if jpl.JplFile != nil {
		// Validate header and retrieve necessary parameters
		irecsz, ncoeffs, err = jpl.headerValidator.readJplHeader()
		if err != nil {
			return fmt.Errorf("failed to read JPL header: %w", err)
		}

		// Check integrity of the JPL file
		if err := jpl.headerValidator.checkJplFileIntegrity(irecsz); err != nil {
			return fmt.Errorf("JPL file integrity check failed: %w", err)
		}
	}

	if list == nil {
		return nil
	}

	s = et - 0.5
	etMn = math.Floor(s)
	etFr = s - etMn
	etMn += 0.5

	// Ensure the provided Julian Date is within the valid range
	if et < jpl.Constants.SS[0] || et > jpl.Constants.SS[1] {
		return fmt.Errorf("Julian Date %f is outside the valid JPL ephemeris range: %.2f to %.2f", et, jpl.Constants.SS[0], jpl.Constants.SS[1])
	}

	nr := int32(((etMn - jpl.Constants.SS[0]) / jpl.Constants.SS[2]) + 2)
	if etMn == jpl.Constants.SS[1] {
		nr--
	}

	t = (etMn - ((float64(nr)-2)*jpl.Constants.SS[2] + jpl.Constants.SS[0]) + etFr) / jpl.Constants.SS[2]

	// Load the JPL record and handle any errors
	if err := jpl.fileManager.loadJPLRecord(nr, irecsz, ncoeffs); err != nil {
		return fmt.Errorf("failed to load JPL record for record number %d: %w", nr, err)
	}

	if jpl.DoKm {
		intv = jpl.Constants.SS[2] * 86400.0
		aufac = 1.0
	} else {
		intv = jpl.Constants.SS[2]
		aufac = 1.0 / jpl.Constants.AU
	}

	jpl.ephemerisHandler.interpolateSunPosition(t, intv, pvsun)
	jpl.ephemerisHandler.interpolateBodies(list, doBary, aufac, t, intv, pv, pvsun)
	jpl.ephemerisHandler.interpolateNutations(list, t, intv, nut)
	jpl.ephemerisHandler.interpolateLibrations(list, t, intv, pv)

	return nil
}

// EphemerisLookup returns the state vectors (position and velocity) of a celestial body relative to another
// celestial body for a given ephemeris time (et). The function handles special cases like nutations and librations,
// and computes relative positions while adjusting for Earth-Moon interactions. It returns the computed state vectors
// or an error if the calculation fails.
func (jpl *JPL) EphemerisLookup(et float64, ntarg, ncent CelestialBody) ([]float64, error) {
	if ntarg == ncent {
		// Return an empty result when the target and center are the same
		return []float64{}, nil
	}

	switch ntarg {
	case Nutations:
		// Handle nutation lookup
		rrd, err := jpl.celestialBodyHandler.handleNutation(et, jpl.PV[:], jpl.PVSun[:])
		if err != nil {
			return nil, fmt.Errorf("failed to lookup nutations: %w", err)
		}
		return rrd, nil
	case Librations:
		// Handle libration lookup
		rrd, err := jpl.celestialBodyHandler.handleLibration(et, jpl.PV[:], jpl.PVSun[:])
		if err != nil {
			return nil, fmt.Errorf("failed to lookup librations: %w", err)
		}
		return rrd, nil
	}

	rrd := make([]float64, 6)
	list := jpl.interpolationHelper.setupListForState(ntarg, ncent)

	// Calculate the state of the celestial body
	if err := jpl.ephemerisHandler.state(et, list, true, jpl.PV[:], jpl.PVSun[:], rrd); err != nil {
		return nil, fmt.Errorf("failed to calculate state for target %d and center %d: %w", ntarg, ncent, err)
	}

	// Adjust positions for the Sun and Earth-Moon barycenter
	jpl.interpolationHelper.adjustPositionsForSunEMBBary(ntarg, ncent, jpl.PV[:], jpl.PVSun[:])

	// Handle specific Earth-Moon interactions
	if !(ntarg == Earth && ncent == Moon) && !(ntarg == Moon && ncent == Earth) {
		jpl.celestialBodyHandler.handleEarthMoonInteraction(list, jpl.PV[:])
	}

	// Calculate the relative positions
	for i := int32(0); i < 6; i++ {
		rrd[i] = jpl.PV[i+int32(ntarg)*6] - jpl.PV[i+int32(ncent)*6]
	}

	return rrd, nil
}

// interpolation handles the interpolation of state vectors (position, velocity, acceleration, and jerk)
// for celestial bodies at a given time (t). The method evaluates Chebyshev polynomials based on the time interval
// and uses these to interpolate the state vectors, storing the results in the provided position vector (pv).
// The ifl parameter controls the depth of interpolation (position, velocity, acceleration, jerk).
func (jpl *JPL) interpolation(buf []float64, t, intv float64, ncfin, ncmin, nain, ifl int32, pv []float64) error {
	// Static-like variables
	var nv, nac, njk int32
	var twot float64

	pc := jpl.Chebyshev.PC
	vc := jpl.Chebyshev.VC
	ac := jpl.Chebyshev.AC
	jc := jpl.Chebyshev.JC
	ncf := ncfin
	ncm := ncmin
	na := nain

	// Calculate the sub-interval number and normalized Chebyshev time
	tc, ni := jpl.interpolationHelper.computeSubInterval(t, na)

	// Evaluate Chebyshev polynomials if time has changed
	if tc != pc[1] {
		nv = 3
		nac = 4
		njk = 5
		twot = jpl.interpolationHelper.evaluatePolynomials(pc[:], tc, ncf)
	}

	// Interpolate position
	jpl.interpolationHelper.interpolatePosition(pv, pc[:], buf, ncf, ncm, ni)

	// If velocity is needed
	if ifl > 1 {
		bma := (float64(na) + float64(na)) / intv
		vc[2] = twot + twot
		if nv < ncf {
			for i := nv; i < ncf; i++ {
				vc[i] = twot*vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2]
			}
		}
		jpl.interpolationHelper.interpolateVelocity(pv, vc, pc, buf, ncf, ncm, ni, bma)
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
		jpl.interpolationHelper.interpolateAcceleration(pv, ac, buf, ncf, ncm, ni, bma2)
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
		jpl.interpolationHelper.interpolateJerk(pv, jc, buf, ncf, ncm, ni, bma3)
	}

	return nil
}

// interpolateBodies computes the state vectors for a set of celestial bodies, adjusting the positions
// based on whether the barycentric (doBary) flag is set. The method handles the interpolation of positions
// for each body in the list and applies any necessary adjustments for the Sun or other barycenters.
func (jpl *JPL) interpolateBodies(list []int32, doBary bool, aufac, t, intv float64, pv, pvsun []float64) {
	for i := 0; i < 10; i++ {
		if list[i] > 0 {
			jpl.ephemerisHandler.interpolation(jpl.Buf[jpl.Constants.IPT[i*3]-1:], t, intv, jpl.Constants.IPT[i*3+1], 3, jpl.Constants.IPT[i*3+2], list[i], pv[i*6:])
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

// interpolateSunPosition interpolates the position of the Sun at a given time (t),
// storing the results in the provided position vector (pvsun). The method uses the
// ephemeris data for the Sun to calculate its position relative to the Solar System barycenter.
func (jpl *JPL) interpolateSunPosition(t, intv float64, pvsun []float64) {
	jpl.ephemerisHandler.interpolation(jpl.Buf[jpl.Constants.IPT[30]-1:], t, intv, jpl.Constants.IPT[31], 3, jpl.Constants.IPT[32], 2, pvsun)
	for i := 0; i < 6; i++ {
		pvsun[i] *= 1.0 / jpl.Constants.AU
	}
}

// interpolateLibrations handles the interpolation of libration values if they are requested
// (as indicated by the list parameter). The method uses the provided ephemeris time (t) and
// time interval (intv) to compute the librations and stores the results in the position vector (pv).
func (jpl *JPL) interpolateLibrations(list []int32, t, intv float64, pv []float64) {
	if list[11] > 0 && jpl.Constants.IPT[37] > 0 {
		jpl.ephemerisHandler.interpolation(jpl.Buf[jpl.Constants.IPT[36]-1:], t, intv, jpl.Constants.IPT[37], 3, jpl.Constants.IPT[38], list[11], pv[60:])
	}
}

// interpolateNutations handles the interpolation of nutation values if they are requested
// (as indicated by the list parameter). The method uses the provided ephemeris time (t) and
// time interval (intv) to compute the nutations and stores the results in the nutation vector (nut).
func (jpl *JPL) interpolateNutations(list []int32, t, intv float64, nut []float64) {
	if list[10] > 0 && jpl.Constants.IPT[34] > 0 {
		jpl.ephemerisHandler.interpolation(jpl.Buf[jpl.Constants.IPT[33]-1:], t, intv, jpl.Constants.IPT[34], 2, jpl.Constants.IPT[35], list[10], nut)
	}
}

// GetDenum returns the ephemeris number (Denum) from the JPL constants.
// The Denum value identifies the specific version of the JPL ephemeris being used.
func (jpl *JPL) GetDenum() int32 {
	return jpl.Constants.Denum
}

// --- jplFileManager Interface Implementation ---

// openJplFile opens the JPL ephemeris file and reads the necessary constants.
// It ensures the JPL instance and file are properly initialized before attempting to read constants
// from the file. The method returns the array of constants (ss) and any error encountered during the process.
func (jpl *JPL) openJplFile() ([]float64, error) {
	// Check if the JPL instance and file are properly initialized
	if jpl == nil || jpl.JplFile == nil {
		return nil, fmt.Errorf("failed to open JPL file: JPL instance is uninitialized or JPL file is nil")
	}

	// Attempt to read constants from the JPL file
	ss, err := jpl.constantsReader.readConstJpl()
	if err != nil {
		return nil, fmt.Errorf("failed to read constants from JPL file: %w", err)
	}

	// Initialize Chebyshev interpolation coefficients
	jpl.Chebyshev.PC[0] = 1
	jpl.Chebyshev.PC[1] = 2
	jpl.Chebyshev.VC[1] = 1
	jpl.Chebyshev.AC[2] = 4
	jpl.Chebyshev.JC[3] = 24

	return ss, nil
}

// loadJPLRecord reads the specified record from the JPL ephemeris file into the buffer.
// It seeks to the correct position in the file based on the record number (nr) and record size (irecsz),
// then reads the coefficients into the buffer. The method returns any error encountered during reading.
func (jpl *JPL) loadJPLRecord(nr int32, irecsz int32, ncoeffs int32) error {
	// Seek to the correct position in the file
	offset := int64(nr) * int64(irecsz)
	if _, err := jpl.JplFile.Seek(offset, 0); err != nil {
		return fmt.Errorf("Read error in JPL eph. at record %d\n", nr)
	}

	// Read the coefficients
	for k := int32(0); k < ncoeffs; k++ {
		if err := binary.Read(jpl.JplFile, jpl.Endian, &jpl.Buf[k]); err != nil {
			return fmt.Errorf("Read error in JPL eph. at record %d\n", nr)
		}
	}

	return nil
}

// fsizer reads and validates the JPL file header, and calculates the ksize (record size).
// It ensures the file header is correctly parsed and that constants and parameters are
// within expected ranges. The method returns the calculated ksize and any error encountered.
func (jpl *JPL) fsizer() (int32, error) {
	var au, emrat float64
	var lpt [3]int32
	var numde, ncon int32

	// Attempt to read and validate the JPL file header
	if err := jpl.headerValidator.readAndValidateHeader(); err != nil {
		return NotAvailable, fmt.Errorf("error reading and validating JPL file header: %w", err)
	}

	// Attempt to read constants and parameters from the JPL file
	if err := jpl.constantsReader.readConstantsAndParameters(&au, &emrat, lpt[:], &numde, &ncon); err != nil {
		return NotAvailable, fmt.Errorf("error reading constants and parameters: %w", err)
	}

	// Calculate the ksize and check if it is within the expected range
	ksize := calculateKsize(jpl.Constants.IPT[:])
	if ksize < 1000 || ksize > 5000 {
		return NotAvailable, fmt.Errorf("invalid ksize value in JPL ephemeris file: %d (expected between 1000 and 5000)", ksize)
	}

	return ksize, nil
}

// readJplHeader reads and processes the header information from the JPL ephemeris file.
// It extracts important metadata such as the record size (irecsz), number of coefficients (ncoeffs),
// and other critical parameters needed for subsequent file operations. The method returns
// these values along with any error encountered during reading.
func (jpl *JPL) readJplHeader() (irecsz int32, ncoeffs int32, err error) {
	var ksize int32
	var nrecl int32 = 4
	var ch_ttl [252]byte
	var lpt [3]int32
	var nrd int

	// Retrieve file size (ksize)
	ksize, err = jpl.fileManager.fsizer()
	if err != nil {
		return 0, 0, fmt.Errorf("failed to get file size: %w", err)
	}
	if ksize == NotAvailable {
		return 0, 0, fmt.Errorf("file size not available")
	}

	irecsz = nrecl * ksize // record size in bytes
	ncoeffs = ksize / 2    // number of coefficients, doubles

	// Read header title
	nrd, err = jpl.JplFile.Read(ch_ttl[:])
	if err != nil {
		return 0, 0, fmt.Errorf("failed to read header title: %w", err)
	}
	if nrd != len(ch_ttl) {
		return 0, 0, fmt.Errorf("incomplete header title read: expected %d bytes, got %d", len(ch_ttl), nrd)
	}

	// Read constant names (ChCnam)
	nrd, err = jpl.JplFile.Read(jpl.ChCnam[:])
	if err != nil {
		return 0, 0, fmt.Errorf("failed to read constant names: %w", err)
	}
	if nrd != len(jpl.ChCnam) {
		return 0, 0, fmt.Errorf("incomplete constant names read: expected %d bytes, got %d", len(jpl.ChCnam), nrd)
	}

	// Read SS constants
	if err = binary.Read(jpl.JplFile, jpl.Endian, &jpl.Constants.SS); err != nil {
		return 0, 0, fmt.Errorf("failed to read SS constants: %w", err)
	}

	// Read Ncon
	if err = binary.Read(jpl.JplFile, jpl.Endian, &jpl.Constants.Ncon); err != nil {
		return 0, 0, fmt.Errorf("failed to read Ncon: %w", err)
	}

	// Read AU
	if err = binary.Read(jpl.JplFile, jpl.Endian, &jpl.Constants.AU); err != nil {
		return 0, 0, fmt.Errorf("failed to read AU: %w", err)
	}

	// Read Emrat
	if err = binary.Read(jpl.JplFile, jpl.Endian, &jpl.Constants.Emrat); err != nil {
		return 0, 0, fmt.Errorf("failed to read Emrat: %w", err)
	}

	// Read IPT constants (36 values)
	var ehIpt [36]int32
	if err = binary.Read(jpl.JplFile, jpl.Endian, &ehIpt); err != nil {
		return 0, 0, fmt.Errorf("failed to read IPT constants: %w", err)
	}
	copy(jpl.Constants.IPT[:], ehIpt[:])

	// Read Denum
	if err = binary.Read(jpl.JplFile, jpl.Endian, &jpl.Constants.Denum); err != nil {
		return 0, 0, fmt.Errorf("failed to read Denum: %w", err)
	}

	// Read LPT
	if err = binary.Read(jpl.JplFile, jpl.Endian, &lpt); err != nil {
		return 0, 0, fmt.Errorf("failed to read LPT: %w", err)
	}

	// Seek to constants position and read Cval
	if _, err := jpl.JplFile.Seek(int64(irecsz), 0); err != nil {
		return 0, 0, fmt.Errorf("failed to seek to constants position: %w", err)
	}
	if err = binary.Read(jpl.JplFile, jpl.Endian, &jpl.Constants.Cval); err != nil {
		return 0, 0, fmt.Errorf("failed to read constants (Cval): %w", err)
	}

	// Update IPT with LPT values
	for i := 0; i < 3; i++ {
		jpl.Constants.IPT[i+36] = lpt[i]
	}

	return irecsz, ncoeffs, nil
}

// checkJplFileIntegrity verifies the integrity of the JPL ephemeris file by checking its length.
// It compares the calculated length of the file with its actual size to ensure the file has not been corrupted.
// The method returns an error if the file size does not match the expected length.
func (jpl *JPL) checkJplFileIntegrity(irecsz int32) error {
	// Seek to the end of the file to determine the file length
	fileSize, err := jpl.JplFile.Seek(0, io.SeekEnd)
	if err != nil {
		return fmt.Errorf("Error determining file length: %v", err)
	}

	nseg := int32((jpl.Constants.SS[1] - jpl.Constants.SS[0]) / jpl.Constants.SS[2])
	nb := int64(0)

	for i := int32(0); i < 13; i++ {
		k := int32(3)
		if i == 11 {
			k = 2
		}
		nb += int64(jpl.Constants.IPT[i*3+1] * jpl.Constants.IPT[i*3+2] * k * nseg)
	}
	nb += int64(2 * nseg)
	nb *= 8
	nb += int64(2 * irecsz)

	if fileSize != nb && fileSize-nb != int64(irecsz) {
		return fmt.Errorf("JPL ephemeris file is mutilated; length = %d instead of %d.", fileSize, nb)
	}

	return nil
}

// readAndValidateHeader reads and validates the header information of the JPL file.
// It ensures that the file's header is correctly formatted and that the constants (SS array)
// are within the expected range. The method returns an error if any validation checks fail.
func (jpl *JPL) readAndValidateHeader() error {
	var ttl [6 * 14 * 3]byte
	var ss [3]float64

	_, err := jpl.JplFile.Seek(0, io.SeekStart)
	if err != nil {
		return fmt.Errorf("failed to reset file pointer to start for reading: %v", err)
	}

	// Read the TTL array
	nrd, err := jpl.JplFile.Read(ttl[:])
	if err != nil || nrd != 252 {
		return fmt.Errorf("failed to read TTL array, %v", err)
	}

	// Read the CNAM array
	nrd, err = jpl.JplFile.Read(jpl.ChCnam[:])
	if err != nil || nrd != 6*400 {
		return errors.New("failed to read CNAM array")
	}

	// Read the SS array using the determined endianness
	for i := 0; i < 3; i++ {
		if err = binary.Read(jpl.JplFile, jpl.Endian, &ss[i]); err != nil {
			return errors.New("failed to read SS array")
		}
	}

	// Store the values in SS
	for i := 0; i < 3; i++ {
		jpl.Constants.SS[i] = ss[i]
	}

	// Validate the SS values
	if jpl.Constants.SS[0] < -5583942 || jpl.Constants.SS[1] > 9025909 || jpl.Constants.SS[2] < 1 || jpl.Constants.SS[2] > 200 {
		return fmt.Errorf("alleged ephemeris file has an invalid format")
	}

	return nil
}

// readConstJpl reads constants and initializes the JPL ephemeris.
// It performs initial setup by calling state and checking if debug mode is enabled, printing detailed
// information about the ephemeris if necessary. The method returns the constants array and any error encountered.
func (jpl *JPL) readConstJpl() ([]float64, error) {
	// Initialize the ephemeris and read in the constants
	if err := jpl.ephemerisHandler.state(0.0, nil, false, nil, nil, nil); err != nil {
		return nil, fmt.Errorf("failed to initialize ephemeris and read constants: %w", err)
	}

	if jpl.DebugMode {
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
			nb := jpl.Constants.IPT[j+1] * jpl.Constants.IPT[j+2] * int32(k)
			nc := float64(nb*36525) / jpl.Constants.SS[2] * 8
			fmt.Printf("%s\t%d\tipt[%d]\t%3d %2d %2d,\t",
				bname[i], i, j, jpl.Constants.IPT[j], jpl.Constants.IPT[j+1], jpl.Constants.IPT[j+2])
			fmt.Printf("%3d double, bytes per century = %6d\n", nb, int32(nc))
		}

		fmt.Printf("%16.2f %16.2f %16.2f\n", jpl.Constants.SS[0], jpl.Constants.SS[1], jpl.Constants.SS[2])

		n := int(jpl.Constants.Ncon)
		if n > len(jpl.Constants.Cval) {
			n = len(jpl.Constants.Cval)
		}
		for i := 0; i < n; i++ {
			fmt.Printf("%.6s\t%24.16f\n", strings.TrimSpace(string(jpl.ChCnam[i*6:(i+1)*6])), jpl.Constants.Cval[i])
		}
	}

	return jpl.Constants.SS[:3], nil
}

// readConstantsAndParameters reads the constants and parameters from the JPL file.
// It reads essential ephemeris constants such as the astronomical unit (AU), Earth-Moon mass ratio (EMRAT),
// and planetary interpolation parameters (IPT). The method returns an error if any of these values cannot be read.
func (jpl *JPL) readConstantsAndParameters(au *float64, emrat *float64, lpt []int32, numde *int32, ncon *int32) error {
	var buf4 [4]byte // Buffer for 4-byte reads
	var buf8 [8]byte // Buffer for 8-byte reads

	// Read ncon (4 bytes)
	if _, err := jpl.JplFile.Read(buf4[:]); err != nil {
		return fmt.Errorf("failed to read 'ncon' value: %w", err)
	}
	*ncon = int32(jpl.Endian.Uint32(buf4[:]))

	// Read au (8 bytes)
	if _, err := jpl.JplFile.Read(buf8[:]); err != nil {
		return fmt.Errorf("failed to read 'au' value: %w", err)
	}
	*au = math.Float64frombits(jpl.Endian.Uint64(buf8[:]))

	// Read emrat (8 bytes)
	if _, err := jpl.JplFile.Read(buf8[:]); err != nil {
		return fmt.Errorf("failed to read 'emrat' value: %w", err)
	}
	*emrat = math.Float64frombits(jpl.Endian.Uint64(buf8[:]))

	// Read IPT (36 values, 4 bytes each)
	for i := 0; i < 36; i++ {
		if _, err := jpl.JplFile.Read(buf4[:]); err != nil {
			return fmt.Errorf("failed to read 'IPT[%d]' value: %w", i, err)
		}
		jpl.Constants.IPT[i] = int32(jpl.Endian.Uint32(buf4[:]))
	}

	// Read numde (4 bytes)
	if _, err := jpl.JplFile.Read(buf4[:]); err != nil {
		return fmt.Errorf("failed to read 'numde' value: %w", err)
	}
	*numde = int32(jpl.Endian.Uint32(buf4[:]))

	// Read lpt (3 values, 4 bytes each)
	for i := 0; i < 3; i++ {
		if _, err := jpl.JplFile.Read(buf4[:]); err != nil {
			return fmt.Errorf("failed to read 'lpt[%d]' value: %w", i, err)
		}
		lpt[i] = int32(jpl.Endian.Uint32(buf4[:]))
		jpl.Constants.IPT[i+36] = lpt[i]
	}

	// Reset file pointer to the beginning
	if _, err := jpl.JplFile.Seek(0, io.SeekStart); err != nil {
		return fmt.Errorf("failed to reset file pointer to the beginning: %w", err)
	}

	return nil
}

// calculateKsize calculates the record size (ksize) based on the planetary interpolation parameters (IPT).
// The ksize value is used to determine the size of the ephemeris records within the JPL file.
// The method returns the calculated ksize or panics if the IPT slice is not properly sized.
func calculateKsize(ipt []int32) int32 {
	if len(ipt) < 39 {
		panic("ipt slice must have at least 39 elements")
	}

	kmx, khi := int32(-1), int32(0)
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

	// The de102 files have an incorrect ksize due to 424 empty bytes per record.
	// This issue was manually corrected.
	if ksize == 1546 {
		ksize = 1652
	}

	return ksize
}
