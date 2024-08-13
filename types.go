package jpl

import (
	"encoding/binary"
	"io"
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

// FileReader defines an interface for file operations used in JPL.
type FileReader interface {
	io.Reader
	io.Seeker
	io.Closer
}

// headerValidator defines an interface for validating and reading file headers.
type headerValidator interface {
	readAndValidateHeader() error
	readJplHeader() (irecsz int32, ncoeffs int32, err error)
	checkJplFileIntegrity(irecsz int32) error
}

// constantsReader defines an interface for reading constants and parameters.
type constantsReader interface {
	readConstantsAndParameters(au *float64, emrat *float64, lpt []int32, numde *int32, ncon *int32) error
	readConstJpl() ([]float64, error)
}

// celestialBodyHandler defines an interface for handling specific celestial body interactions.
type celestialBodyHandler interface {
	handleNutation(et float64, pv, pvsun []float64) ([]float64, error)
	handleLibration(et float64, pv, pvsun []float64) ([]float64, error)
	handleEarthMoonInteraction(list []int32, pv []float64)
}

// interpolationHelper defines an interface for functions related to interpolation calculations.
type interpolationHelper interface {
	setupListForState(ntarg, ncent CelestialBody) []int32
	adjustPositionsForSunEMBBary(ntarg, ncent CelestialBody, pv, pvsun []float64)
	computeSubInterval(t float64, na int32) (tc float64, ni int32)
	evaluatePolynomials(pc []float64, tc float64, ncf int32) (twot float64)
	interpolatePosition(pv []float64, pc []float64, buf []float64, ncf, ncm, ni int32)
	interpolateVelocity(pv []float64, vc, pc [18]float64, buf []float64, ncf, ncm, ni int32, bma float64)
	interpolateAcceleration(pv []float64, ac [18]float64, buf []float64, ncf, ncm, ni int32, bma2 float64)
	interpolateJerk(pv []float64, jc [18]float64, buf []float64, ncf, ncm, ni int32, bma3 float64)
}

// ephemerisHandler defines a combined interface for handling state calculations,
type ephemerisHandler interface {
	state(et float64, list []int32, doBary bool, pv, pvsun, rrd []float64) error
	EphemerisLookup(et float64, ntarg, ncent CelestialBody) ([]float64, error)
	interpolation(buf []float64, t, intv float64, ncfin, ncmin, nain, ifl int32, pv []float64) error
	interpolateBodies(list []int32, doBary bool, aufac, t, intv float64, pv, pvsun []float64)
	interpolateSunPosition(t, intv float64, pvsun []float64)
	interpolateLibrations(list []int32, t, intv float64, pv []float64)
	interpolateNutations(list []int32, t, intv float64, nut []float64)
}

// jplFileManager defines an interface for file management operations.
type jplFileManager interface {
	openJplFile() ([]float64, error)
	loadJPLRecord(nr int32, irecsz int32, ncoeffs int32) error
	fsizer() (int32, error)
}

// JPL is the main struct that implements all the interfaces.
type JPL struct {
	JplFile   FileReader
	Endian    binary.ByteOrder
	Constants struct {
		Cval  [400]float64
		SS    [3]float64
		AU    float64
		Emrat float64
		Denum int32
		Ncon  int32
		IPT   [39]int32
	}
	ChCnam    [6 * 400]byte
	PV        [78]float64
	PVSun     [6]float64
	Buf       [1500]float64
	Chebyshev struct {
		PC [18]float64
		VC [18]float64
		AC [18]float64
		JC [18]float64
	}
	DoKm                 bool
	DebugMode            bool
	ephemerisHandler     ephemerisHandler
	headerValidator      headerValidator
	constantsReader      constantsReader
	fileManager          jplFileManager
	celestialBodyHandler celestialBodyHandler
	interpolationHelper  interpolationHelper
}

// Separate structs for each interface implementation

type jplFileManagerImpl struct {
	jpl *JPL
}

func (fm *jplFileManagerImpl) openJplFile() ([]float64, error) {
	return fm.jpl.openJplFile()
}

func (fm *jplFileManagerImpl) loadJPLRecord(nr int32, irecsz int32, ncoeffs int32) error {
	return fm.jpl.loadJPLRecord(nr, irecsz, ncoeffs)
}

func (fm *jplFileManagerImpl) fsizer() (int32, error) {
	return fm.jpl.fsizer()
}

type headerValidatorImpl struct {
	jpl *JPL
}

func (hv *headerValidatorImpl) readAndValidateHeader() error {
	return hv.jpl.readAndValidateHeader()
}

func (hv *headerValidatorImpl) readJplHeader() (int32, int32, error) {
	return hv.jpl.readJplHeader()
}

func (hv *headerValidatorImpl) checkJplFileIntegrity(irecsz int32) error {
	return hv.jpl.checkJplFileIntegrity(irecsz)
}

type constantsReaderImpl struct {
	jpl *JPL
}

func (cr *constantsReaderImpl) readConstantsAndParameters(au *float64, emrat *float64, lpt []int32, numde *int32, ncon *int32) error {
	return cr.jpl.readConstantsAndParameters(au, emrat, lpt, numde, ncon)
}

func (cr *constantsReaderImpl) readConstJpl() ([]float64, error) {
	return cr.jpl.readConstJpl()
}

type celestialBodyHandlerImpl struct {
	jpl *JPL
}

func (cb *celestialBodyHandlerImpl) handleNutation(et float64, pv, pvsun []float64) ([]float64, error) {
	return cb.jpl.handleNutation(et, pv, pvsun)
}

func (cb *celestialBodyHandlerImpl) handleLibration(et float64, pv, pvsun []float64) ([]float64, error) {
	return cb.jpl.handleLibration(et, pv, pvsun)
}

func (cb *celestialBodyHandlerImpl) handleEarthMoonInteraction(list []int32, pv []float64) {
	cb.jpl.handleEarthMoonInteraction(list, pv)
}

type interpolationHelperImpl struct {
	jpl *JPL
}

func (ih *interpolationHelperImpl) setupListForState(ntarg, ncent CelestialBody) []int32 {
	return ih.jpl.setupListForState(ntarg, ncent)
}

func (ih *interpolationHelperImpl) adjustPositionsForSunEMBBary(ntarg, ncent CelestialBody, pv, pvsun []float64) {
	ih.jpl.adjustPositionsForSunEMBBary(ntarg, ncent, pv, pvsun)
}

func (ih *interpolationHelperImpl) computeSubInterval(t float64, na int32) (tc float64, ni int32) {
	return ih.jpl.computeSubInterval(t, na)
}

func (ih *interpolationHelperImpl) evaluatePolynomials(pc []float64, tc float64, ncf int32) (twot float64) {
	return ih.jpl.evaluatePolynomials(pc, tc, ncf)
}

func (ih *interpolationHelperImpl) interpolatePosition(pv []float64, pc []float64, buf []float64, ncf, ncm, ni int32) {
	ih.jpl.interpolatePosition(pv, pc, buf, ncf, ncm, ni)
}

func (ih *interpolationHelperImpl) interpolateVelocity(pv []float64, vc, pc [18]float64, buf []float64, ncf, ncm, ni int32, bma float64) {
	ih.jpl.interpolateVelocity(pv, vc, pc, buf, ncf, ncm, ni, bma)
}

func (ih *interpolationHelperImpl) interpolateAcceleration(pv []float64, ac [18]float64, buf []float64, ncf, ncm, ni int32, bma2 float64) {
	ih.jpl.interpolateAcceleration(pv, ac, buf, ncf, ncm, ni, bma2)
}

func (ih *interpolationHelperImpl) interpolateJerk(pv []float64, jc [18]float64, buf []float64, ncf, ncm, ni int32, bma3 float64) {
	ih.jpl.interpolateJerk(pv, jc, buf, ncf, ncm, ni, bma3)
}

type ephemerisHandlerImpl struct {
	jpl *JPL
}

func (eh *ephemerisHandlerImpl) state(et float64, list []int32, doBary bool, pv, pvsun, rrd []float64) error {
	return eh.jpl.state(et, list, doBary, pv, pvsun, rrd)
}

func (eh *ephemerisHandlerImpl) EphemerisLookup(et float64, ntarg, ncent CelestialBody) ([]float64, error) {
	return eh.jpl.EphemerisLookup(et, ntarg, ncent)
}

func (eh *ephemerisHandlerImpl) interpolation(buf []float64, t, intv float64, ncfin, ncmin, nain, ifl int32, pv []float64) error {
	return eh.jpl.interpolation(buf, t, intv, ncfin, ncmin, nain, ifl, pv)
}

func (eh *ephemerisHandlerImpl) interpolateBodies(list []int32, doBary bool, aufac, t, intv float64, pv, pvsun []float64) {
	eh.jpl.interpolateBodies(list, doBary, aufac, t, intv, pv, pvsun)
}

func (eh *ephemerisHandlerImpl) interpolateSunPosition(t, intv float64, pvsun []float64) {
	eh.jpl.interpolateSunPosition(t, intv, pvsun)
}

func (eh *ephemerisHandlerImpl) interpolateLibrations(list []int32, t, intv float64, pv []float64) {
	eh.jpl.interpolateLibrations(list, t, intv, pv)
}

func (eh *ephemerisHandlerImpl) interpolateNutations(list []int32, t, intv float64, nut []float64) {
	eh.jpl.interpolateNutations(list, t, intv, nut)
}
