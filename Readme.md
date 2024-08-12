# JPL Ephemeris Reader for Go

This library provides an implementation to read and interpret the Jet Propulsion Laboratory (JPL) planetary ephemerides binary data files in Go. It supports calculations of positions, velocities, and interactions of celestial bodies within the solar system.

## Features

- Read and process JPL ephemeris files.
- Calculate the position and velocity of planets with respect to each other.
- Handle specific celestial interactions like Earth-Moon.
- Supports different coordinate frames and units (kilometers or astronomical units).
- Debug mode for detailed processing logs.

## Installation

To use this library in your project, import it as follows:

```go
import "github.com/mshafiee/jpl"
```

Ensure you have Go installed on your system. You can install the library directly using:

```bash
go get github.com/mshafiee/jpl
```

## Usage

### Setting Up

First, you need a JPL ephemeris file, which can be obtained from the [NASA JPL website](https://ssd.jpl.nasa.gov/?ephemerides). The library takes the file path as an argument to initialize the JPL struct:

```go
filepath := "path/to/your/linux_m13000p17000.441"
js := jpl.NewJPL(filepath)
```

### Opening the Ephemeris File

Before performing any calculations, open the ephemeris file:

```go
ss, err := js.OpenJplFile()
if err != nil {
    fmt.Println("Error opening JPL file:", err)
    return
}
fmt.Println("Ephemeris opened successfully.")
```

### Ephemeris Lookups

You can calculate the position and velocity of a planet (e.g., Earth) with respect to another (e.g., the Sun) at a specific Julian Date:

```go
et := 2451545.0 // Example Julian Date
ntarg := jpl.CelestialBodyEarth
ncent := jpl.CelestialBodySun

rrd, err := js.EphemerisLookup(et, ntarg, ncent)
if err != nil {
    fmt.Println("Error calculating positions and velocities:", err)
    return
}

fmt.Println("Position and Velocity of Earth relative to Sun:")
fmt.Printf("X: %f AU, Y: %f AU, Z: %f AU\n", rrd[0], rrd[1], rrd[2])
fmt.Printf("dX: %f AU/day, dY: %f AU/day, dZ: %f AU/day\n", rrd[3], rrd[4], rrd[5])
```

### Closing the File

After you are done, close the ephemeris file to free up resources:

```go
js.CloseJplFile()
```

## Contributing

Contributions to the library are welcome. Please feel free to fork the repository, make changes, and submit pull requests. You can also open issues if you find bugs or have feature suggestions.

## License

This project is licensed under the GNU General Public License Version 3 (GPLv3) - see the [LICENSE](LICENSE) file for details. You can also view the full text of the license [here](https://www.gnu.org/licenses/gpl-3.0.en.html).
