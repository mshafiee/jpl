[![Go Reference](https://pkg.go.dev/badge/github.com/mshafiee/jpl.svg)](https://pkg.go.dev/github.com/mshafiee/jpl)
![Coverage](https://img.shields.io/badge/Coverage-48.1%25-yellow)
[![Go Report Card](https://goreportcard.com/badge/github.com/mshafiee/jpl)](https://goreportcard.com/report/github.com/mshafiee/jpl)
# JPL Ephemeris Reader for Go

This Go library allows you to read and interpret the Jet Propulsion Laboratory (JPL) planetary ephemerides binary data files. It supports calculations of positions and velocities of celestial bodies within the solar system.

## Features

- Open and process JPL ephemeris files.
- Calculate the position and velocity of planets relative to other celestial bodies.
- Supports precise calculations using JPL ephemeris data.
- Provides results in astronomical units (AU) and AU/day.
- Supports different coordinate frames and units (kilometers or astronomical units).

## Installation

To use this library in your project, import it as follows:

```go
import "github.com/mshafiee/jpl"
```

Ensure Go is installed on your system. Install the library using:

```bash
go get github.com/mshafiee/jpl
```

## Usage

### Running the Example Program

The provided `main.go` demonstrates how to use the library to calculate the position and velocity of Earth relative to the Sun.

To run the program, you need to pass the path to a JPL ephemeris file as a command-line argument:

```bash
go run main.go path/to/your/ephemeris_file
```

### Example Code

Here’s a quick example of how to use the library:

```go
package main

import (
    "fmt"
    "github.com/mshafiee/jpl"
    "os"
)

func main() {
    if len(os.Args) < 2 {
        fmt.Printf("Usage: %s <ephemeris_file>\n", os.Args[0])
        return
    }

    filepath := os.Args[1]

    jplFileReader, err := os.Open(filepath)
    if err != nil {
        panic(err)
    }
    defer jplFileReader.Close()

    js, ss, err := jpl.NewJPL(jplFileReader)
    if err != nil {
        fmt.Printf("Error opening JPL file: %s\n", err)
        return
    }

    fmt.Println("Ephemeris opened successfully.")
    fmt.Printf("Start Epoch: %.2f, End Epoch: %.2f, Granule Size: %.2f\n", ss[0], ss[1], ss[2])

    denum := js.GetDenum()
    fmt.Printf("DE Number: %d\n", denum)

    et := 2451545.0
    ntarg := jpl.Earth
    ncent := jpl.Sun

    rrd, err := js.EphemerisLookup(et, ntarg, ncent)
    if err != nil {
        fmt.Printf("Error calculating positions and velocities: %s\n", err)
        return
    }

    fmt.Println("Position and Velocity of Earth relative to Sun:")
    fmt.Printf("X: %f AU, Y: %f AU, Z: %f AU\n", rrd[0], rrd[1], rrd[2])
    fmt.Printf("dX: %f AU/day, dY: %f AU/day, dZ: %f AU/day\n", rrd[3], rrd[4], rrd[5])
}
```

### Understanding the Output

```plaintext
Ephemeris opened successfully.
Start Epoch: -3027215.50, End Epoch: 7930192.50, Granule Size: 32.00
DE Number: 441
Position and Velocity of Earth relative to Sun:
X: -0.177135 AU, Y: 0.887429 AU, Z: 0.384743 AU
dX: -0.017208 AU/day, dY: -0.002898 AU/day, dZ: -0.001256 AU/day
```

The example program prints the following details:

- **Start Epoch, End Epoch, Granule Size:** Information about the ephemeris file.
- **DE Number:** The DE number of the ephemeris data set.
- **Position and Velocity:** The position (in AU) and velocity (in AU/day) of Earth relative to the Sun on the specified Julian date.

## Contributing

Contributions are welcome! Feel free to fork the repository, make changes, and submit a pull request. You can also open issues for bugs or feature requests.

## License

This project is licensed under the GNU General Public License Version 3 (GPLv3). See the [LICENSE](LICENSE) file for details or view the full license text [here](https://www.gnu.org/licenses/gpl-3.0.en.html).
