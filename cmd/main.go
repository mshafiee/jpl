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

	js := jpl.NewJPL(filepath)

	ss, err := js.OpenJplFile()
	if err != nil {
		fmt.Printf("Error opening JPL file: %s\n", err)
		return
	}

	fmt.Println("Ephemeris opened successfully.")
	fmt.Printf("Start Epoch: %.2f, End Epoch: %.2f, Granule Size: %.2f\n", ss[0], ss[1], ss[2])

	denum := js.GetDenum()
	fmt.Printf("DE Number: %d\n", denum)

	et := 2451545.0
	ntarg := jpl.Earth // J_EARTH
	ncent := jpl.Sun   // J_SUN

	rrd, err := js.EphemerisLookup(et, ntarg, ncent)
	if err != nil {
		fmt.Printf("Error calculating positions and velocities: %s\n", err)
		js.CloseJplFile()
		return
	}

	fmt.Println("Position and Velocity of Earth relative to Sun:")
	fmt.Printf("X: %f AU, Y: %f AU, Z: %f AU\n", rrd[0], rrd[1], rrd[2])
	fmt.Printf("dX: %f AU/day, dY: %f AU/day, dZ: %f AU/day\n", rrd[3], rrd[4], rrd[5])

}
