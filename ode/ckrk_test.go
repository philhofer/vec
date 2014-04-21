package ode

import (
	"testing"
	"math"
)

func globalODE(xy []float64, t float64) (dxy []float64) {
	dxy = make([]float64, 2)
	x := xy[0]
	y := xy[1]
	dxy[0] = 2.0 - y
	dxy[1] = 4.0 - x
	return
}

func TestCkrkStep(t *testing.T) {
	initial := []float64{0.0, 0.0}
	nxs, maxerr := CKRKStep(globalODE, initial, 0.0, 0.01)
	xActual := 0.0198003
	yActual := 0.0399007
	t.Log("Got x, y = ", xActual, yActual)
	t.Log("Got Max Err = ", maxerr)
	if math.Abs(xActual-nxs[0]) > 10E-6 {
		t.Error("ckrkStep() Failed. Expected x =", xActual, "got", nxs[0])
	}
	if math.Abs(yActual-nxs[1]) > 10E-6 {
		t.Error("ckrkStep() Failed. Expected y =", yActual, "got", nxs[1])
	}
	return
}
