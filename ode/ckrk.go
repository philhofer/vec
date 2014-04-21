package ode

/* Cash-Karp Explicit Runge-Kutta method
(This Runge-Kutta method uses a fourth- and
fifth-order approximation to generate error approximations
for each step. Solve() creates a wrapper for each step of the
function, setting 'h' at each step to keep the estimated error below
the specified error tolerance.)
*/

import (
	"math"
)

const (
	//Butcher Tableau
	C2 = 0.2
	C3 = 0.3
	C4 = 0.6
	C5 = 1.0
	C6 = 0.875
	A21 = 0.2
	A31 = 0.075
	A32 = 0.225
	A41 = 0.3
	A42 = -0.9
	A43 = 1.2
	A51 = -0.203703703703703703703
	A52 = 2.5
	A53 = -2.592592592592592592592
	A54 = 1.2962962962962962962962
	A61 = 0.0294958043981481481481
	A62 = 0.3141796875
	A63 = 0.0415943287037037037037
	A64 = 0.4003454137731481481481
	A65 = 0.061767578125
	B11 = 0.0978835978835978835978
	B13 = 0.4025764895330112721417
	B14 = 0.2104377104377104377104
	B16 = 0.2891022021456804065499
	B21 = 0.1021773726851851851851
	B23 = 0.3839079034391534391534
	B24 = 0.2445927372685185185185
	B25 = 0.0193219866071428571485
	B26 = 0.25
)

type ODEFunc func(xs []float64, t float64) (nxs []float64)

func CKRKStep(df ODEFunc, x0 []float64, t float64, h float64) (nxs []float64, err float64) {
	currX := make([]float64, len(x0))
	k1 := df(x0, t)
	for i,_ := range k1 {
		k1[i] *= h
		currX[i] = x0[i] + A21*k1[i]
	}
	k2 := df(currX, t + C2*h)
	for i,_ := range k2 {
		k2[i] *= h
		currX[i] = x0[i] + A31*k1[i] + A32*k2[i]
	}
	k3 := df(currX, t + C3*h)
	for i,_ := range k3 {
		k3[i] *= h
		currX[i] = x0[i] + A41*k1[i] + A42*k2[i] + A43*k3[i]
	}
	k4 := df(currX, t + C4*h)
	for i,_ := range k4 {
		k4[i] *= h
		currX[i] = x0[i] + A51*k1[i] + A52*k2[i] + A53*k3[i] + A54*k4[i]
	}
	k5 := df(currX, t + C5*h)
	for i,_ := range k5 {
		k5[i] *= h
		currX[i] = x0[i] + A61*k1[i] + A62*k2[i] + A63*k3[i] + A64*k4[i] + A65*k5[i]
	}
	k6 := df(currX, t + C6*h)
	for i,_ := range k6 {
		k6[i] *= h
	}
	est1 := make([]float64, len(currX))
	nxs = make([]float64, len(currX))
	err = 0.0
	for i,_ := range est1 {
		est1[i] = x0[i] + B11*k1[i] + B13*k3[i] + B14*k4[i] + B16*k6[i]
		nxs[i] = x0[i] + B21*k1[i] + B23*k3[i] + B24*k4[i] + B25*k5[i] + B26*k6[i]
		err += math.Abs(est1[i] - nxs[i]) 
	}
	err = err/float64(len(currX))
	return
}

func Solve(fn ODEFunc, x0 []float64, t0 float64, tfinal float64, globalerr float64) (out []float64, t float64) {
	out = make([]float64, len(x0))
	copy(x0, out)
	if globalerr < 0.0 { panic("Error needs to be positive.") }
	if tfinal < t0 { panic("tfinal must be less than t0.") }
	t = t0
	h := 0.01

	//main loop - go until t>t0
	for t < tfinal {
		thisf, thiserr := CKRKStep(fn, out, t, h)
		h = 0.99*h*math.Pow(math.Abs(globalerr/thiserr), 0.2)
		if math.Abs(thiserr) < globalerr {
			t += h
			out = thisf
		}
	}
	return 
}
