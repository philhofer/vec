package stat

import(
	"math"
)

type LinearFit struct {
	Slope float64
	Intercept float64
}

func (lf *LinearFit) Eval(x float64) float64 {
	return lf.Slope*x + lf.Intercept
}

func (lf *LinearFit) SumSquares(xs []float64, ys []float64) float64 {
	if len(xs) != len(ys) {
		return math.NaN()
	}
	ss := 0.0
	for i:=0; i<len(xs); i++ {
		lfval := lf.Eval(xs[i])
		dif := ys[i]-lfval
		ss += dif*dif
	}
	return ss
}

//Perform a simple linear regression of the form y = a + b*x
//'xs' and 'ys' must be of the same length; otherwise nil is returned
func LinFit(xs []float64, ys []float64) *LinearFit {
	if len(xs) != len(ys) {
		return nil
	}

	b := Covariance(xs, ys)/Variance(xs)
	a := Mean(ys) - b*Mean(xs)
	return &LinearFit{Slope: b, Intercept: a}
}

//Performs a Deming regression when 'delta' is known
func DemingReg(xs []float64, ys []float64, delta float64) *LinearFit {
	if len(xs) != len(ys) {
		return nil
	}
	N := len(xs)
	n := float64(N-1)
	sxx := 0.0
	sxy := 0.0
	syy := 0.0
	xbar := Mean(xs)
	ybar := Mean(ys)
	for i:=0; i<N; i++ {
		sxx += (xs[i]-xbar)*(xs[i]-xbar)
		sxy += (xs[i]-xbar)*(ys[i]-ybar)
		syy += (ys[i]-ybar)*(ys[i]-ybar)
	}
	sxx, sxy, syy = sxx/n, sxy/n, syy/n
	b := (syy - delta*sxx + math.Sqrt((syy-delta*sxx)*(syy-delta*sxx) + 4.0*delta*sxy*sxy))/(2.0*sxy)
	a := ybar - b*xbar
	return &LinearFit{Slope: b, Intercept: a}
}
