package stat

import(
	"math"
	"math/rand"
)

type LinearFit struct {
	lambda float64
	xs []float64
	ys []float64
	Slope float64
	Intercept float64
}

func (lf *LinearFit) Eval(x float64) float64 {
	return lf.Slope*x + lf.Intercept
}

/* Linear Fit Bootstraps
Creates N Linear Fits by bootstrap resampling of underlying data,
returns an array of length N of slopes and intercepts, respectively.
*/
func (fit *LinearFit) Bootstraps(N int) (slopes []float64, intercepts []float64) {
	slopes, intercepts = make([]float64, N), make([]float64, N)
	// make fit
	for i:=0; i<N; i++ {
		this_x, this_y := resample(fit.xs, fit.ys)
		this_fit := DemingReg(this_x, this_y, fit.lambda)
		slopes[i], intercepts[i] = this_fit.Slope, this_fit.Intercept
	}
	return
}

func resample(xs []float64, ys []float64) ([]float64, []float64) {
	nsamp := len(xs)
	indexes := make([]int, nsamp)
	for i, _ := range indexes {
		indexes[i] = rand.Intn(nsamp)
	}
	return elements(xs, indexes), elements(ys, indexes)
}

func elements(arr []float64, indexes []int) []float64 {
	out := make([]float64, len(indexes))
	for i, x := range indexes {
		out[i] = arr[x]
	}
	return out
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

/*Orthogonal Least-Squares Linear Fit
Performs an orthogonal least-squares fit of a straight line
through the data points 'xs' and 'ys.' The arrangement of 'xs'
and 'ys' is irrelevant, as orthogonal least-squares fits are robust
to inversion. (In other words, LinFit(xs, ys) is equivalent to Linfit(ys, xs)).
*/
func LinFit(xs []float64, ys []float64) *LinearFit {
	return DemingReg(xs, ys, 1.0)
}

/*Deming regression
Performs a least-squares straight-line fit when the ratio of the
variances of errors between 'x' and 'y' ('delta,' in this case) is known.
Note that delta can be calculated, given error slices 'y_errors' and 'x_errors,' by
Variance(y_errors)/Variance(x_errors)
*/
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
	return &LinearFit{lambda: delta, xs: xs, ys: ys, Slope: b, Intercept: a}
}
