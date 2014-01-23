package vec

/*
Discreet Data Operations




*/

type BiVariateData struct {
	x []float64
	y []float64
}

type Interpolation interface {
	//Interpolation value
	f(x float64)
	
	//First derivative
	df(x float64)

	//Third derivative
	ddf(x float64)
}


/* 
func CubicSpline(data BiVariateData) Interpolation {}

Returns a CubicSplineInterpolation interface on 'data'



*/

/*
func DiscreetConvolve(dat []float64, conv []float64) []float64 {}

Returns a discreet convolution of 'conv' on 'dat'


*/
