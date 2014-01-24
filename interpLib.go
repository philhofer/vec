package vec 


func Arange(start float64, stop float64, N int) []float64 {
	f := make([]float64, N)
	h := (stop - start) / float64(N)
	for i := range f {
		f[i] = start + (float64(i)*h)
	}
	return f
}

//Type + Container for Cubic Spline Interpolation
type CubicSplineInterpolation struct {
	data *BiVariateData
	coeffs []float64
}

//'Natural' Cubic Spline
func CubicSpline(d BiVariateData) CubicSplineInterpolation {
	spline := CubicSplineInterpolation{data: &d}
	N := len(d.Ys)
	fs := make([]float64, N)
	fs[0] = 3*(d.Ys[1] - d.Ys[0])
	for i:=1; i<(N-1); i++ {
		fs[i] = 3*(d.Ys[i+1] - d.Ys[i-1])
	}
	fs[N-1] = 3*(d.Ys[N-1] - d.Ys[N-2])

	getCubicSplinCoeffs(fs)
	spline.coeffs = fs

	return spline
}

//Find interpolated f(x)
func (s CubicSplineInterpolation) F(x float64) float64{
	i, i1 := s.data.findXBounds(x)
	t := (x - s.data.Xs[i])/(s.data.Xs[i1] - s.data.Xs[i])
	yi, yi1 := s.data.Ys[i], s.data.Ys[i1]
	Di, Di1 := s.coeffs[i], s.coeffs[i+1]
	a := yi
	b := Di
	c := 3*(yi1 - yi) - 2*Di - Di1
	d := 2*(yi - yi1) + Di + Di1

	return a + b*t + c*t*t + d*t*t*t
}

//Find interpolated df(x)/dx
func (s CubicSplineInterpolation) DF(x float64) float64 {
	i, i1 := s.data.findXBounds(x)
	fin := s.data.Xs[i1]
	start := s.data.Xs[i]
	t := (x - start)/(fin - start)
	yi, yi1 := s.data.Ys[i], s.data.Ys[i1]
	Di, Di1 := s.coeffs[i], s.coeffs[i+1]
	b := Di
	c := 3*(yi1 - yi) - 2*Di - Di1
	d := 2*(yi - yi1) + Di + Di1

	return (b + 2*c*t + 3*d*t*t)/(fin-start)
}

func (s CubicSplineInterpolation) DDF(x float64) float64 {
	i, i1 := s.data.findXBounds(x)
	fin := s.data.Xs[i1]
	start := s.data.Xs[i]
	t := (x - start)/(fin - start)
	yi, yi1 := s.data.Ys[i], s.data.Ys[i1]
	Di, Di1 := s.coeffs[i], s.coeffs[i+1]
	c := 3*(yi1 - yi) - 2*Di - Di1
	d := 2*(yi - yi1) + Di + Di1

	return (2*c + 6*d*t)/((fin-start)*(fin-start))
}

func makeConstVec(a float64, N int) []float64{
	out := make([]float64, N)
	for i := range out {
		x := a
		out[i] = x
	}
	return out
}


func getCubicSplinCoeffs(x []float64) {
	l := len(x)

	//subdiagonal
	a := makeConstVec(1.0, l)

	//main diag
	b := makeConstVec(4.0, l)
	b[0] = 2.0
	b[l-1] = 2.0

	//superdiagonal
	c := makeConstVec(1.0, l)

	c[0] = (c[0] / b[0])
	x[0] = (x[0] / b[0])

	//gaussian elimination
	for i := 1; i<l; i++ {
		m := 1.0 / (b[i] - (a[i] * c[i-1]))
		c[i] = c[i] * m
		x[i] = (x[i] - (a[i] * x[i-1]))*m
	}

	//backsubstitution
	for i := l-2; i >= 0; i-- {
		x[i] = x[i] - c[i]*x[i+1]
	}

}
