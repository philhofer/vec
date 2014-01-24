package vec

//Mathop - Univariate math function
type Mathop func(float64) float64

//BiMathop - Bivariate math function
type BiMathop func(float64, float64) float64

//Condenser - turns an array into a single element
type Condenser func([]float64) float64

//Linop - takes an array, yields an array of the same size
type Linop func([]float64) []float64

//Vecop - operates on two arrays of the same size to yield a new array
type Vecop func([]float64, []float64) []float64

//FtMathop - a future Mathop (returned on a channel)
type FtMathop func(float64) chan float64

//FtBiMathop - a future BiMathop (returned on a channel)
type FtBiMathop func(float64, float64) chan float64

//FtCondenser - a future Condenser (returned on a channel)
type FtCondenser func([]float64) chan float64

//FtLinop - a future Linop (returned on a channel)
type FtLinop func([]float64) chan []float64

//MakeFtMathop - makes a future from a Mathop
func MakeFtMathop(f Mathop) FtMathop {
	g := func(x float64) chan float64 {
		c := make(chan float64, 1)
		go func() {
			c <- f(x)
			return
		}()
		return c
	}
	return g
}

//Makes a future from a BiMathop; returns a channel
func MakeFtBiMathop(f BiMathop) FtBiMathop {
	g := func(x float64, y float64) chan float64 {
		c := make(chan float64, 1)
		go func() {
			c <- f(x, y)
			return
		}()
		return c
	}
	return g
}

//Makes a future from a Linop; returns a channel
func MakeFtLinop(f Linop) FtLinop {
	g := func(x []float64) chan []float64 {
		c := make(chan []float64, 1)
		go func() {
			c <- f(x)
			return
		}()
		return c
	}
	return g
}

//Canonical foldl - folds from vec[0] to vec[n]
func Fold(f BiMathop, vec []float64) float64 {
	if len(vec) <= 1 {
		return vec[0]
	} else {
		val := vec[0]
		nvec := vec[1:]
		for _, x := range nvec {
			val = f(val, x)
		}
		return val
	}
}


//Romberg Integration
//func romberg_integrate(f Mathop, start float64, end float64) float64 {}

//Functional Integral
//func Integral(f Mathop) BiMathop {}

//Convolution
//func Convolution(f Mathop, g Mathop) Mathop {}

//Numerical 1st Derivative
//func NDeriv(f Mathop) float64 {}

type VectorField struct {
	rule Linop
}

func (v VectorField) Val(pts []float64) []float64{
	out := make([]float64, len(pts))
	out = v.rule(pts)
	return out
}
