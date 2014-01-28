package vec

import "testing"
import "math"
import "math/rand"
import "reflect"

//For computing relative error
func relativeError(out float64, expct float64) float64 {
	return math.Abs((out - expct) / expct)
}

//Comparing interpolation(s) with expected data
func areSimilar(a float64, b float64) bool {
	const MaxAllowedError float64 = 1E-12
	if relativeError(a, b) > MaxAllowedError {
		return false
	} else {
		return true
	}
}

/*
Arange - yields N equally-spaced elements from a to b

EDGE CASES:
- N = 0 -> return []float64{a}
*/
func TestArange(t *testing.T) {
	zero := Arange(0, 1, 0)
	if !reflect.DeepEqual(zero, []float64{}) {
		t.Error("Arange(x, y, 0) doesn't return []float64")
	}
	a := -rand.Float64() * 1000
	b := rand.Float64() * 1000
	n := rand.Intn(100000) + 1
	h := (b - a) / float64(n)
	many := Arange(a, b, n)
	if areSimilar(many[0], a) {
		//pass
	} else {
		t.Error("Arange doesn't start with the right value. Got", many[0], "Should be:", a)
	}

	if areSimilar(many[n-1], b-h) {
		//pass
	} else {
		t.Error("Arange doesn't end with the right value. Got:", many[n-1], "Should be:", b-h)
	}

}

/*
Test initialization of CubicSpline object

EDGE CASES:
- check that evaluation of empty struct literal or nil yields
another empty struct literal
*/
func TestCubicSpline(t *testing.T) {
	bvd := BiVariateData{}
	spl := CubicSpline(bvd)
	if !reflect.DeepEqual(spl, CubicSplineInterpolation{data: &bvd}) {
		t.Error("CubicSpline(BiVariateData{}) didn't return CubicSplineInterpolation{data: &bvd}")
	}
}

//verify data operations
func TestXBounds(t *testing.T) {
	xs := Arange(0, 20, 500)
	ys := Arange(0, 20, 500)
	PPmap(math.J1, ys)

	bvd := MakeBiVariateData(xs, ys)
	spl := CubicSpline(*bvd)

	zero, one := spl.data.findXBounds(-1.0)
	if one != 1 || zero != 0 {
		t.Error("findXBounds failed on x < x_first")
	}

	slast, last := spl.data.findXBounds(21.0)
	if slast != 498 || last != 499 {
		t.Error("findXBounds failed on x > x_last")
	}
}

/*
Test evaluation of interpolation at x

EDGE CASES:
- test for linear extrapolation outside data bounds
*/
func TestSplineOuts(t *testing.T) {
	//Data represents Bessel function of 1st kind
	xs := Arange(0, 20, 500)
	ys := Arange(0, 20, 500)
	PPmap(math.J1, ys)

	bvd := MakeBiVariateData(xs, ys)
	spl := CubicSpline(*bvd)

	for i, _ := range xs {
		yi := ys[i]
		xi := xs[i]
		if !areSimilar(yi, spl.F(xi)) {
			t.Error("Spline does not conform to yi = f(xi) rule")
		}
		if !areSimilar(spl.coeffs[i], spl.DF(xi)*(20.0/500.0)) {
			t.Error("Spline does not conform to Di*h = df(xi) rule")
		}
	}

	//Check integration over (0, 20) with closed-form answer
	relError := relativeError(1.0-math.J0(20.0), spl.Integral(0.0, 20.0))
	if relError > 1E-5 {
		t.Error("Spline did not integrate accurately enough. Error:", relError)
	}

}

func BenchmarkSpline(b *testing.B) {
	xs := Arange(0, 10, 1000)
	ys := Arange(0, 10, 1000)

	PPmap(math.Cos, ys)

	b.ResetTimer()
	bvd := MakeBiVariateData(xs, ys)
	for i := 0; i < b.N; i++ {
		spl := CubicSpline(*bvd)
		spl.F(5.46827)
		spl.DF(5.2908)
		spl.DDF(2.2908)
		spl.Integral(0.0, 10.0)
	}

}
