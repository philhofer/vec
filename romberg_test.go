package vec

import "testing"
import "math"

func myFunc(x float64) float64 {
	return x * x
}

func myInfFunc(x float64) float64 {
	return math.Exp(-5.0 * x)
}

func TestTrapezoids(t *testing.T) {
	out := trap(myFunc, 0.0, 3.0, 256)
	if out-9.0 > 1E-4 {
		t.Error("Trapezoids didn't integrate accurately enough.")
		t.Error("Expected: 9.0")
		t.Error("Got:", out)
	}
}

func TestRomberg(t *testing.T) {
	out, conv := Integral(myFunc, 0.0, 3.0)
	if !conv {
		t.Error("Romberg Integration failed to converge.")
		t.Error("Expected: 9.0")
		t.Error("Got:", out)
	}
}

func TestInfRomberg(t *testing.T) {
	out, conv := Integral(myInfFunc, 0, math.Inf(1))
	if !conv {
		t.Error("Romberg Integration didn't converge.")
		t.Error("Testing exp(-5x) from 0 to +Inf.")
		t.Error("Expecting: 0.2")
		t.Error("Got", out)
	}
}
