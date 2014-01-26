package vec

import "testing"

func myFunc(x float64) float64 {
	return x*x
}

func TestTrapezoids(t *testing.T) {
	out := trap(myFunc, 0.0, 3.0, 500)
	if out - 9.0 > 1E-4 {
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
