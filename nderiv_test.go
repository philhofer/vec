package vec

import (
	"math"
	"testing"
)

func TestDerivBasic(t *testing.T) {
	out, conv := Deriv(math.Cos, math.Pi/2.0)
	if !conv {
		t.Error("Deriv() did not converge on Cos()")
		t.Error("Output:", out)
		t.Error("Expected -1.0")
	} else {
		if math.Abs(out+1.0) > 10E-12 {
			t.Error("Deriv was not accurate enough.")
			t.Error("Expected -1.0")
			t.Error("Got:", out)
		}
	}
}
