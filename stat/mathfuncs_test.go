package stat

import (
	"testing"
	"math"
)

func TestBeta(t *testing.T) {
	b1 := Beta(0.5, 0.5)
	if math.Abs(b1 - math.Pi) > 10E-15 {
		t.Error("Beta did not achieve desired accuracy.")
		t.Error("Expected pi. Got:", b1)
	}

	b2 := Beta(12.22, 1.384)
	if math.Abs(b2 - 0.02721500351095505) > 10E-11 {
		t.Error("Beta did not achieve desired accuracy.")
		t.Error("Expected 0.02721500351095505")
		t.Error("Got:", b2)
	}
}

func TestIncBeta(t *testing.T) {
	b1 := IncBeta(0.25, 0.5, 1.5)
	if math.Abs(b1 - 0.608997781044229) > 10E-11 {
		t.Error("IncBeta did not achieve desired accuracy.")
		t.Error("Expected 0.608997781044229")
		t.Error("Got:", b1)
	}

	b2 := IncBeta(8, 2, 4)
	if math.Abs(b2 +  79232) > 10E-11 {
		t.Error("IncBeta did not achieve desired accuracy.")
		t.Error("Expected 79232")
		t.Error("Got:", b2)
	}
}
