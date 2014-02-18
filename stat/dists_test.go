package stat

import(
	"testing"
	"math"
)

func TestStudentT(t *testing.T) {
	dist := StudentTDist(3)
	v1 := dist.CDF(1.8)
	if math.Abs(v1 - 0.91516) > 10E-4 {
		t.Error("Student T distribution CDF did not pass.")
		t.Error("Expected 0.91516.")
		t.Error("Got:", v1)
	}

	v2 := dist.CDF(-.4)
	if math.Abs(v2 - 0.357968) > 10E-3 {
		t.Error("Student T dist CDF did not pass.")
		t.Error("Expected 0.357968")
		t.Error("Got:", v2)
	}
}
