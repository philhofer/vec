package vec

import "testing"
import "math"

func TestFindRoot1(t *testing.T) {
	out, conv := FindRoot(math.J0, 0.0, 4.0)
	if !conv {
		t.Error("FindRoot failed to converge.")
	} else {
		if math.Abs(out-2.4048255576957727) > 1E-14 {
			t.Error("FindRoot didn't compute the right value.")
			t.Error("Expected:", 2.4048255576957727)
			t.Error("Got:", out)
			t.Error("Evaluates to:", math.J0(out))
		}
	}
}
