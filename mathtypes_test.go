package vec

import "testing"
import "math"

//BiMathops
func add(x float64, y float64) float64 {
	return x + y
}

func mult(x float64, y float64) float64 {
	return x * y
}

//Mathops
func square(x float64) float64 {
	return x*x
}

func divbytwo(x float64) float64 {
	return x/2
}

//Condenser
func sum(x []float64) float64 {
	out := 0.0
	for _, xi := range x {
		out += xi
	}
	return out
}

var arrOne []float64 = []float64{1, 2, 3, 4}

func compareFloat(x float64, y float64) bool {
	if math.Abs(x - y) > 1E-15 {
		return true
	} else {
		return false
	}
}

func TestFold(t *testing.T) {
	if Fold(add, arrOne) != 10 {
		t.Error("Fold failed on 'add'")
		t.Error("Expected 10, got", Fold(add, arrOne))
	} else if Fold(mult, arrOne) != 24 {
		t.Error("Fold failed on 'mult'")
		t.Error("Expected 24, got", Fold(mult, arrOne))
	} else {
		return
	}
}

func TestFtMathop(t *testing.T) {
	sq_future := MakeFtMathop(square)
	divbytwo_future := MakeFtMathop(divbytwo)

	c4 := sq_future(2.0)
	c1 := divbytwo_future(2.0)

	a := <- c4
	b := <- c1
	if compareFloat(a, 4.0) {
		t.Error("MakeFtMathop failed on 'square'")
	} else if compareFloat(b, 1.0) {
		t.Error("MakeFtMathop failed on 'divbytwo'")
	} else {
		return
	}
}


func TestFtBiMathop(t *testing.T) {
	add_future := MakeFtBiMathop(add)
	mult_future := MakeFtBiMathop(mult)

	c1 := add_future(2.0, -5.6)
	c2 := mult_future(12.0, 12.0)

	a := <- c1
	b := <- c2

	if compareFloat(a, -3.6) {
		t.Error("MakeFtBiMathop failed on 'add' -- Yielded", a)
	} else if compareFloat(b, 144.0) {
		t.Error("MakeFtBiMathop failed on 'mult' -- Yielded", b)
	} else {
		return
	}
}
