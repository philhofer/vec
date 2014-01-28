package vec

import "testing"
import "math"
import "math/rand"

type testpair struct {
	in  []float64
	out []float64
}

type triple struct {
	a1  []float64
	a2  []float64
	out []float64
}

func (s testpair) match() bool {
	var mPrec float64 = 1E-15
	if len(s.in) != len(s.out) {
		return false
	} else {
		for i, x := range s.in {
			if math.Abs(x-s.out[i]) > mPrec {
				return false
			}
		}
		return true
	}
}

func (t triple) yields(f BiMathop) bool {
	var mPrec float64 = 1E-15
	if len(t.out) != len(t.a1) {
		return false
	} else {
		for i, x := range t.out {
			if math.Abs(x-f(t.a1[i], t.a2[i])) > mPrec {
				return false
			}
		}
	}
	return true
}

func genPair(length int, f Mathop) testpair {
	in := make([]float64, length)
	out := make([]float64, length)

	for i := 0; i < length; i++ {
		in[i] = rand.Float64() * 100
		out[i] = f(in[i])

	}

	var tpair testpair = testpair{in, out}
	return tpair
}

func genTriple(length int, f BiMathop) triple {
	a1 := make([]float64, length)
	a2 := make([]float64, length)
	out := make([]float64, length)

	for i := 0; i < length; i++ {
		a1[i] = rand.Float64()
		a2[i] = rand.Float64() * 500
		out[i] = f(a1[i], a2[i])

	}

	var genout triple = triple{a1, a2, out}
	return genout

}

func BenchmarkFold(b *testing.B) {
	g := genTriple(10000, add)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Fold(add, g.a1)
	}
}

//Test MPmap with Cos()
func TestMPmapCos(t *testing.T) {
	localarr := genPair(10000, math.Cos)
	MPmap(math.Cos, localarr.in)
	if localarr.match() {
		return
	} else {
		t.Error("Cos() failed on MPmap. Produced:", localarr.in, "Expected:", localarr.out)
	}
}

//Test Smap with Cos()
func TestSmapCos(t *testing.T) {
	localarr := genPair(10000, math.Cos)
	Smap(math.Cos, localarr.in, 0, len(localarr.in))
	if localarr.match() {
		return
	} else {
		t.Error("Cos() failed on Smap. Produced:", localarr.in, "Expected:", localarr.out)
	}
}

//Test PPmap with Cos()
func TestPPmapCos(t *testing.T) {
	localarr := genPair(10000, math.Cos)
	PPmap(math.Cos, localarr.in)
	if localarr.match() {
		return
	} else {
		t.Error("Cos() failed on PPmap. Produced:", localarr.in, "Expected:", localarr.out)
	}
}
