
package stat

import (
	"github.com/philhofer/vec"
	"math/rand"
)

// pseudo-random number generation by CDF

type Generator struct {
	interp *vec.CubicSplineInterpolation
	rgen *rand.Rand
}

func (g *Generator) Next() float64 {
	return g.interp.F(g.rgen.Float64())
}

func (g *Generator) Gen(N int) []float64 {
	if N <= 0 {
		return []float64{}
	}
	out := make([]float64, N)
	for i:=0;i<N;i++ {
		p := g.rgen.Float64()
		out[i] = g.interp.F(p)
	}
	return out
}

func (g *Generator) SetSeed(seed int64) {
	src := rand.NewSource(seed)
	g.rgen = rand.New(src)
	return
}

func NewGenerator(CDF vec.Mathop, mean float64, std float64) *Generator {
	pts := vec.Arange(mean-(10*std), mean+(10*std), 10000)
	cdf := vec.Arange(mean-(10*std), mean+(10*std), 10000)
	vec.PPmap(CDF, cdf)
	bvd := vec.MakeBiVariateData(cdf, pts)
	spl := vec.CubicSpline(bvd)
	src := rand.NewSource(0)
	rnd := rand.New(src)
	return &Generator{interp: spl, rgen: rnd}
}
