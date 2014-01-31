package stat

import (
	"math"
)
/* Distributions (to be) defined in this source:

  -- Normal
  -- Unit Normal
  -- Exponential
  -- Poisson
  -- Binomial
  -- Negative Binomial
  -- Cauchy
  -- Geometric
  -- Laplace
  -- Chi-Squared
  -- Beta
  -- Student-T

*/


// INTERFACE FOR UNIVARIATE DISTRIBUTIONS
type UnivariateDistribution interface {
	PDF(x float64) float64
	CDF(x float64) float64
}

type DiscreetDistribution interface {
	PDF(k int) float64
	CDF(x float64) float64
}

// Unit Normal distribution type
type UnitNormal struct {

}

// Unit Normal PDF
func (u *UnitNormal) PDF(x float64) float64 {
	return (1.0/math.Sqrt(2.0*math.Pi))*(math.Exp(-0.5*x*x))
}

//Unit Normal CDF
func (u *UnitNormal) CDF(x float64) float64 {
	return 0.5*(1.0 + math.Erf(x/(math.Sqrt2)))
}

//To create a UnitNormal struct
func UnitNormalDist() *UnitNormal {
	return &UnitNormal{}
}

//Definition of Normal Distribution struct
type Normal struct {
	mean float64
	sigma float64
}

//Normal Distribution PDF
func (n *Normal) PDF(x float64) float64 {
	return (1.0/(n.sigma*math.Sqrt(2.0*math.Pi)))*math.Exp(-(x-n.mean)*(x-n.mean)/(2*n.sigma*n.sigma))
}

//Normal Distribution CDF
func (n *Normal) CDF(x float64) float64 {
	return 0.5*(1.0 + math.Erf((x - n.mean)/(math.Sqrt(2*n.sigma*n.sigma))))
}

//To create a Normal struct
func NormalDist(mean float64, sigma float64) *Normal {
	if sigma <= 0 {
		return nil
	}
	return &Normal{mean, sigma}
}

//Definition of Exponential struct
type Exponential struct {
	lambda float64
}

//Exponential PDF
func (e *Exponential) PDF(x float64) float64 {
	if x < 0 {
		return 0
	}
	return e.lambda*math.Exp(-e.lambda*x)
}

//Exponential CDF
func (e *Exponential) CDF(x float64) float64 {
	if x < 0 {
		return 0
	}
	return 1.0 - math.Exp(-e.lambda*x)
}

//To create an Exponential struct
func ExponentialDist(lambda float64) *Exponential {
	if lambda >= 0 {
		return nil
	}
	return &Exponential{lambda}
}


//Poisson type struct
type Poisson struct {
	lambda float64
}

//Stirling -  De Moivre Series - accurate to O(n^-7)
func delta(n float64) float64 {
	return (1.0/(12.0*n)) - (1.0/(360.0*n*n*n)) + (1.0/(1260.0*n*n*n*n*n))
}

func D0(x float64) float64 {
	return x*math.Log(x) + 1.0 - x
}

func (p *Poisson) PDF(k int) float64 {
	//Use:
	//http://projects.scipy.org/scipy/raw-attachment/ticket/620/loader2000Fast.pdf
	
	x := float64(k)
	return (1.0/math.Sqrt(2.0*math.Pi*x))*math.Exp(-delta(x) - p.lambda*D0(x/p.lambda))
}

func (p *Poisson) CDF(x float64) float64 {
	//Use normalized incomplete gamma function from specialfuncs
	return G(math.Floor(x), p.lambda)
}

func PoissonDist(lambda float64) *Poisson {
	if lambda<=0 {
		return nil
	} 
	return &Poisson{lambda}
}

type Gamma struct {
	alpha float64
	beta float64
}

//Gamma Distribution PDF
func (g *Gamma) PDF(x float64) float64 {
	if x <= 0 { return math.NaN() }
	logp := g.alpha*math.Log(g.beta) + (g.alpha - 1.0)*math.Log(x) - x*g.beta - math.Lgamma(g.alpha)
	return math.Exp(logp)
}

//Gamma Distribution CDF
func (g *Gamma) CDF(x float64) float64 {
	if x <= 0 { return math.NaN() }
	return P(g.alpha, g.beta*x)
}

//Gamma Distribution Constructor
func GammaDist(alpha float64, beta float64) *Gamma {
	if alpha < 0 || beta < 0 {
		return nil
	}
	return &Gamma{alpha, beta}
}

//Chi-squared Distribution Constructor
//(Special case of Gamma Distribution Constructor)
func ChiSquaredDist(v int) {
	if v <= 0 {
		return nil
	}
	alpha := float64(v)/2.0
	beta := 0.5
	return &Gamma{alpha, beta}
}
