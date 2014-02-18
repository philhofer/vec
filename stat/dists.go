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
func UnitNormalDist() UnivariateDistribution {
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
	return 0.5*(1.0 + math.Erf((x - n.mean)/(math.Sqrt2*n.sigma)))
}

//To create a Normal object
func NormalDist(mean float64, sigma float64) UnivariateDistribution {
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
	if math.IsInf(x, 1) {
		return 1.0
        }

	if x <= 0 {return 0}
	return e.lambda*math.Exp(-e.lambda*x)
}

//Exponential CDF
func (e *Exponential) CDF(x float64) float64 {
	if x <= 0 { return 0 }
	if math.IsInf(x, 1) {
		return 1.0
        }
	return 1.0 - math.Exp(-e.lambda*x)
}

//To create an Exponential struct
func ExponentialDist(lambda float64) UnivariateDistribution {
	if lambda <= 0 {
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
	return (1.0/(12.0*n)) - (1.0/(360.0*n*n*n)) + (1.0/(1260.0*n*n*n*n*n)) - (1.0/(1680.0*math.Pow(n, 7)))
}

func d0(x float64) float64 {
	return x*math.Log(x) + 1.0 - x
}

func (p *Poisson) PDF(k int) float64 {
	//Use:
	//http://projects.scipy.org/scipy/raw-attachment/ticket/620/loader2000Fast.pdf

	x := float64(k)
	return (1.0/math.Sqrt(2.0*math.Pi*x))*math.Exp(-delta(x) - p.lambda*d0(x/p.lambda))
}

func (p *Poisson) CDF(x float64) float64 {
	if math.IsInf(x, 1) {
		return 1.0
	} else if math.IsInf(x, -1) {
		return 0.0
	}

	//Use regularized incomplete gamma function from specialfuncs
	k := math.Floor(x)+1.0
	return G(k, p.lambda)
}

func PoissonDist(lambda float64) DiscreetDistribution {
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
	if math.IsInf(x, 1) {
		return 1.0
	} else if math.IsInf(x, -1) {
		return 0.0
	}
	lga, _ := math.Lgamma(g.alpha)
	logp := (g.alpha * math.Log(g.beta)) + ((g.alpha - 1.0) * math.Log(x)) - (x * g.beta) - lga
	return math.Exp(logp)
}

//Gamma Distribution CDF
func (g *Gamma) CDF(x float64) float64 {
	if x <= 0 { return math.NaN() }
	if math.IsInf(x, 1) {
		return 1.0
	} else if math.IsInf(x, -1) {
		return 0.0
	}
	return P(g.alpha, g.beta*x)
}

//Gamma Distribution Constructor - uses 'shape' and 'rate' parameters
func GammaDist(alpha float64, beta float64) UnivariateDistribution {
	if alpha < 0 || beta < 0 {
		return nil
	}
	return &Gamma{alpha, beta}
}

//Chi-squared Distribution Constructor
//(Special case of Gamma Distribution Constructor)
func ChiSquaredDist(v int) UnivariateDistribution {
	if v <= 0 {
		return nil
	}
	alpha := float64(v)/2.0
	beta := 0.5
	return &Gamma{alpha, beta}
}

//Student-T distribution
type StudentT struct {
	nu float64
}

func (s *StudentT) PDF(x float64) float64 {
	if math.IsNaN(x) { return math.NaN() }
	if math.IsInf(x, 1) {
		return 1.0
	} else if math.IsInf(x, -1) {
		return 0.0
	}

	return 1.0/(math.Sqrt(s.nu) * Beta(0.5, s.nu/2.0)) * math.Pow(1.0 + (x*x)/s.nu, 2.0)
}

func (s *StudentT) CDF(t float64) float64 {
	if math.IsNaN(t) { return math.NaN() }
	if math.IsInf(t, 1) {
		return 1.0
	} else if math.IsInf(t, -1) {
		return 0.0
	}

	if t < 0 {
		return 1.0 - s.CDF(-t)
	}

	x := s.nu/(t*t + s.nu)
	return 1.0 - 0.5*IncBeta(x, s.nu/2.0, 0.5)
}

func StudentTDist(nu float64) UnivariateDistribution {
	if nu < 0 {
		return nil
	} else if math.IsInf(nu, 1) {
		return &UnitNormal{}
	}
	return &StudentT{nu}
}

//CDF of KS-distribution
func PKS(z float64) float64 {
	if math.IsNaN(z) {
		return math.NaN()
	} else if math.IsInf(z, 1) {
		return 1.0
	} else if math.IsInf(z, -1) {
		return 0.0
	}

	out := math.Sqrt(2.0*math.Pi)/z
	sum := 0.0
	for j:=1; j<15; j++ {
		j := float64(j)
		sum += math.Exp(-1.0*(2*j-1.0)*(2*j-1.0)*math.Pi*math.Pi / (8.0*z*z))
	}
	return out*sum
}

//1-CDF of KS-distribution 
func QKS(z float64) float64 {
	return 1.0-PKS(z)
}
