package stat

import (
	"math"
)

/*Beta function
Special cases:
z<0: NaN
w<0: NaN
*/
func Beta(z float64, w float64) float64 {
	if z<0 || w<0 {
		return math.NaN()
	}
	a, _ := math.Lgamma(z)
	b, _ := math.Lgamma(w)
	c, _ := math.Lgamma(z+w)
	return math.Exp(a+b-c)
}


/*Normalized Gamma Function 
(or Complementary Incomplete Gamma Function)

Equal to Gamma(a, x)/Gamma(a)
Evaluated by Legendre's continued fraction

G E [0, 1)
*/
func G(a float64, x float64) float64 {
	if x==0 {
		if a==0 {
			return math.Inf(1)
		} else if a > 0 {
			return 1.0
		} else if a < 0 {
			return 1.0/math.Abs(a)
		}
	}

	//Evaluate Legendre's continued fraction
	//using Lentz's algorithm
	//Shift from Thomson and Barnett
	//Continued fraction:
	//http://functions.wolfram.com/GammaBetaErf/GammaRegularized/10/0003/
	
	b0 := x+1.0-a
	C := 1.0/(10E-30)
	D := 1.0/b0
	if b0 == 0 { D = 10E30 }
	f := D
	//numerator
	an := func(n int) float64 {
		if n==0 { return 1.0 }
		return -1.0*float64(n)*(float64(n)-a)
	}
	//denominator
	bn := func(n int) float64 {
		return x+(float64(2*n+1))-a
	}
	//Lentz's algorithm until machine precision or 1,000 iterations
	for j:=1; j<1000; j++ {
		D = bn(j) + an(j)*D
		if math.Abs(D)<10E-20 {D = 10E-30}
		C = bn(j) + an(j)/C
		if math.Abs(C)<10E-20 {C = 10E-30}
		D = 1.0/D
		del := D*C
		f *= del
		if math.Abs(del-1.0)<10E-15 {
			break
		}
	} 
	lnGa, _ := math.Lgamma(a)
	return f*math.Exp(-x+a*math.Log(x)-lnGa)
}

/*Incomplete Gamma Function P(a, x)
Complement of G(a, x)

P(a, x) E [1, 0)
*/
func P(a float64, x float64) float64 {
	if x < 0.0 {
		return math.NaN()
	} else if a <= 0.0 {
		return math.NaN()
	} else if x == 0.0 {
		return 0.0 
	}
	return 1.0 - G(a, x)
}

//Upper Incomplete Gamma Function
func IncGamma(a float64, x float64) float64 {
	return math.Gamma(a)*G(a, x)
}

//Lower Incomplete Gamma Function
func CompGamma(a float64, x float64) float64 {
	return math.Pow(x, -a)*(P(a, x))
}
