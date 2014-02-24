package stat

import (
	"math"
	"sort"
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

func IncBeta(z float64, a float64, b float64) float64 {
	if z < 0 {
		return 0.0
	}
	r := func(k int) float64 {
		if k % 2 == 0 {
			k = k/2
			kf := float64(k)
			return (kf*(b-kf)*z)/((a + 2*kf - 1.0)*(a + 2.0*kf))
		} else {
			k = k/2
			kf := float64(k)
			return -1.0*((a+kf)*(a+b+kf)*z)/((a+2.0*kf)*(a+2.0*kf+1.0))
		}
	}
	mult := r(20)
	for k:=19; k>0; k-- {
		mult = 1.0 + r(k)/(mult)

	}
	return math.Pow(z, a)*math.Pow(1.0-z, b)/(a*Beta(a, b)*mult)

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

func Mean(arr []float64) float64 {
	sum := 0.0
	if len(arr) == 1 {
		return arr[0]
	}
	for _, x := range arr {
		sum += x
	} 
	return sum/float64(len(arr))
}

func Variance(arr []float64) float64 {
	xbar := Mean(arr)
	if len(arr) <= 1 {
		return 0
	}
	s := 0.0
	for _, x := range arr {
		xsqr := x - xbar
		s += (xsqr*xsqr)
	}
	return s/float64(len(arr)-1)
}

func StDev(arr []float64) float64 {
	return math.Sqrt(Variance(arr))
}

func Covariance(xs []float64, ys []float64) float64 {
	if len(xs) != len(ys) {
		return 0.0
	} 
	N := len(xs)
	difs := make([]float64, N)
	xbar := Mean(xs)
	ybar := Mean(ys)
	for i, x := range xs {
		difs[i] = (x-xbar)*(ys[i]-ybar)
	}
	return Mean(difs)
}

func Skewness(arr []float64) float64 {
	dnm := math.Pow(Variance(arr), 1.5)
	sum := 0.0
	xbar := Mean(arr)
	for _, x := range arr {
		xdif := x - xbar
		sum += (xdif*xdif*xdif)
	}
	return (sum/float64(len(arr)))/dnm
}

func Kurtosis(arr []float64) float64 {
	mu := Mean(arr)
	sigsq := Variance(arr)
	sum := 0.0
	for _, x := range arr {
		xdif := x - mu
		sum += (xdif*xdif*xdif*xdif)
	}
	return (sum/float64(len(arr)))/(sigsq*sigsq)
}

type Condition func(float64) bool

func Which(arr []float64, cond Condition) []float64 {
	out := []float64{}
	for _, x := range arr {
		if cond(x) {
			out = append(out, x)
		}
	}
	return out
}

func Median(arr []float64) float64 {
	n := len(arr)
	var thisarr []float64
	if !sort.Float64sAreSorted(arr) {
		thisarr = make([]float64, n)
		copy(arr, thisarr)
		sort.Float64s(thisarr)
	} else {
		thisarr = arr
	}
	if n%2 == 0 {
		j := n/2
		return (thisarr[j]+thisarr[j+1])/2
	} else {
		return thisarr[(n+1)/2]
	}
}

func ECDF(arr []float64) (func(float64) float64) {
	lessthan := func(n float64) Condition {
		return func(x float64) bool {
			if x < n {
				return true
			} else {
				return false 
			}
		}
	}

	return func(x float64) float64 {
		n := len(Which(arr, lessthan(x)))
		nt := len(arr)
		return float64(n)/float64(nt)
	}
}
