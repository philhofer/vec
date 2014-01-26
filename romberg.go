package vec

import "math"

//Trapezoidal rule - 'f(x)' from 'a' to 'b' with N steps
func trap(f Mathop, a float64, b float64, N int) float64 {
	if N <= 0 { return a }
	if N == 1 { return (f(b) - f(a))/2 }
	
	out := 0.0
	h := (b-a)/float64(N)
	x1s := Arange(a, b, N)
	x2s := Arange(a+h, b+h, N)
	PPmap(f, x1s)
	PPmap(f, x2s)
	for i := range x1s {
		out += x1s[i] + x2s[i]
	}
	out -= (h/2)*(f(a) + f(b))
	return out*(h/2.0)
}

/*
Performs Romberg integration on a Mathop

Returns the integral evaluated from 'a' to 'b', and error

If error >> 0, then it's possible the integral does not converge
*/
func Integral(f Mathop, a float64, b float64) (out float64, err float64) {
	err = math.Inf(1)
	out = 0.0
	const K int = 10
	if a == b { return }
	if math.IsNaN(a) || math.IsNaN(b) {
		return
	}

	if math.IsInf(a, 0) || math.IsInf(b, 0) {
		//TODO: transform f, a, b
	}

	var Ip []float64
	for k:=0; k<K; k++{
		Ik := make([]float64, K-k)
		
		if k == 0 {
			for i := range Ik {
				Ik[i] = trap(f, a, b, int(math.Pow(2, float64(i))))
			}
		} else {
			for i := range Ik {
				j := i+1
				m := math.Pow(4.0, float64(k))
				Ik[i] = (m*Ip[j] - Ip[j-1])/(m-1.0)
			}
		}

		if len(Ik)>1 {
			err = math.Abs(Ik[K-k-2] - Ik[K-k-1])
		}
		out = Ik[K-k-1]
		
		if err < 2E-16 {
			return out, err
		}

		Ip = Ik
	}

	return out, err
}
