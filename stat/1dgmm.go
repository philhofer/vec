package stat

import "math"
// 1-dimensional Gaussian Mixture Modelling

var EMAcc float64 = 10E-8

type Gaussian struct {
	Mean float64
	StDev float64
}

func (g Gaussian) PDF(x float64) float64 {
	n := &Normal{mean: g.Mean, sigma: g.StDev}
	return n.PDF(x)
}

func UnivariateGMM(arr []float64, K int) ([]Gaussian, []float64) {
	N := len(arr)
	out := make([]Gaussian, K)
	dataSig := StDev(arr)
	dataMean := Mean(arr)
	estSig := dataSig/float64(K)
	ps := make([]float64, K) //probabilites
	const MAXITER = 100000 //set maximum iterations to prevent hanging

	//set initial ps
	for i:=0; i<K; i++ {
		ps[i] = 1.0/float64(K)
	}

	//set starting points - spaced equally along dataset
	for i:=0; i<K; i++ {
		thisS := float64(i - K/2)
		out[i] = Gaussian{Mean: dataMean+(thisS*estSig), StDev: estSig}
	}
	
	//mixture pdf
	mix := func(i int) (p float64) {
		for k:=0; k<K; k++ {
			p += ps[k]*out[k].PDF(arr[i])
		}
		return
	}

	//'responsibility' func
	res := func(i int, k int) float64 {
		return out[k].PDF(arr[i])*ps[k]/mix(i)
	}

	//iterate
	niter := 0
	nef := make([]float64, K)
	for niter < MAXITER {

		oldmeans := make([]float64, K)
		oldsigms := make([]float64, K)
		for k:=0; k<K; k++ {
			oldmeans[k] = out[k].Mean
			oldsigms[k] = out[k].StDev
		}

		//get (n | k) and new mean estimate
		for k:=0; k<K; k++ {
			nef[k] = 0.0
			uk := 0.0
			for i:=0;i<N;i++{
				rik := res(i, k)
				nef[k] += rik
				uk += rik*arr[i]
			}
			out[k].Mean = uk/nef[k]
			ps[k] = nef[k]/float64(N)
		}

		//get new sigma estimates
		for k:=0; k<K; k++ {
			sk := 0.0
			nef[k] = 0.0
			for i:=0;i<N;i++ {
				rik := res(i, k)
				ndif := arr[i] - out[k].Mean
				sk += rik*ndif*ndif
				nef[k] += rik
			}
			out[k].StDev = math.Sqrt(sk/nef[k])
			ps[k] = nef[k]/float64(N)
		}

		niter++
		//check for convergence
		conv := 0
		for k:=0; k<K; k++ {
			deltaM := math.Abs(oldmeans[k] - out[k].Mean)
			deltaS := math.Abs(oldsigms[k] - out[k].StDev)
			if deltaM > EMAcc || deltaS > EMAcc {
				conv++
			}
		}
		if conv == 0 {
			return out, ps
		}
	}
	return out, ps
}
