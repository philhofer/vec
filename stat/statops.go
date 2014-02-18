package stat

import (
	"math"
	"sort"
)

//Welch's t-test
//Unequal sample size, unequal variance, 2-tailed
func WelchTtest(A []float64, B []float64) (p float64, t float64) {
	//test for proper initial conditions
	if len(A) <= 1 || len(B) <= 1 {
		return 1.0, 0.0
	}

	xbarA := Mean(A)
	xbarB := Mean(B)
	sA := Variance(A)/float64(len(A))
	sB := Variance(B)/float64(len(B))
	s := math.Sqrt(sA + sB)
	t = math.Abs(xbarA - xbarB)/s
	df := math.Pow(sA + sB, 2)/((sA*sA)/float64(len(A)-1) + (sB*sB)/float64(len(B)-1))
	Dist := StudentTDist(df)
	if Dist == nil {
		panic("Problem in Welch's t-test...")
	}
	p = (1.0 - Dist.CDF(t))*2.0
	return
}

func Ttest(mean float64, A []float64) (p float64, t float64) {
	t = math.Abs(Mean(A) - mean)/(StDev(A)/math.Sqrt(float64(len(A))))
	df := len(A) - 1
	Dist := StudentTDist(float64(df))
	if Dist == nil {
		panic("Problem in Ttest...")
	}
	p = (1.0 - Dist.CDF(t))*2.0
	return
}

//2-sample Kolmogorov-Smirnov test
func KSTest2(A []float64, B []float64) (p float64, D float64) {
	//test for proper initial conditions
	if len(A) <= 1 || len(B) <= 1 {
		return 1.0, 0.0
	}

	Ne := (len(A)*len(B))/(len(A) + len(B))
	EA := ECDF(A)
	EB := ECDF(B)
	pts := make([]float64, len(A)+len(B))
	pts = append(pts, A...)
	pts = append(pts, B...)
	sort.Float64s(pts)
	D = 0.0
	for _, x := range pts {
		thisd := math.Abs(EA(x)-EB(x))
		if thisd > D { D = thisd }
	}
	p = QKS(D*(math.Sqrt(float64(Ne))+0.12+0.11/math.Sqrt(float64(Ne))))
	return
}
