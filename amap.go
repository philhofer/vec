package vec

/* Massively Parallel Function-Slice Mapping
Should be used for mapping 'expensive' functions
onto short(ish) arrays in-place. Creates as many goroutines
as there are floats in the array, so this bogs down for 
long arrays.
*/ 
func MPmap(fm Mathop, arr []float64) {
	sem := NewSem(len(arr))
	for i := range arr {
		go func (i int){
			arr[i] = fm(arr[i])
			sem.Signal()
		}(i)
	}
	sem.Wait(len(arr))
}
