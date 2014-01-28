package vec

type empty struct{}
type semaphore chan empty

func NewSem(N int) semaphore {
	s := make(semaphore, N)
	return s
}

//semaphore.Signal()
//sends a signal on the 'semaphore' channel
func (s semaphore) Signal() {
	e := empty{}
	s <- e
}

//semaphore.Wait(n int)
//blocks until 'n' Signal()s are sent
func (s semaphore) Wait(n int) {
	for i := 0; i < n; i++ {
		<-s
	}
}
