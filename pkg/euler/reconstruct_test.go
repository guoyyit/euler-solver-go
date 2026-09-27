package euler

import "testing"

func TestReconstructLinearExact(t *testing.T) {
	eos := NewEOS(1.4)
	n := 8
	w := make([]Primitive, n)
	for i := range w {
		w[i] = Primitive{Rho: float64(i), U: 1, P: 1}
	}
	lo, hi := 2, n-2
	r := &Reconstructor{EOS: eos, Limiter: LimitNone}
	qL, qR := r.Reconstruct(w, lo, hi)
	for j := lo - 1; j < hi; j++ {
		if d := abs(qL[j].Rho - (float64(j) + 0.5)); d > 1e-12 {
			t.Fatalf("qL[%d].rho=%v want %v", j, qL[j].Rho, float64(j)+0.5)
		}
		if d := abs(qR[j].Rho - (float64(j) + 0.5)); d > 1e-12 {
			t.Fatalf("qR[%d].rho=%v want %v", j, qR[j].Rho, float64(j)+0.5)
		}
	}
}

func TestReconstructMinmodMonotone(t *testing.T) {
	eos := NewEOS(1.4)
	w := []Primitive{
		{Rho: 1, P: 1}, {Rho: 1, P: 1}, {Rho: 1, P: 1},
		{Rho: 2, P: 2}, {Rho: 2, P: 2}, {Rho: 2, P: 2},
	}
	r := &Reconstructor{EOS: eos, Limiter: LimitMinmod}
	qL, qR := r.Reconstruct(w, 2, 4)
	if qL[2].Rho != 1 || qR[2].Rho != 2 {
		t.Fatalf("minmod overshoot at jump: qL=%v qR=%v", qL[2], qR[2])
	}
}
