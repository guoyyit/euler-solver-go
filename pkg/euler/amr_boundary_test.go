package euler

import (
	"math"
	"testing"
)

func TestAMRReflectingWallConservative(t *testing.T) {
	cfg := Config{
		X0: 0, X1: 1, NBase: 160,
		Boundary: BoundaryReflecting,
		Flux:     FluxHLL, Limiter: LimitMC, CFL: 0.7, Gamma: 1.4,
		AMR: DefaultAMRConfig(),
	}
	cfg.AMR.MaxLevels = 2
	cfg.AMR.RefineThreshold = 0.08
	cfg.AMR.CoarsenThreshold = 0.01
	ic := func(x float64) Primitive {
		q := Primitive{Rho: 1, U: 0, P: 1}
		if x > 0.2 && x < 0.4 {
			q.P = 2.5
		}
		return q
	}
	s, err := NewSolver(cfg, ic)
	if err != nil {
		t.Fatal(err)
	}
	r0, _, e0 := s.TotalConserved()
	if err := s.Advance(1.2); err != nil {
		t.Fatal(err)
	}
	r1, _, e1 := s.TotalConserved()
	// At impermeable walls MASS and ENERGY are conserved (zero mass and zero
	// advective energy flux); momentum is not, because the wall supplies an
	// impulse.
	tol := 1e-9
	if math.Abs(r1-r0) > tol || math.Abs(e1-e0) > 1e-8 {
		t.Fatalf("AMR wall run not mass/energy conservative: dR=%.2e dE=%.2e", r1-r0, e1-e0)
	}
	for _, lv := range s.Mesh().Levels {
		for i := NGhost; i < NGhost+lv.NCells(); i++ {
			if lv.Q[i].Rho <= 0 || lv.Q[i].P <= 0 {
				t.Fatalf("negative state dx=%v: %+v", lv.DX, lv.Q[i])
			}
		}
	}
}
