package euler

import (
	"math"
	"testing"
)

func TestDoubleRarefactionPositive(t *testing.T) {
	p := GetProblem(ProblemDoubleRarefaction)
	cfg := Config{
		X0: 0, X1: 1, NBase: 256,
		Boundary: BoundaryReflecting,
		Flux:     FluxRusanov, Limiter: LimitMinmod, CFL: 0.8, Gamma: 1.4,
	}
	s, err := NewSolver(cfg, p.Init)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Advance(0.15); err != nil {
		t.Fatal(err)
	}
	lv := s.Mesh().Levels[0]
	minP, minR := math.Inf(1), math.Inf(1)
	for i := NGhost; i < NGhost+lv.NCells(); i++ {
		if lv.Q[i].P < minP {
			minP = lv.Q[i].P
		}
		if lv.Q[i].Rho < minR {
			minR = lv.Q[i].Rho
		}
	}
	if minR <= 0 || minP <= 0 {
		t.Fatalf("non-positive state: minRho=%v minP=%v", minR, minP)
	}
	// Total energy conserved even through the near-vacuum region.
	r0, _, e0 := totalsAt(t, cfg, p.Init)
	r1, _, e1 := s.TotalConserved()
	tol := 1e-10 * math.Max(1, math.Abs(e0))
	if math.Abs(r1-r0) > tol || math.Abs(e1-e0) > tol {
		t.Fatalf("conservation lost: r %v->%v e %v->%v", r0, r1, e0, e1)
	}
}

func totalsAt(t *testing.T, cfg Config, ic func(x float64) Primitive) (r, m, e float64) {
	t.Helper()
	s, err := NewSolver(cfg, ic)
	if err != nil {
		t.Fatal(err)
	}
	return s.TotalConserved()
}

// A wave must not gain or lose energy bouncing off reflecting walls over many
// wave-transit times.
func TestReflectingLongRunStable(t *testing.T) {
	cfg := Config{
		X0: 0, X1: 1, NBase: 128,
		Boundary: BoundaryReflecting,
		Flux:     FluxHLL, Limiter: LimitMC, CFL: 0.7, Gamma: 1.4,
	}
	// A compact high-pressure blob in the middle.
	ic := func(x float64) Primitive {
		q := Primitive{Rho: 1, U: 0, P: 1}
		if math.Abs(x-0.5) < 0.08 {
			q.P = 2
		}
		return q
	}
	s, err := NewSolver(cfg, ic)
	if err != nil {
		t.Fatal(err)
	}
	r0, m0, e0 := s.TotalConserved()
	if err := s.Advance(2.0); err != nil {
		t.Fatal(err)
	}
	r1, m1, e1 := s.TotalConserved()
	tol := 1e-9
	check := func(name string, a, b float64) {
		if math.Abs(a-b) > tol {
			t.Fatalf("%s drift %.3e", name, math.Abs(a-b))
		}
	}
	check("mass", r0, r1)
	check("momentum", m0, m1)
	check("energy", e0, e1)
	lv := s.Mesh().Levels[0]
	for i := NGhost; i < NGhost+lv.NCells(); i++ {
		if lv.Q[i].Rho <= 0 || lv.Q[i].P <= 0 {
			t.Fatalf("negative after long run: %+v", lv.Q[i])
		}
	}
}
