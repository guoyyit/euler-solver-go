package euler

import (
	"math"
	"testing"
)

func sodCfg(n int) Config {
	return Config{
		X0: 0, X1: 1, NBase: n,
		Boundary: BoundaryPeriodic,
		Flux:     FluxHLL,
		Limiter:  LimitMinmod,
		CFL:      0.8,
		Gamma:    1.4,
	}
}

func TestUniformSodConservationAndPositivity(t *testing.T) {
	p := GetProblem(ProblemSod)
	s, err := NewSolver(sodCfg(200), p.Init)
	if err != nil {
		t.Fatal(err)
	}
	r0, m0, e0 := s.TotalConserved()
	if err := s.Advance(0.2); err != nil {
		t.Fatal(err)
	}
	r1, m1, e1 := s.TotalConserved()
	tol := 1e-12 * math.Max(1, math.Max(math.Abs(e0), math.Abs(r0)))
	if math.Abs(r1-r0) > tol {
		t.Fatalf("mass drift: %.3e > %.3e", math.Abs(r1-r0), tol)
	}
	if math.Abs(m1-m0) > tol {
		t.Fatalf("momentum drift: %.3e", math.Abs(m1-m0))
	}
	if math.Abs(e1-e0) > tol {
		t.Fatalf("energy drift: %.3e", math.Abs(e1-e0))
	}
	// No negative states; densities/pressures bounded by the two initial
	// states (up to a small tolerance for MUSCL roundoff).
	lv := s.Mesh().Levels[0]
	for i := NGhost; i < NGhost+lv.NCells(); i++ {
		q := lv.Q[i]
		if q.Rho <= 0 || q.P <= 0 {
			t.Fatalf("negative state at %d: %+v", i, q)
		}
		if q.Rho > 1.0+1e-10 || q.Rho < 0.125-1e-10 {
			t.Fatalf("density out of range: %v", q.Rho)
		}
	}
}

// shockIndex returns the packed-cell position of the right-going shock. Across
// a right-going shock the density jumps UP, whereas across the contact it jumps
// down, so locating the strongest positive density rise in the right half
// selects the shock rather than the contact.
func shockIndex(lv *Level, mesh *Mesh) int {
	best, bestAt := 0.0, NGhost
	start := NGhost + lv.NCells()/2
	for i := start; i < NGhost+lv.NCells(); i++ {
		if d := lv.Q[i].Rho - lv.Q[i-1].Rho; d > best {
			best, bestAt = d, i
		}
	}
	return bestAt
}

func TestUniformSodShockConvergence(t *testing.T) {
	p := GetProblem(ProblemSod)
	// High-resolution reference sampled on its own grid; compare coarser
	// solutions against it by piecewise-average L1 density error.
	solve := func(n int) []float64 {
		s, err := NewSolver(sodCfg(n), p.Init)
		if err != nil {
			t.Fatal(err)
		}
		if err := s.Advance(0.2); err != nil {
			t.Fatal(err)
		}
		lv := s.Mesh().Levels[0]
		out := make([]float64, n)
		for i := 0; i < n; i++ {
			out[i] = lv.Q[NGhost+i].Rho
		}
		return out
	}
	ref := solve(800)
	l1 := func(n int) float64 {
		q := solve(n)
		factor := len(ref) / n // integer for 100/200 vs 800
		var sum float64
		for i := 0; i < n; i++ {
			// average reference over the coarse cell
			var ravg float64
			for k := 0; k < factor; k++ {
				ravg += ref[i*factor+k]
			}
			ravg /= float64(factor)
			sum += math.Abs(q[i]-ravg) / float64(n)
		}
		return sum
	}
	e100 := l1(100)
	e200 := l1(200)
	if !(e200 < e100) {
		t.Fatalf("density L1 error not decreasing: e100=%v e200=%v", e100, e200)
	}
	// A shock-capturing TVD scheme converges at roughly first order in L1;
	// require a clear reduction (rate > ~0.5).
	if rate := math.Log(e100/e200) / math.Log(2); rate < 0.5 {
		t.Fatalf("convergence rate too low: %.3f", rate)
	}
}

// Running the same configuration twice must give bit-identical totals.
func TestDeterminism(t *testing.T) {
	p := GetProblem(ProblemLax)
	run := func() (r, m, e float64) {
		s, err := NewSolver(sodCfg(160), p.Init)
		if err != nil {
			t.Fatal(err)
		}
		if err := s.Advance(0.16); err != nil {
			t.Fatal(err)
		}
		return s.TotalConserved()
	}
	r1, m1, e1 := run()
	r2, m2, e2 := run()
	if r1 != r2 || m1 != m2 || e1 != e2 {
		t.Fatalf("runs differ: (%.17g,%.17g,%.17g) vs (%.17g,%.17g,%.17g)",
			r1, m1, e1, r2, m2, e2)
	}
}
