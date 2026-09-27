package euler

import (
	"math"
	"testing"
)

var sodExact = &RiemannData{
	PL:    Primitive{Rho: 1, U: 0, P: 1},
	PR:    Primitive{Rho: 0.125, U: 0, P: 0.1},
	Gamma: 1.4,
	X0:    0.5,
}

func runSod(t *testing.T, n int) *Simulation {
	t.Helper()
	sim, err := NewSimulation(testConfig(n, BoundaryReflecting))
	if err != nil {
		t.Fatalf("init: %v", err)
	}
	if err := sim.Run(); err != nil {
		t.Fatalf("run: %v", err)
	}
	return sim
}

// L1 errors are integrated over cells, comparing cell averages to the exact
// solution sampled at each cell center (away from discontinuities this is
// effectively the exact cell average).
func sodL1(sim *Simulation) (rho, p float64) {
	g := sim.Cfg.gas()
	var dom float64
	for _, c := range sim.Grid.Cells {
		pr, _ := g.ConservativeToPrimitive(c.U)
		ex := sodExact.ExactState(c.Center(), sim.Time)
		w := c.Width()
		rho += math.Abs(pr.Rho-ex.Rho) * w
		p += math.Abs(pr.P-ex.P) * w
		dom += w
	}
	return rho / dom, p / dom
}

func TestSodMatchesExact(t *testing.T) {
	sim := runSod(t, 512)
	eRho, eP := sodL1(sim)
	if eRho > 0.01 || eP > 0.01 {
		t.Fatalf("L1 error too large: rho=%g p=%g", eRho, eP)
	}
}

func TestSodConvergesUnderRefinement(t *testing.T) {
	var prev float64
	for _, n := range []int{128, 256, 512} {
		sim := runSod(t, n)
		eRho, _ := sodL1(sim)
		if prev > 0 && eRho >= prev {
			t.Fatalf("density L1 did not decrease: %g -> %g", prev, eRho)
		}
		prev = eRho
	}
}

// The minmod limiter is TVD: no cell value may exceed the initial extrema
// (1.0 density high, 0.125 low; etc.). This rules out Gibbs-style ringing.
func TestSodNoSpuriousExtrema(t *testing.T) {
	sim := runSod(t, 512)
	g := sim.Cfg.gas()
	rhoMin, rhoMax := math.Inf(1), math.Inf(-1)
	pMin, pMax := math.Inf(1), math.Inf(-1)
	for _, c := range sim.Grid.Cells {
		pr, err := g.ConservativeToPrimitive(c.U)
		if err != nil {
			t.Fatalf("negative state: %v", err)
		}
		rhoMin = math.Min(rhoMin, pr.Rho)
		rhoMax = math.Max(rhoMax, pr.Rho)
		pMin = math.Min(pMin, pr.P)
		pMax = math.Max(pMax, pr.P)
	}
	if rhoMax > 1+1e-9 || rhoMin < 0.125-1e-9 {
		t.Fatalf("density overshoot [%.6f,%.6f]", rhoMin, rhoMax)
	}
	if pMax > 1+1e-9 || pMin < 0.1-1e-9 {
		t.Fatalf("pressure overshoot [%.6f,%.6f]", pMin, pMax)
	}
}
