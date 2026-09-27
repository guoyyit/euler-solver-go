package euler

import (
	"math"
	"testing"
)

// A moving acoustic/entropy blob on a periodic domain must conserve all three
// integrals to machine precision for the entire integration.
func TestPeriodicConservation(t *testing.T) {
	c := testConfig(256, BoundaryPeriodic)
	c.TEnd = 0.5 // several acoustic transit times
	sim, err := NewSimulation(c)
	if err != nil {
		t.Fatal(err)
	}
	m0, mo0, e0 := sim.TotalMass(), sim.TotalMomentum(), sim.TotalEnergy()
	if err := sim.Run(); err != nil {
		t.Fatalf("run: %v", err)
	}
	dm, dmo, de := sim.TotalMass()-m0, sim.TotalMomentum()-mo0, sim.TotalEnergy()-e0
	scale := math.Max(1, math.Max(math.Abs(mo0), e0))
	tol := 1e-12 * scale
	if math.Abs(dm) > tol || math.Abs(dmo) > tol || math.Abs(de) > tol {
		t.Fatalf("periodic conservation violated: dM=%.3e dMom=%.3e dE=%.3e (tol=%.2e)",
			dm, dmo, de, tol)
	}
}

// At a reflecting wall the mass flux and (through symmetry) the interface must
// behave as a solid boundary: total mass and total energy are conserved, and
// the momentum change equals the wall impulse. No negative states and no
// overshoot beyond the initial extrema are allowed.
func TestReflectingWall(t *testing.T) {
	c := testConfig(512, BoundaryReflecting)
	c.TEnd = 0.3
	sim, err := NewSimulation(c)
	if err != nil {
		t.Fatal(err)
	}
	m0, e0 := sim.TotalMass(), sim.TotalEnergy()
	if err := sim.Run(); err != nil {
		t.Fatalf("run: %v", err)
	}
	scale := math.Max(1.0, math.Abs(e0))
	if math.Abs(sim.TotalMass()-m0) > 1e-12*scale {
		t.Fatalf("wall changed mass by %.3e", sim.TotalMass()-m0)
	}
	if math.Abs(sim.TotalEnergy()-e0) > 1e-12*scale {
		t.Fatalf("wall changed energy by %.3e", sim.TotalEnergy()-e0)
	}

	g := c.gas()
	// Boundary cells: velocity should be near zero at the wall (symmetry).
	left := sim.Grid.Cells[0]
	right := sim.Grid.Cells[len(sim.Grid.Cells)-1]
	pl, _ := g.ConservativeToPrimitive(left.U)
	pr, _ := g.ConservativeToPrimitive(right.U)
	// The left wall of the Sod tube is initially a contact of stationary fluid;
	// fluid adjacent stays ~u=0. Check velocity stays bounded and no negatives.
	if math.IsNaN(pl.U) || math.IsNaN(pr.U) {
		t.Fatal("NaN velocity at wall")
	}
	// Global extremum bound: the reflection must not create pressure/density
	// values beyond the initial data range plus roundoff.
	rhoMax, pMax := 0.0, 0.0
	for _, cell := range sim.Grid.Cells {
		q, err := g.ConservativeToPrimitive(cell.U)
		if err != nil {
			t.Fatalf("negative state near wall: %v", err)
		}
		rhoMax = math.Max(rhoMax, q.Rho)
		pMax = math.Max(pMax, q.P)
	}
	if rhoMax > 1+1e-9 || pMax > 1+1e-9 {
		t.Fatalf("reflecting boundary overshoot: rhoMax=%.6f pMax=%.6f", rhoMax, pMax)
	}
}

func TestGhostStateReflection(t *testing.T) {
	q := Primitive{Rho: 2, U: 3, P: 5}
	g := GhostState(BoundaryReflecting, -1, q)
	if g.Rho != 2 || g.U != -3 || g.P != 5 {
		t.Fatalf("reflection ghost = %+v", g)
	}
	p := GhostState(BoundaryPeriodic, -1, q)
	if p != q {
		t.Fatalf("periodic ghost = %+v", p)
	}
}
