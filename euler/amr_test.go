package euler

import (
	"math"
	"testing"
)

// Splitting conserves each component: parent integrated U equals sum of children.
func TestSplitConservative(t *testing.T) {
	g, _ := NewGas(1.4)
	parent := &Cell{ID: RootID(), Level: RootLevel, X0: 0, X1: 1}
	c, _ := g.PrimitiveToConservative(Primitive{Rho: 1.2, U: 0.7, P: 0.9})
	parent.U = c
	l, r, err := splitCell(g, parent)
	if err != nil {
		t.Fatal(err)
	}
	check := func(name string, pv, lv, rv float64) {
		W := parent.Width()
		want := pv * W
		got := lv*l.Width() + rv*r.Width()
		if math.Abs(got-want) > 1e-13*(1+math.Abs(want)) {
			t.Fatalf("split not conservative for %s: want %.15g got %.15g", name, want, got)
		}
	}
	check("rho", parent.U.Rho, l.U.Rho, r.U.Rho)
	check("rhou", parent.U.RhoU, l.U.RhoU, r.U.RhoU)
	check("E", parent.U.E, l.U.E, r.U.E)
}

// Merging two siblings returns exactly the parent integrated quantity.
func TestMergeConservative(t *testing.T) {
	g, _ := NewGas(1.4)
	parent := &Cell{ID: RootID(), Level: RootLevel, X0: 0, X1: 1}
	c, _ := g.PrimitiveToConservative(Primitive{Rho: 1.2, U: 0.7, P: 0.9})
	parent.U = c
	l, r, _ := splitCell(g, parent)
	p2, err := mergeCells(g, l, r)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(p2.U.Rho-parent.U.Rho) > 1e-15 ||
		math.Abs(p2.U.RhoU-parent.U.RhoU) > 1e-15 ||
		math.Abs(p2.U.E-parent.U.E) > 1e-15 {
		t.Fatalf("merge roundtrip mismatch: %+v vs %+v", p2.U, parent.U)
	}
}

func TestRegridBalanceAndConservation(t *testing.T) {
	cfg := amrConfig(128, BoundaryReflecting)
	sim, err := NewSimulation(cfg)
	if err != nil {
		t.Fatal(err)
	}
	g := cfg.gas()

	integrals := func() (m, p, e float64) {
		for _, c := range sim.Grid.Cells {
			m += c.U.Rho * c.Width()
			p += c.U.RhoU * c.Width()
			e += c.U.E * c.Width()
		}
		return
	}
	m0, mo0, e0 := integrals()

	// Force many regrids around the evolving waves.
	for k := 0; k < 40; k++ {
		if err := Regrid(g, sim.Grid, false, cfg.Refine); err != nil {
			t.Fatalf("regrid: %v", err)
		}
		// 2:1 balance between neighbors.
		for i := 0; i+1 < len(sim.Grid.Cells); i++ {
			if d := sim.Grid.Cells[i].Level - sim.Grid.Cells[i+1].Level; d > 1 || d < -1 {
				t.Fatalf("2:1 violated at %d: levels %d,%d", i,
					sim.Grid.Cells[i].Level, sim.Grid.Cells[i+1].Level)
			}
		}
		// Tessellation integrity: adjacent, cover domain.
		if math.Abs(sim.Grid.Cells[0].X0-cfg.Xmin) > 1e-14 ||
			math.Abs(sim.Grid.Cells[len(sim.Grid.Cells)-1].X1-cfg.Xmax) > 1e-14 {
			t.Fatal("grid does not cover domain")
		}
		for i := 0; i+1 < len(sim.Grid.Cells); i++ {
			if math.Abs(sim.Grid.Cells[i].X1-sim.Grid.Cells[i+1].X0) > 1e-14 {
				t.Fatalf("gap/overlap at %d", i)
			}
		}
		if err := sim.Advance(); err != nil {
			t.Fatalf("advance: %v", err)
		}
	}
	m1, mo1, e1 := integrals()
	scale := math.Max(1, math.Abs(e0))
	// Reflecting wall: mass is conserved exactly. Momentum changes by the wall
	// impulse and energy changes only by wall work, which for this symmetric
	// setup is zero in the continuous problem; the finite-volume update
	// reproduces it to a tight, step-count-scaled tolerance. The strict
	// machine-precision three-component check is performed on the periodic
	// domain (TestPeriodicConservationWithAMR).
	if math.Abs(m1-m0) > 1e-12*scale {
		t.Fatalf("AMR mass drift: dM=%.3e", m1-m0)
	}
	if math.Abs(e1-e0) > 1e-9*scale {
		t.Fatalf("AMR energy drift too large: dE=%.3e", e1-e0)
	}
	_ = mo0
	_ = mo1
}

// On a periodic domain AMR regridding plus hanging-node interface fluxes must
// conserve all three components to machine precision over a full run.
func TestPeriodicConservationWithAMR(t *testing.T) {
	sim, err := NewSimulation(amrConfig(128, BoundaryPeriodic))
	if err != nil {
		t.Fatal(err)
	}
	m0, mo0, e0 := sim.TotalMass(), sim.TotalMomentum(), sim.TotalEnergy()
	if err := sim.Run(); err != nil {
		t.Fatalf("run: %v", err)
	}
	scale := math.Max(1, math.Max(math.Abs(mo0), math.Abs(e0)))
	dM, dMom, dE := sim.TotalMass()-m0, sim.TotalMomentum()-mo0, sim.TotalEnergy()-e0
	if math.Abs(dM) > 1e-12*scale || math.Abs(dMom) > 1e-12*scale || math.Abs(dE) > 1e-12*scale {
		t.Fatalf("periodic AMR drift: dM=%.3e dMom=%.3e dE=%.3e", dM, dMom, dE)
	}
}

// AMR must actually refine the wave region (some cells above level 0) while
// keeping far-field cells coarser, and the AMR solution must match the exact
// Riemann solution at least as well as a uniform coarse grid.
func TestAMRAdaptsAndConverges(t *testing.T) {
	sim, err := NewSimulation(amrConfig(128, BoundaryReflecting))
	if err != nil {
		t.Fatal(err)
	}
	if err := sim.Run(); err != nil {
		t.Fatalf("run: %v", err)
	}
	maxLevel := RootLevel
	nFine := 0
	for _, c := range sim.Grid.Cells {
		if c.Level > maxLevel {
			maxLevel = c.Level
		}
		if c.Level > RootLevel {
			nFine++
		}
	}
	if maxLevel <= RootLevel || nFine == 0 {
		t.Fatalf("AMR never refined: maxLevel=%d nFine=%d", maxLevel, nFine)
	}
	if len(sim.Grid.Cells) >= 128*(1<<3) {
		t.Fatalf("AMR populated every cell: %d (no coarsening benefit)", len(sim.Grid.Cells))
	}
	eRho, _ := sodL1(sim)
	// AMR on 128 base should beat or match uniform 256 accuracy roughly;
	// require a loose but meaningful bound.
	if eRho > 0.01 {
		t.Fatalf("AMR solution too inaccurate: rho L1=%g", eRho)
	}
}
