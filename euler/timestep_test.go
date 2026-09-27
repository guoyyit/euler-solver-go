package euler

import "testing"

// Every accepted step must satisfy the CFL condition on every cell:
//
//	dt * (|u|+c) / dx <= CFL
//
// at the beginning of the step (the estimate that controls the update).
func TestCFLNeverViolated(t *testing.T) {
	cfg := amrConfig(256, BoundaryReflecting)
	cfg.CFL = 0.8
	sim, err := NewSimulation(cfg)
	if err != nil {
		t.Fatal(err)
	}
	g := cfg.gas()
	check := func(dt float64) {
		for _, c := range sim.Grid.Cells {
			p, err := g.ConservativeToPrimitive(c.U)
			if err != nil {
				t.Fatalf("negative state: %v", err)
			}
			a := absF(p.U) + g.SoundSpeed(p)
			if dt*a/c.Width() > cfg.CFL+1e-12 {
				t.Fatalf("CFL violated: dt*a/dx=%.6f > %g (dx=%.4g a=%.4g)",
					dt*a/c.Width(), cfg.CFL, c.Width(), a)
			}
		}
	}
	for sim.Step < 400 {
		dt := sim.AdaptiveDT()
		if rem := cfg.TEnd - sim.Time; dt > rem {
			dt = rem
		}
		check(dt)
		if err := sim.Advance(); err != nil {
			t.Fatalf("advance: %v", err)
		}
		if sim.Time >= cfg.TEnd {
			break
		}
	}
}

// A periodic problem integrated for many acoustic transit times must remain
// finite and strictly positive: no blow-up, no negative density/pressure.
func TestLongTimeStability(t *testing.T) {
	cfg := testConfig(256, BoundaryPeriodic)
	cfg.TEnd = 5.0 // ~ many domain transits (c ~ 1)
	cfg.CFL = 0.8
	sim, err := NewSimulation(cfg)
	if err != nil {
		t.Fatal(err)
	}
	m0 := sim.TotalMass()
	if err := sim.Run(); err != nil {
		t.Fatalf("long run failed: %v", err)
	}
	g := cfg.gas()
	for i, c := range sim.Grid.Cells {
		if _, err := g.ConservativeToPrimitive(c.U); err != nil {
			t.Fatalf("negative state at cell %d: %v", i, err)
		}
	}
	if d := sim.TotalMass() - m0; absF(d) > 1e-11 {
		t.Fatalf("mass changed over long run: %.3e", d)
	}
}

// AdaptiveDT must shrink as the mesh refines and never be non-positive.
func TestAdaptiveDTPositive(t *testing.T) {
	coarse, _ := NewSimulation(testConfig(128, BoundaryReflecting))
	fine, _ := NewSimulation(testConfig(512, BoundaryReflecting))
	if dt := coarse.AdaptiveDT(); dt <= 0 {
		t.Fatalf("nonpositive dt coarse: %g", dt)
	}
	if dtF := fine.AdaptiveDT(); dtF >= coarse.AdaptiveDT() {
		t.Fatalf("finer mesh did not reduce dt: fine=%g coarse=%g", dtF, coarse.AdaptiveDT())
	}
}
