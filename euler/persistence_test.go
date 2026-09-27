package euler

import (
	"os"
	"path/filepath"
	"testing"
)

// A snapshot must round-trip every bit (JSON uses shortest round-trip decimal).
func TestSnapshotRoundTrip(t *testing.T) {
	sim, err := NewSimulation(amrConfig(128, BoundaryPeriodic))
	if err != nil {
		t.Fatal(err)
	}
	for k := 0; k < 10; k++ {
		if err := sim.Advance(); err != nil {
			t.Fatal(err)
		}
	}
	snap := sim.Snapshot()
	dir := t.TempDir()
	path := filepath.Join(dir, "snap.json")
	if err := SaveSnapshot(path, snap); err != nil {
		t.Fatal(err)
	}
	got, err := LoadSnapshot(path)
	if err != nil {
		t.Fatal(err)
	}
	if got.Time != sim.Time || got.Step != sim.Step || len(got.Cells) != len(sim.Grid.Cells) {
		t.Fatal("scalar snapshot mismatch")
	}
	for i := range got.Cells {
		a, b := got.Cells[i], snap.Cells[i]
		if a != b {
			t.Fatalf("cell %d not bit-identical: %+v vs %+v", i, a, b)
		}
	}
}

// Continuing from a checkpoint must give the IDENTICAL trajectory as an
// uninterrupted run (the core crash-recovery guarantee).
func TestResumeMatchesContinuous(t *testing.T) {
	cfg := amrConfig(128, BoundaryReflecting)
	continuous, err := NewSimulation(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := continuous.Run(); err != nil {
		t.Fatal(err)
	}

	// Interrupted run: stop partway, snapshot, resume to the end.
	interrupted, err := NewSimulation(cfg)
	if err != nil {
		t.Fatal(err)
	}
	cut := int64(37)
	for interrupted.Step < cut {
		if err := interrupted.Advance(); err != nil {
			t.Fatal(err)
		}
	}
	dir := t.TempDir()
	path := filepath.Join(dir, "crash.json")
	if err := interrupted.SaveSnapshot(path); err != nil {
		t.Fatal(err)
	}
	snap, err := LoadSnapshot(path)
	if err != nil {
		t.Fatal(err)
	}
	resumed, err := Resume(snap)
	if err != nil {
		t.Fatal(err)
	}
	if err := resumed.Run(); err != nil {
		t.Fatal(err)
	}

	if resumed.Time != continuous.Time || resumed.Step != continuous.Step {
		t.Fatalf("scalar mismatch: resumed t=%v n=%v vs continuous t=%v n=%v",
			resumed.Time, resumed.Step, continuous.Time, continuous.Step)
	}
	if len(resumed.Grid.Cells) != len(continuous.Grid.Cells) {
		t.Fatalf("cell count mismatch %d vs %d",
			len(resumed.Grid.Cells), len(continuous.Grid.Cells))
	}
	for i, c := range continuous.Grid.Cells {
		r := resumed.Grid.Cells[i]
		if c.ID != r.ID || c.Level != r.Level || c.U != r.U {
			t.Fatalf("cell %d diverged after resume: %+v vs %+v", i, c.U, r.U)
		}
	}
}

// Two independent runs with identical config must be bit-for-bit identical,
// and the checkpoint file bytes must be stable across rewrites.
func TestDeterministicOutput(t *testing.T) {
	cfg := amrConfig(128, BoundaryPeriodic)
	run := func() *Simulation {
		s, err := NewSimulation(cfg)
		if err != nil {
			t.Fatal(err)
		}
		if err := s.Run(); err != nil {
			t.Fatal(err)
		}
		return s
	}
	a, b := run(), run()
	if len(a.Grid.Cells) != len(b.Grid.Cells) {
		t.Fatal("nondeterministic cell count")
	}
	for i := range a.Grid.Cells {
		if a.Grid.Cells[i].U != b.Grid.Cells[i].U {
			t.Fatalf("nondeterministic state at cell %d", i)
		}
	}
	dir := t.TempDir()
	p1 := filepath.Join(dir, "a.json")
	p2 := filepath.Join(dir, "b.json")
	if err := a.SaveSnapshot(p1); err != nil {
		t.Fatal(err)
	}
	if err := b.SaveSnapshot(p2); err != nil {
		t.Fatal(err)
	}
	d1, _ := os.ReadFile(p1)
	d2, _ := os.ReadFile(p2)
	if string(d1) != string(d2) {
		t.Fatal("snapshot file bytes differ between identical runs")
	}
}
