package euler

import (
	"bytes"
	"math"
	"testing"
)

func compareSolvers(t *testing.T, a, b *Solver) {
	t.Helper()
	if len(a.Mesh().Levels) != len(b.Mesh().Levels) {
		t.Fatalf("level count differs: %d vs %d", len(a.Mesh().Levels), len(b.Mesh().Levels))
	}
	for l, lv := range a.Mesh().Levels {
		lv2 := b.Mesh().Levels[l]
		if lv.IV != lv2.IV || len(lv.U) != len(lv2.U) {
			t.Fatalf("level %d geometry differs: %v vs %v", l, lv.IV, lv2.IV)
		}
		if len(lv.Active) != len(lv2.Active) {
			t.Fatalf("level %d active intervals differ", l)
		}
		for i := range lv.Active {
			if lv.Active[i] != lv2.Active[i] {
				t.Fatalf("level %d active %v differs", l, lv.Active)
			}
		}
		for i := range lv.U {
			if lv.U[i] != lv2.U[i] {
				t.Fatalf("level %d cell %d differs: %v vs %v", l, i, lv.U[i], lv2.U[i])
			}
		}
	}
	if a.Time() != b.Time() || a.Step() != b.Step() {
		t.Fatalf("runtime differs: t %v/%v step %d/%d",
			a.Time(), b.Time(), a.Step(), b.Step())
	}
}

// Crash-recovery contract: checkpoints are written after completed steps.
// We run the uninterrupted reference to tEnd, then run a second instance that
// crashes after a fixed number of steps, restores from the latest checkpoint,
// and continues; the final states must be bit-identical.
func TestSnapshotRestartBitExact(t *testing.T) {
	p := GetProblem(ProblemSod)
	cfg := amrSodCfg(120, 2, 0.05, 0.005)
	const tEnd = 0.2
	const crashAfter = 37

	cont, err := NewSolver(cfg, p.Init)
	if err != nil {
		t.Fatal(err)
	}
	for {
		done, err := cont.StepOne(tEnd)
		if err != nil {
			t.Fatal(err)
		}
		if done {
			break
		}
	}

	crashed, err := NewSolver(cfg, p.Init)
	if err != nil {
		t.Fatal(err)
	}
	var latest []byte
	steps := 0
	for {
		done, err := crashed.StepOne(tEnd)
		if err != nil {
			t.Fatal(err)
		}
		if done {
			break
		}
		if latest == nil {
			latest, _ = crashed.Snapshot()
		}
		steps++
		if steps == crashAfter {
			// Simulate a crash: discard the running instance and restore.
			break
		}
		latest, _ = crashed.Snapshot()
	}
	r, err := Restore(latest)
	if err != nil {
		t.Fatal(err)
	}
	for {
		done, err := r.StepOne(tEnd)
		if err != nil {
			t.Fatal(err)
		}
		if done {
			break
		}
	}
	compareSolvers(t, cont, r)
}

func TestSnapshotRoundTripAndDeterminism(t *testing.T) {
	p := GetProblem(ProblemLax)
	cfg := sodCfg(96)
	s, err := NewSolver(cfg, p.Init)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Advance(0.1); err != nil {
		t.Fatal(err)
	}
	d1, err := s.Snapshot()
	if err != nil {
		t.Fatal(err)
	}
	d2, err := s.Snapshot()
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(d1, d2) {
		t.Fatal("snapshot bytes not stable")
	}
	r, err := Restore(d1)
	if err != nil {
		t.Fatal(err)
	}
	d3, _ := r.Snapshot()
	if !bytes.Equal(d1, d3) {
		t.Fatal("snapshot changes on restore round-trip")
	}
}

func TestSnapshotFileIO(t *testing.T) {
	p := GetProblem(ProblemSod)
	cfg := amrSodCfg(80, 2, 0.06, 0.01)
	s, err := NewSolver(cfg, p.Init)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Advance(0.12); err != nil {
		t.Fatal(err)
	}
	var buf bytes.Buffer
	if err := WriteSnapshot(&buf, s); err != nil {
		t.Fatal(err)
	}
	r, err := ReadSnapshot(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatal(err)
	}
	compareSolvers(t, s, r)
	a1, b1, c1 := s.TotalConserved()
	a2, b2, c2 := r.TotalConserved()
	if math.Abs(a1-a2) > 1e-14 || math.Abs(b1-b2) > 1e-14 || math.Abs(c1-c2) > 1e-14 {
		t.Fatalf("totals changed: %.17g vs %.17g", a1, a2)
	}
	// Restored solver must continue identically.
	if err := r.Advance(0.2); err != nil {
		t.Fatal(err)
	}
	s2, _ := Restore(append([]byte(nil), buf.Bytes()...))
	if err := s2.Advance(0.2); err != nil {
		t.Fatal(err)
	}
	if math.Abs(r.Time()-s2.Time()) > 0 {
		t.Fatal("restart time differs")
	}
}
