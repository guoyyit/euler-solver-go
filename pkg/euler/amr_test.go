package euler

import (
	"math"
	"testing"
)

func amrSodCfg(n, maxL int, rt, ct float64) Config {
	c := sodCfg(n)
	c.AMR = DefaultAMRConfig()
	c.AMR.MaxLevels = maxL
	c.AMR.RefineThreshold = rt
	c.AMR.CoarsenThreshold = ct
	c.AMR.MergeBufferCells = 1
	return c
}

func TestAMRConservationSod(t *testing.T) {
	p := GetProblem(ProblemSod)
	cfg := amrSodCfg(100, 2, 0.05, 0.005)
	s, err := NewSolver(cfg, p.Init)
	if err != nil {
		t.Fatal(err)
	}
	r0, m0, e0 := s.TotalConserved()
	if err := s.Advance(0.2); err != nil {
		t.Fatal(err)
	}
	if len(s.Mesh().Levels) < 2 {
		t.Fatal("expected refinement to create a fine level")
	}
	r1, m1, e1 := s.TotalConserved()
	tol := 1e-11 * math.Max(1, math.Abs(e0))
	if d := math.Abs(r1 - r0); d > tol {
		t.Fatalf("AMR mass drift %.3e", d)
	}
	if d := math.Abs(m1 - m0); d > tol {
		t.Fatalf("AMR momentum drift %.3e", d)
	}
	if d := math.Abs(e1 - e0); d > tol {
		t.Fatalf("AMR energy drift %.3e", d)
	}
	for _, lv := range s.Mesh().Levels {
		for i := NGhost; i < NGhost+lv.NCells(); i++ {
			if lv.Q[i].Rho <= 0 || lv.Q[i].P <= 0 {
				t.Fatalf("negative AMR state dx=%v: %+v", lv.DX, lv.Q[i])
			}
		}
	}
}
