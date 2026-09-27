package euler

import "testing"

// TestAMRImprovesOverCoarse verifies that refining around the waves lowers the
// density L1 error relative to a high-resolution reference, using the same base
// cell count. The composite profile is sampled at base resolution with fine
// cells averaged back into their parents.
func TestAMRImprovesOverCoarse(t *testing.T) {
	p := GetProblem(ProblemSod)
	ref := profileDensity(t, amrSodCfg(800, 1, 0, 0), 0.2, p.Init) // fine reference
	errOf := func(cfg Config) float64 {
		q := profileDensity(t, cfg, 0.2, p.Init)
		factor := 800 / cfg.NBase
		var sum float64
		for i := 0; i < cfg.NBase; i++ {
			var ravg float64
			for k := 0; k < factor; k++ {
				ravg += ref[i*factor+k]
			}
			ravg /= float64(factor)
			d := q[i] - ravg
			if d < 0 {
				d = -d
			}
			sum += d / float64(cfg.NBase)
		}
		return sum
	}
	coarse := errOf(amrSodCfg(200, 1, 0, 0))
	amr := errOf(amrSodCfg(200, 2, 0.05, 0.005))
	if !(amr < coarse) {
		t.Fatalf("AMR did not reduce error: coarse=%v amr=%v", coarse, amr)
	}
}

// profileDensity returns the base-grid density, using fine-cell averages in
// refined regions (a conservative down-sampling).
func profileDensity(t *testing.T, cfg Config, T float64, ic func(x float64) Primitive) []float64 {
	t.Helper()
	s, err := NewSolver(cfg, ic)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Advance(T); err != nil {
		t.Fatal(err)
	}
	n := cfg.NBase
	out := make([]float64, n)
	base := s.Mesh().Levels[0]
	for i := 0; i < n; i++ {
		out[i] = base.Q[NGhost+i].Rho
	}
	if len(s.Mesh().Levels) > 1 {
		fine := s.Mesh().Levels[1]
		r := cfg.AMR.Ratio
		for _, iv := range fine.Active {
			for fi := iv.Lo; fi < iv.Hi; fi += r {
				var rho float64
				for k := 0; k < r; k++ {
					rho += fine.Q[fine.at(fi+k)].Rho
				}
				out[fi/r] = rho / float64(r)
			}
		}
	}
	return out
}
