// Command euler1d runs a 1-D compressible Euler simulation with shock
// capturing, adaptive mesh refinement, deterministic time stepping and
// crash-recovery snapshots.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"

	"github.com/guoyyit/euler-solver-go/pkg/euler"
)

func main() {
	var (
		problem     = flag.String("problem", "sod", "built-in problem: sod|lax|2raref")
		n           = flag.Int("n", 200, "number of base cells")
		cfl         = flag.Float64("cfl", 0.8, "CFL number (0,1]")
		tEnd        = flag.Float64("t", 0.2, "final simulation time")
		flux        = flag.String("flux", "hll", "numerical flux: rusanov|hll")
		limiter     = flag.String("limiter", "minmod", "slope limiter: none|minmod|mc")
		bc          = flag.String("bc", "", "boundary: periodic|reflect (default: problem)")
		maxLevels   = flag.Int("levels", 1, "number of AMR levels (1 = uniform)")
		refine      = flag.Float64("refine", 0.05, "refinement threshold")
		coarsen     = flag.Float64("coarsen", 0.005, "coarsening threshold")
		snapshotIn  = flag.String("resume", "", "restart from this snapshot file")
		snapshotOut = flag.String("snapshot", "", "write final snapshot to this file")
		every       = flag.Int("checkpoint-every", 0, "if >0, write <snapshot>.step after each N steps")
		out         = flag.String("out", "", "write final primitive profile (x rho u p per line)")
	)
	flag.Parse()

	p := parseProblem(*problem)
	boundary := p.Boundary
	if *bc == "periodic" {
		boundary = euler.BoundaryPeriodic
	} else if *bc == "reflect" {
		boundary = euler.BoundaryReflecting
	}

	cfg := euler.Config{
		X0:       0,
		X1:       1,
		NBase:    *n,
		Boundary: boundary,
		Flux:     parseFlux(*flux),
		Limiter:  parseLimiter(*limiter),
		CFL:      *cfl,
		Gamma:    p.Gamma,
		AMR: euler.AMRConfig{
			Ratio:            2,
			MaxLevels:        *maxLevels,
			NGhost:           euler.NGhost,
			MergeBufferCells: 1,
			RefineThreshold:  *refine,
			CoarsenThreshold: *coarsen,
			CriterionEps:     1e-12,
		},
	}

	var s *euler.Solver
	if *snapshotIn != "" {
		data, err := os.ReadFile(*snapshotIn)
		if err != nil {
			log.Fatalf("read snapshot: %v", err)
		}
		s, err = euler.Restore(data)
		if err != nil {
			log.Fatalf("restore: %v", err)
		}
	} else {
		var err error
		s, err = euler.NewSolver(cfg, p.Init)
		if err != nil {
			log.Fatalf("init: %v", err)
		}
	}

	if *every > 0 && *snapshotOut != "" {
		count := 0
		s.CheckpointFunc = func(data []byte, t float64, step int64) {
			count++
			if count%*every == 0 {
				name := *snapshotOut + "." + strconv.FormatInt(step, 10)
				if err := os.WriteFile(name, data, 0o644); err != nil {
					log.Printf("checkpoint: %v", err)
				}
			}
		}
	}

	if err := s.Advance(*tEnd); err != nil {
		log.Fatalf("advance: %v", err)
	}

	if *snapshotOut != "" {
		if err := writeSnapshot(*snapshotOut, s); err != nil {
			log.Fatal(err)
		}
	}
	if *out != "" {
		if err := writeProfile(*out, s); err != nil {
			log.Fatal(err)
		}
	}

	r, m, e := s.TotalConserved()
	fmt.Printf("done t=%.6f steps=%d levels=%d total(rho,rho*u,E)=%.12g %.12g %.12g\n",
		s.Time(), s.Step(), len(s.Mesh().Levels), r, m, e)
}

func parseProblem(name string) euler.Problem {
	switch name {
	case "lax":
		return euler.GetProblem(euler.ProblemLax)
	case "2raref", "double-rarefaction":
		return euler.GetProblem(euler.ProblemDoubleRarefaction)
	default:
		return euler.GetProblem(euler.ProblemSod)
	}
}

func parseFlux(name string) euler.FluxScheme {
	if name == "rusanov" {
		return euler.FluxRusanov
	}
	return euler.FluxHLL
}

func parseLimiter(name string) euler.Limiter {
	switch name {
	case "none":
		return euler.LimitNone
	case "mc":
		return euler.LimitMC
	default:
		return euler.LimitMinmod
	}
}

func writeSnapshot(path string, s *euler.Solver) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	return euler.WriteSnapshot(f, s)
}

// writeProfile writes the composite (fine-preferred) primitive solution sampled
// on the base grid for easy inspection.
func writeProfile(path string, s *euler.Solver) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	mesh := s.Mesh()
	base := mesh.Levels[0]
	for i := 0; i < base.NCells(); i++ {
		x := mesh.CellCenter(i, base.DX)
		q := base.Q[euler.NGhost+i]
		fmt.Fprintf(f, "%.10g %.10g %.10g %.10g\n", x, q.Rho, q.U, q.P)
	}
	return nil
}
