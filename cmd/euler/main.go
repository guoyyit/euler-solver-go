package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"os"
	"strconv"

	"github.com/guoyyit/euler-solver-go/euler"
)

// sodIC is the classic Sod shock tube:
//
//	(rho,u,p)_L = (1.0,0,1.0), (rho,u,p)_R = (0.125,0,0.1), x0=0.5.
func sodIC(x float64) euler.Primitive {
	if x < 0.5 {
		return euler.Primitive{Rho: 1.0, U: 0, P: 1.0}
	}
	return euler.Primitive{Rho: 0.125, U: 0, P: 0.1}
}

func main() {
	var (
		ic        = flag.String("ic", "sod", "named initial condition (built-in: sod)")
		xmin      = flag.Float64("xmin", 0, "domain left end")
		xmax      = flag.Float64("xmax", 1, "domain right end")
		n0        = flag.Int("n0", 256, "initial number of cells (power of two)")
		gamma     = flag.Float64("gamma", 1.4, "ratio of specific heats")
		boundary  = flag.String("boundary", "reflecting", "periodic|reflecting")
		flux      = flag.String("flux", "rusanov", "rusanov|hllc")
		limiter   = flag.String("limiter", "minmod", "minmod|vanleer")
		cfl       = flag.Float64("cfl", 0.8, "CFL number in (0,1]")
		tend      = flag.Float64("tend", 0.2, "final physical time")
		maxSteps  = flag.Int64("maxsteps", 0, "hard cap on steps (0 = unlimited)")
		amr       = flag.Bool("amr", true, "enable adaptive refinement")
		minLevel  = flag.Int("minlevel", 0, "coarsest allowed level (relative)")
		maxLevel  = flag.Int("maxlevel", 3, "finest allowed level (relative)")
		refineAt  = flag.Float64("refine-at", 0.05, "score above which a cell refines")
		coarsenAt = flag.Float64("coarsen-at", 0.01, "score below which calm siblings merge")
		criterion = flag.String("refine-criterion", "density-gradient",
			"density-gradient|pressure-gradient")
		regridEvery = flag.Int("regrid-every", 4, "regrid every N steps (0 disables)")
		snapEvery   = flag.Int("snapshot-every", 0, "checkpoint every N steps (0 = end only)")
		snapPath    = flag.String("snapshot", "snapshot.json", "checkpoint file path")
		noSnap      = flag.Bool("no-snapshot", false, "do not write a checkpoint")
		resume      = flag.String("resume", "", "resume from this checkpoint instead of starting fresh")
		out         = flag.String("out", "solution.csv", "final solution CSV (empty disables)")
	)
	flag.Parse()

	euler.RegisterProblem(euler.Problem{Name: "sod", IC: sodIC})

	cfg := euler.Config{
		ICName:        *ic,
		Xmin:          *xmin,
		Xmax:          *xmax,
		N0:            *n0,
		Gamma:         *gamma,
		Boundary:      euler.BoundaryType(*boundary),
		Flux:          euler.FluxScheme(*flux),
		Limiter:       euler.Limiter(*limiter),
		CFL:           *cfl,
		TEnd:          *tend,
		MaxSteps:      *maxSteps,
		RegridEvery:   *regridEvery,
		SnapshotEvery: *snapEvery,
		OutputPath:    *out,
		Refine: euler.RefineConfig{
			Enabled:   *amr,
			MinLevel:  *minLevel,
			MaxLevel:  *maxLevel,
			RefineAt:  *refineAt,
			CoarsenAt: *coarsenAt,
			Criterion: euler.RefineCriterion(*criterion),
		},
	}
	if !*noSnap {
		cfg.SnapshotPath = *snapPath
	}

	var sim *euler.Simulation
	if *resume != "" {
		snap, err := euler.LoadSnapshot(*resume)
		if err != nil {
			fatal("load snapshot: %v", err)
		}
		sim, err = euler.Resume(snap)
		if err != nil {
			fatal("resume: %v", err)
		}
		fmt.Printf("resumed at t=%.12g step=%d cells=%d\n", sim.Time, sim.Step, len(sim.Grid.Cells))
		// Output destinations belong to this invocation; override the paths
		// carried inside the checkpoint without touching physics settings.
		sim.Cfg.SnapshotPath = cfg.SnapshotPath
		sim.Cfg.SnapshotEvery = cfg.SnapshotEvery
		sim.Cfg.OutputPath = cfg.OutputPath
	} else {
		var err error
		sim, err = euler.NewSimulation(cfg)
		if err != nil {
			fatal("init: %v", err)
		}
	}

	m0, p0, e0 := sim.TotalMass(), sim.TotalMomentum(), sim.TotalEnergy()
	if err := sim.Run(); err != nil {
		fatal("run: %v", err)
	}
	m1, p1, e1 := sim.TotalMass(), sim.TotalMomentum(), sim.TotalEnergy()

	fmt.Printf("done t=%.12g steps=%d cells=%d\n", sim.Time, sim.Step, len(sim.Grid.Cells))
	fmt.Printf("conservation drift |dM|=%.3e |dP|=%.3e |dE|=%.3e\n",
		absD(m1-m0), absD(p1-p0), absD(e1-e0))

	if *out != "" {
		if err := writeCSV(*out, sim); err != nil {
			fatal("write output: %v", err)
		}
		fmt.Printf("wrote %s\n", *out)
	}
}

func writeCSV(path string, sim *euler.Simulation) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	defer w.Flush()
	if err := w.Write([]string{"x", "dx", "level", "rho", "u", "p", "rhou", "E"}); err != nil {
		return err
	}
	for _, c := range sim.Grid.Cells {
		pr, err := sim.Cfg.Gas().ConservativeToPrimitive(c.U)
		if err != nil {
			return err
		}
		rec := []string{
			strconv.FormatFloat(c.Center(), 'g', 17, 64),
			strconv.FormatFloat(c.Width(), 'g', 17, 64),
			strconv.Itoa(c.Level),
			strconv.FormatFloat(pr.Rho, 'g', 17, 64),
			strconv.FormatFloat(pr.U, 'g', 17, 64),
			strconv.FormatFloat(pr.P, 'g', 17, 64),
			strconv.FormatFloat(c.U.RhoU, 'g', 17, 64),
			strconv.FormatFloat(c.U.E, 'g', 17, 64),
		}
		if err := w.Write(rec); err != nil {
			return err
		}
	}
	return w.Error()
}

func absD(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func fatal(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "euler: "+format+"\n", args...)
	os.Exit(1)
}
