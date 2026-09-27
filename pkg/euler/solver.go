package euler

import (
	"fmt"
)

// Config bundles all user-facing solver options.
type Config struct {
	X0, X1   float64 // domain [X0,X1)
	NBase    int     // number of level-0 cells
	Boundary BoundaryKind
	Flux     FluxScheme
	Limiter  Limiter
	CFL      float64
	AMR      AMRConfig
	Gamma    float64
}

// Solver marches the Euler equations on an adaptive hierarchy.
type Solver struct {
	cfg  Config
	eos  *EOS
	mesh *Mesh
	step int64
	t    float64

	// CheckpointFunc, if set, is invoked by StepOne after every completed time
	// step with a serialisable snapshot. A solver restored from such a snapshot
	// and resumed reproduces the run bit-for-bit because checkpoints are only
	// written at step boundaries and each dt is a deterministic function of
	// the restored state.
	CheckpointFunc func(data []byte, t float64, step int64)
}

// NewSolver constructs a solver and samples the initial condition at cell
// centres. Every initial cell must be a physical state.
func NewSolver(cfg Config, ic func(x float64) Primitive) (*Solver, error) {
	if err := cfg.validate(); err != nil {
		return nil, err
	}
	eos := NewEOS(cfg.Gamma)
	cfg.AMR = cfg.AMR.withDefaults()
	m := &Mesh{
		Ratio: cfg.AMR.Ratio,
		X0:    cfg.X0,
		X1:    cfg.X1,
		NBase: cfg.NBase,
	}
	lv := NewLevel(0, cfg.NBase, m.BaseDX())
	for i := 0; i < cfg.NBase; i++ {
		x := m.CellCenter(i, m.BaseDX())
		q := ic(x)
		if !(q.Rho > 0) || !(q.P > 0) {
			return nil, fmt.Errorf("non-physical initial state at x=%v: %+v", x, q)
		}
		lv.U[NGhost+i] = eos.PrimToCons(q)
		lv.Q[NGhost+i] = q
	}
	FillGhostCons(lv.U, NGhost, NGhost+cfg.NBase, NGhost, cfg.Boundary)
	FillGhost(lv.Q, NGhost, NGhost+cfg.NBase, NGhost, cfg.Boundary)
	m.Levels = []*Level{lv}
	sol := &Solver{cfg: cfg, eos: eos, mesh: m}
	if cfg.AMR.MaxLevels > 1 {
		if err := sol.adapt(); err != nil {
			return nil, err
		}
	}
	return sol, nil
}

func (c Config) validate() error {
	if !(c.X1 > c.X0) {
		return fmt.Errorf("empty domain: X0=%v X1=%v", c.X0, c.X1)
	}
	if c.NBase < 4 {
		return fmt.Errorf("NBase must be >= 4, got %d", c.NBase)
	}
	if c.Gamma <= 1 {
		return fmt.Errorf("gamma must be > 1, got %v", c.Gamma)
	}
	if c.CFL <= 0 || c.CFL > 1 {
		return fmt.Errorf("CFL must be in (0,1], got %v", c.CFL)
	}
	return nil
}

// Time returns the current simulation time.
func (s *Solver) Time() float64 { return s.t }

// Step returns the number of completed time steps.
func (s *Solver) Step() int64 { return s.step }

// Mesh returns the current hierarchy (read-only use outside the package).
func (s *Solver) Mesh() *Mesh { return s.mesh }

// EOS returns the equation-of-state object.
func (s *Solver) EOS() *EOS { return s.eos }

func (s *Solver) stepConfig() StepConfig {
	return StepConfig{Flux: s.cfg.Flux, Limiter: s.cfg.Limiter, Boundary: s.cfg.Boundary}
}

// Advance integrates to exactly target. Each step size is computed
// deterministically from the current state:
//
//	dt = min(safety*cfl*dx/smax, target-t)
//
// Because the state at a given time is byte-identical for a run resumed from a
// snapshot, the recomputed dt and the resulting step sequence are identical to
// an uninterrupted run, while the safety factor guarantees the step stays
// CFL-legal even as shock wave speeds evolve. A final short step lands exactly
// on target.
func (s *Solver) Advance(target float64) error {
	for s.t < target {
		done, err := s.StepOne(target)
		if err != nil {
			return err
		}
		if done {
			break
		}
	}
	return nil
}

// NextDT returns the CFL-limited step that the next time step would use
// (before any target capping). It is deterministic in the current state.
func (s *Solver) NextDT() float64 {
	return StableDT(s.mesh, s.eos, s.cfg.CFL) * dtSafety
}

// StepOne advances by exactly one physical time step of length min(NextDT,
// target-Time). It is the unit at which checkpoints are valid: snapshots may
// only be taken (and resumed) after StepOne returns.
func (s *Solver) StepOne(target float64) (bool, error) {
	if s.t >= target {
		return true, nil
	}
	dt := s.NextDT()
	if s.t+dt > target {
		dt = target - s.t
	}
	if err := s.takeStep(dt, 0); err != nil {
		return false, err
	}
	s.t += dt
	s.step++
	if s.cfg.AMR.MaxLevels > 1 {
		if err := s.adapt(); err != nil {
			return false, err
		}
	}
	if s.CheckpointFunc != nil {
		data, err := s.Snapshot()
		if err != nil {
			return false, err
		}
		s.CheckpointFunc(data, s.t, s.step)
	}
	return s.t >= target, nil
}

// dtSafety is the uniform margin below the instantaneous CFL bound; 0.9 gives
// robust headroom for wave-speed growth during shock formation.
const dtSafety = 0.9

// takeStep advances the hierarchy by dt. With the CFL safety factor the
// nominal step is already conservative; the half-step retry is a defensive
// last resort against roundoff-level positivity violations.
func (s *Solver) takeStep(dt float64, depth int) error {
	const maxRetry = 20
	if err := s.advanceHierarchy(dt); err != nil {
		if depth < maxRetry {
			return s.takeStep(0.5*dt, depth+1)
		}
		return fmt.Errorf("step failed at t=%v even after %d retries: %w", s.t, maxRetry, err)
	}
	return nil
}

// advanceHierarchy advances the single active level (AMR path replaced in
// amr.go integration) while keeping every level time-consistent.
func (s *Solver) advanceHierarchy(dt float64) error {
	if len(s.mesh.Levels) > 1 {
		if err := s.coupledAdvance(dt); err != nil {
			return err
		}
		for _, lv := range s.mesh.Levels {
			lv.Time += dt
		}
		return nil
	}
	lv := s.mesh.Levels[0]
	if err := advanceLevel(lv, s.eos, s.stepConfig(), dt); err != nil {
		return err
	}
	lv.Time += dt
	return nil
}

// TotalConserved returns the global integrals of mass, momentum and energy,
// counting every spatial point exactly once across the hierarchy.
func (s *Solver) TotalConserved() (r, mm, e float64) {
	for _, part := range s.uniqueCells() {
		dx := part.lv.DX
		for _, p := range part.idx {
			u := part.lv.U[p]
			r += u.R * dx
			mm += u.M * dx
			e += u.E * dx
		}
	}
	return r, mm, e
}

type cellPart struct {
	lv  *Level
	idx []int
}

// uniqueCells partitions the domain into packed cell indices counted exactly
// once: fine cells win over coarser regions they cover.
func (s *Solver) uniqueCells() []cellPart {
	if len(s.mesh.Levels) == 1 {
		lv := s.mesh.Levels[0]
		idx := make([]int, lv.NCells())
		for i := range idx {
			idx[i] = NGhost + i
		}
		return []cellPart{{lv: lv, idx: idx}}
	}
	return s.uniqueCellsAMR()
}
