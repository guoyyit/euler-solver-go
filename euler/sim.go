package euler

import (
	"fmt"
	"math"
	"sync"
)

// InitialCondition prescribes primitive data as a function of position.
type InitialCondition func(x float64) Primitive

// Config is the complete configuration of a run. Identical Config together
// with identical initial data produces bit-identical output (determinism).
type Config struct {
	ICName        string
	Xmin          float64
	Xmax          float64
	N0            int // number of level-0 cells covering the domain
	Gamma         float64
	Boundary      BoundaryType
	Flux          FluxScheme
	Limiter       Limiter
	CFL           float64
	TEnd          float64
	MaxSteps      int64
	Refine        RefineConfig
	RegridEvery   int    // 0 disables regridding after initialization
	SnapshotEvery int    // take a snapshot every N steps (0 = at end only)
	SnapshotPath  string // checkpoint file
	OutputPath    string // final state CSV
}

func (c Config) gas() Gas { return Gas{Gamma: c.Gamma} }

// Gas returns the configured equation-of-state handle.
func (c Config) Gas() Gas { return c.gas() }

func (c Config) periodic() bool { return c.Boundary == BoundaryPeriodic }

// Validate checks the configuration for internal consistency.
func (c Config) Validate() error {
	if !(c.Xmax > c.Xmin) {
		return fmt.Errorf("Xmax must exceed Xmin")
	}
	if c.N0 < 1 {
		return fmt.Errorf("N0 must be >= 1")
	}
	if c.N0 > (1<<RootLevel) || (c.N0&(c.N0-1)) != 0 {
		return fmt.Errorf("N0 must be a power of two and <= %d, got %d", 1<<RootLevel, c.N0)
	}
	if _, err := NewGas(c.Gamma); err != nil {
		return err
	}
	if c.Boundary != BoundaryPeriodic && c.Boundary != BoundaryReflecting {
		return fmt.Errorf("unknown boundary type %q", c.Boundary)
	}
	if !(c.CFL > 0) || c.CFL > 1 {
		return fmt.Errorf("CFL must be in (0,1], got %g", c.CFL)
	}
	if c.Refine.Enabled {
		if c.Refine.MaxLevel < c.Refine.MinLevel {
			return fmt.Errorf("Refine.MaxLevel < MinLevel")
		}
		if !(c.Refine.RefineAt >= 0) || c.Refine.RefineAt > 1 {
			return fmt.Errorf("RefineAt must be in [0,1]")
		}
	}
	return nil
}

// Problem is a named built-in initial condition.
type Problem struct {
	Name string
	IC   InitialCondition
}

var (
	problemMu sync.RWMutex
	problems  = map[string]Problem{}
)

// RegisterProblem adds a named initial condition.
func RegisterProblem(p Problem) {
	problemMu.Lock()
	defer problemMu.Unlock()
	problems[p.Name] = p
}

// GetProblem returns a registered named problem.
func GetProblem(name string) (Problem, error) {
	problemMu.RLock()
	defer problemMu.RUnlock()
	if p, ok := problems[name]; ok {
		return p, nil
	}
	return Problem{}, fmt.Errorf("unknown problem %q", name)
}

// Simulation is a running state machine: grid, physical time, step counter.
type Simulation struct {
	Cfg  Config
	Grid *Grid
	Time float64
	Step int64
}

// NewSimulation builds the level-0 tessellation (N0 cells), initializes each
// cell by averaging the initial condition with a two-point Gauss rule, and
// performs the first adaptive regrid when refinement is enabled.
func NewSimulation(cfg Config) (*Simulation, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	prob, err := GetProblem(cfg.ICName)
	if err != nil {
		return nil, err
	}
	g := cfg.gas()
	m := NewGrid(cfg.Xmin, cfg.Xmax)
	w := (cfg.Xmax - cfg.Xmin) / float64(cfg.N0)
	for k := 0; k < cfg.N0; k++ {
		id := BaseID(k)
		x0 := cfg.Xmin + float64(k)*w
		x1 := x0 + w
		u, err := cellAverage(g, prob.IC, x0, x1)
		if err != nil {
			return nil, err
		}
		m.Cells = append(m.Cells, &Cell{ID: id, Level: RootLevel, X0: x0, X1: x1, U: u})
	}
	m.sortCells()
	s := &Simulation{Cfg: cfg, Grid: m}
	if cfg.Refine.Enabled {
		if err := s.regrid(); err != nil {
			return nil, err
		}
	}
	return s, nil
}

// cellAverage integrates the primitive IC over [x0,x1] with a two-point Gauss
// rule, converts to conserved at each node, then averages the conserved
// samples. The result is positive by construction for a positive IC.
func cellAverage(g Gas, ic InitialCondition, x0, x1 float64) (Conserved, error) {
	const s = 0.2886751345948129 // 1/sqrt(12)
	mid := 0.5 * (x0 + x1)
	h := x1 - x0
	qa := ic(mid - s*h)
	qb := ic(mid + s*h)
	ca, err := g.PrimitiveToConservative(qa)
	if err != nil {
		return Conserved{}, err
	}
	cb, err := g.PrimitiveToConservative(qb)
	if err != nil {
		return Conserved{}, err
	}
	return Conserved{
		Rho:  0.5 * (ca.Rho + cb.Rho),
		RhoU: 0.5 * (ca.RhoU + cb.RhoU),
		E:    0.5 * (ca.E + cb.E),
	}, nil
}

// Resume rebuilds a simulation from a checkpoint. All physics/numerical
// configuration (mesh, flux, limiter, CFL, refinement) is taken verbatim from
// the snapshot so the continued trajectory is bit-identical to an
// uninterrupted run. The transient step cap MaxSteps is reset to zero
// (unlimited): a fresh process resuming after a crash should continue to the
// physical end time TEnd, not inherit a "stop after N steps" instruction that
// belonged to the interrupted invocation.
func Resume(snap *Snapshot) (*Simulation, error) {
	cfg := snap.Config
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	cfg.MaxSteps = 0
	m := NewGrid(cfg.Xmin, cfg.Xmax)
	for _, sc := range snap.Cells {
		m.Cells = append(m.Cells, &Cell{
			ID: sc.ID, Level: sc.Level, X0: sc.X0, X1: sc.X1,
			U: Conserved{Rho: sc.Rho, RhoU: sc.RhoU, E: sc.E},
		})
	}
	m.sortCells()
	return &Simulation{Cfg: cfg, Grid: m, Time: snap.Time, Step: snap.Step}, nil
}

func (s *Simulation) regrid() error {
	return Regrid(s.Cfg.gas(), s.Grid, s.Cfg.periodic(), s.Cfg.Refine)
}

// MaxWavespeed returns max_i (|u_i| + c_i), used by the CFL estimate.
func (s *Simulation) MaxWavespeed() float64 {
	g := s.Cfg.gas()
	smax := 0.0
	for _, c := range s.Grid.Cells {
		p, err := g.ConservativeToPrimitive(c.U)
		if err != nil {
			continue
		}
		v := math.Abs(p.U) + g.SoundSpeed(p)
		if v > smax {
			smax = v
		}
	}
	return smax
}

// AdaptiveDT returns the largest dt satisfying the CFL condition on every
// cell: dt = CFL * min_i dx_i / (|u_i|+c_i). The minimum over the full
// ordered cell slice is deterministic.
func (s *Simulation) AdaptiveDT() float64 {
	g := s.Cfg.gas()
	dt := math.Inf(1)
	for _, c := range s.Grid.Cells {
		p, err := g.ConservativeToPrimitive(c.U)
		if err != nil {
			continue
		}
		a := math.Abs(p.U) + g.SoundSpeed(p)
		if a <= 0 {
			continue
		}
		d := s.Cfg.CFL * c.Width() / a
		if d < dt {
			dt = d
		}
	}
	if math.IsInf(dt, 1) {
		dt = s.Cfg.CFL * (s.Cfg.Xmax - s.Cfg.Xmin)
	}
	return dt
}

// rightHandSide returns dU/dt for every leaf cell from interface fluxes.
//
// At an interface the two traces are obtained by MUSCL reconstruction when the
// adjacent cells and their outer stencil neighbors all share one width; at a
// 2:1 jump (and at boundaries) a robust first-order trace is used. This keeps
// hanging-node interfaces conservative and oscillation free.
func (s *Simulation) rightHandSide(dt float64) ([]Conserved, error) {
	g := s.Cfg.gas()
	cells := s.Grid.Cells
	n := len(cells)

	prim := make([]Primitive, n)
	for i, c := range cells {
		p, err := g.ConservativeToPrimitive(c.U)
		if err != nil {
			return nil, fmt.Errorf("step %d cell %d: %w", s.Step, i, err)
		}
		prim[i] = p
	}

	// Every interface -- including the two domain boundaries -- is treated
	// identically. owner(j) is the primitive state of cell j, with out-of-range
	// indices mapped according to the boundary type: periodic indices wrap to
	// the opposite edge, while reflecting indices return the velocity-mirrored
	// edge state (rho,-u,p). width(j) likewise provides the cell width seen at
	// an index (mirror ghosts inherit the edge width).
	//
	// A face gets a MUSCL trace from each owning cell only when that cell's
	// three-point stencil is geometrically uniform; otherwise the plain cell
	// average is used. Reflecting walls are forced first order so the mirror
	// pair presents an exactly symmetric Riemann problem with zero mass and
	// energy flux. On a periodic mesh the left face of cell 0 and the right
	// face of cell n-1 are the identical Riemann problem, hence identical
	// flux: boundary fluxes telescope and global conservation is machine
	// exact even across the domain cut.
	owner := func(j int) Primitive {
		switch {
		case j >= 0 && j < n:
			return prim[j]
		case s.Cfg.periodic():
			return prim[(j+n)%n]
		case j < 0:
			return GhostState(s.Cfg.Boundary, -1, prim[0])
		default:
			return GhostState(s.Cfg.Boundary, +1, prim[n-1])
		}
	}
	width := func(j int) float64 {
		switch {
		case j >= 0 && j < n:
			return cells[j].Width()
		case s.Cfg.periodic():
			return cells[(j+n)%n].Width()
		case j < 0:
			return cells[0].Width()
		default:
			return cells[n-1].Width()
		}
	}
	uniform := func(j int) bool { return width(j-1) == width(j) && width(j+1) == width(j) }

	// faceFlux[k] is the flux on the LEFT face of cell k (face x = x0_k).
	faceFlux := make([]Conserved, n+1)
	for face := 0; face <= n; face++ {
		li, ri := face-1, face // left-owner and right-owner cell indices
		wall := !s.Cfg.periodic() && (li < 0 || ri >= n)
		// A face is "hanging" when its two owner cells have different widths.
		// Such a face is evaluated first order from the two cell averages: this
		// guarantees the single shared face carries ONE common flux to both
		// owners (exact conservation across the coarse/fine boundary) and is
		// monotone (no reconstruction across a size jump -> no oscillation).
		hanging := width(li) != width(ri)
		var qL, qR Primitive
		if wall || hanging || !uniform(li) {
			qL = owner(li)
		} else {
			qL = musclRight(s.Cfg.Limiter, owner(li-1), owner(li), owner(li+1))
		}
		if wall || hanging || !uniform(ri) {
			qR = owner(ri)
		} else {
			qR = musclLeft(s.Cfg.Limiter, owner(ri-1), owner(ri), owner(ri+1))
		}
		faceFlux[face] = NumericalFlux(g, s.Cfg.Flux, qL, qR)
	}

	rhs := make([]Conserved, n)
	for i, c := range cells {
		fR := faceFlux[i+1]
		fL := faceFlux[i]
		inv := 1 / c.Width()
		rhs[i] = Conserved{
			Rho:  (fL.Rho - fR.Rho) * inv,
			RhoU: (fL.RhoU - fR.RhoU) * inv,
			E:    (fL.E - fR.E) * inv,
		}
	}
	return rhs, nil
}

// musclRight reconstructs the right-face trace c+0.5 slope from (l,c,r).
func musclRight(lim Limiter, l, c, r Primitive) Primitive {
	return Primitive{
		Rho: c.Rho + 0.5*Slope(lim, l.Rho, c.Rho, r.Rho),
		U:   c.U + 0.5*Slope(lim, l.U, c.U, r.U),
		P:   c.P + 0.5*Slope(lim, l.P, c.P, r.P),
	}
}

// musclLeft reconstructs the left-face trace c-0.5 slope from (l,c,r).
func musclLeft(lim Limiter, l, c, r Primitive) Primitive {
	return Primitive{
		Rho: c.Rho - 0.5*Slope(lim, l.Rho, c.Rho, r.Rho),
		U:   c.U - 0.5*Slope(lim, l.U, c.U, r.U),
		P:   c.P - 0.5*Slope(lim, l.P, c.P, r.P),
	}
}

// Advance performs one adaptive time step. It optionally regrids, computes a
// CFL-limited dt (shortened to land exactly on TEnd), and updates with the
// positivity-preserving SSP-RK2 (Heun) method.
func (s *Simulation) Advance() error {
	if s.Cfg.RegridEvery > 0 && s.Step > 0 && s.Step%int64(s.Cfg.RegridEvery) == 0 {
		if err := s.regrid(); err != nil {
			return err
		}
	}
	g := s.Cfg.gas()
	dt := s.AdaptiveDT()
	if s.Cfg.TEnd > s.Time {
		if rem := s.Cfg.TEnd - s.Time; dt > rem {
			dt = rem
		}
	}
	if dt <= 0 || math.IsNaN(dt) {
		return fmt.Errorf("invalid time step %g at t=%g", dt, s.Time)
	}

	// Stage 1: U1 = U + dt R(U)
	rhs0, err := s.rightHandSide(dt)
	if err != nil {
		return err
	}
	u1 := make([]Conserved, len(s.Grid.Cells))
	for i, c := range s.Grid.Cells {
		u1[i] = Conserved{
			Rho:  c.U.Rho + dt*rhs0[i].Rho,
			RhoU: c.U.RhoU + dt*rhs0[i].RhoU,
			E:    c.U.E + dt*rhs0[i].E,
		}
		if _, err := g.ConservativeToPrimitive(u1[i]); err != nil {
			return fmt.Errorf("RK stage 1: %w", err)
		}
	}

	// Stage 2: Unew = 0.5 U + 0.5 (U1 + dt R(U1))
	prev := make([]Conserved, len(s.Grid.Cells))
	for i, c := range s.Grid.Cells {
		prev[i] = c.U
		c.U = u1[i]
	}
	rhs1, err := s.rightHandSide(dt)
	if err != nil {
		for i, c := range s.Grid.Cells {
			c.U = prev[i]
		}
		return err
	}
	for i, c := range s.Grid.Cells {
		cand := Conserved{
			Rho:  0.5*prev[i].Rho + 0.5*(u1[i].Rho+dt*rhs1[i].Rho),
			RhoU: 0.5*prev[i].RhoU + 0.5*(u1[i].RhoU+dt*rhs1[i].RhoU),
			E:    0.5*prev[i].E + 0.5*(u1[i].E+dt*rhs1[i].E),
		}
		if _, err := g.ConservativeToPrimitive(cand); err != nil {
			for j, cc := range s.Grid.Cells {
				cc.U = prev[j]
			}
			return fmt.Errorf("RK stage 2: %w", err)
		}
		c.U = cand
	}

	s.Time += dt
	s.Step++
	return nil
}

// Run advances until Time >= TEnd or MaxSteps, persisting snapshots on schedule.
func (s *Simulation) Run() error {
	for {
		if s.Cfg.TEnd > 0 && s.Time >= s.Cfg.TEnd {
			break
		}
		if s.Cfg.MaxSteps > 0 && s.Step >= s.Cfg.MaxSteps {
			break
		}
		if err := s.Advance(); err != nil {
			return err
		}
		if s.Cfg.SnapshotPath != "" && s.Cfg.SnapshotEvery > 0 &&
			s.Step%int64(s.Cfg.SnapshotEvery) == 0 {
			if err := s.SaveSnapshot(s.Cfg.SnapshotPath); err != nil {
				return err
			}
		}
	}
	if s.Cfg.SnapshotPath != "" {
		if err := s.SaveSnapshot(s.Cfg.SnapshotPath); err != nil {
			return err
		}
	}
	return nil
}

// TotalMass/Momentum/Energy integrate conserved variables over the domain.
// Summation follows the left-to-right leaf order (deterministic).
func (s *Simulation) total(comp func(Conserved) float64) float64 {
	var sum float64
	for _, c := range s.Grid.Cells {
		sum += comp(c.U) * c.Width()
	}
	return sum
}

func (s *Simulation) TotalMass() float64 {
	return s.total(func(u Conserved) float64 { return u.Rho })
}
func (s *Simulation) TotalMomentum() float64 {
	return s.total(func(u Conserved) float64 { return u.RhoU })
}
func (s *Simulation) TotalEnergy() float64 {
	return s.total(func(u Conserved) float64 { return u.E })
}
