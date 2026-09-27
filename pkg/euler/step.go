package euler

import "fmt"

// StepConfig groups the spatial-discretisation options shared by all levels.
type StepConfig struct {
	Flux     FluxScheme
	Limiter  Limiter
	Boundary BoundaryKind
}

// faceFluxes computes the numerical flux at every face bounding the interior
// cells. Ghost cells and lv.Q must already be populated for the stage.
func faceFluxes(lv *Level, eos *EOS, sc StepConfig, nf NumericalFlux, rec *Reconstructor) []Conserved {
	n := lv.NCells()
	lo := NGhost
	qL, qR := rec.Reconstruct(lv.Q, lo, lo+n)
	F := make([]Conserved, len(lv.U))
	for j := lo - 1; j <= lo+n-1; j++ {
		F[j] = nf.Flux(eos, qL[j], qR[j])
	}
	return F
}

// fillConservativeGhosts applies the configured physical-boundary ghost fill
// and rebuilds primitive states. It is the standalone-grid ghost callback.
func fillConservativeGhosts(lv *Level, eos *EOS, bc BoundaryKind) error {
	FillGhostCons(lv.U, NGhost, NGhost+lv.NCells(), NGhost, bc)
	return lv.syncQ(eos)
}

// rk2Update computes U -= dt/dx*(F_right-F_left) on masked interior cells.
// Cells for which active[i]==false are copied unchanged (they are coarse-filled
// gaps used only as fine stencil/ghost values).
func rk2Update(lv *Level, eos *EOS, base []Conserved, F []Conserved, dt float64, active []bool) ([]Conserved, error) {
	n := lv.NCells()
	lo := NGhost
	out := make([]Conserved, len(base))
	copy(out, base)
	for i := 0; i < n; i++ {
		p := lo + i
		if active != nil && !active[i] {
			continue
		}
		out[p] = Conserved{
			R: base[p].R - dt*(F[p].R-F[p-1].R)/lv.DX,
			M: base[p].M - dt*(F[p].M-F[p-1].M)/lv.DX,
			E: base[p].E - dt*(F[p].E-F[p-1].E)/lv.DX,
		}
		if _, err := eos.ConsToPrim(out[p]); err != nil {
			return nil, fmt.Errorf("rk stage at cell %d: %w", i, err)
		}
	}
	return out, nil
}

// rk2Combo forms the Heun combination on masked interior cells only.
func rk2Combo(lv *Level, eos *EOS, u0, F0, F1 []Conserved, dt float64, active []bool) ([]Conserved, error) {
	n := lv.NCells()
	lo := NGhost
	out := make([]Conserved, len(u0))
	copy(out, u0)
	for i := 0; i < n; i++ {
		p := lo + i
		if active != nil && !active[i] {
			continue
		}
		c := func(a float64) float64 {
			return a - 0.5*dt*((F0[p].R-F0[p-1].R)+(F1[p].R-F1[p-1].R))/lv.DX
		}
		_ = c
		comb := func(get func(Conserved) float64) float64 {
			d0 := get(F0[p]) - get(F0[p-1])
			d1 := get(F1[p]) - get(F1[p-1])
			return get(u0[p]) - 0.5*dt*(d0+d1)/lv.DX
		}
		out[p] = Conserved{R: comb(func(c Conserved) float64 { return c.R }),
			M: comb(func(c Conserved) float64 { return c.M }),
			E: comb(func(c Conserved) float64 { return c.E })}
		if _, err := eos.ConsToPrim(out[p]); err != nil {
			return nil, fmt.Errorf("rk final at cell %d: %w", i, err)
		}
	}
	return out, nil
}

// rk2StepWithGhost advances lv by one SSP-RK2 (Heun) step. fill0 supplies the
// ghost states (and rebuilds Q) for stage 0 evaluated at time n; fill1 does so
// for stage 1 evaluated at the time-n+1 predictor. On any non-physical stage
// the level is restored to its input state.
//
//	U1 = U0 + dt R(U0)
//	U2 = U0 + 0.5 dt (R(U0) + R(U1))
func rk2StepWithGhost(lv *Level, eos *EOS, sc StepConfig, dt float64,
	fill0, fill1 func(*Level) error, active []bool) (F0, F1 []Conserved, err error) {
	nf := NewFlux(sc.Flux)
	rec := &Reconstructor{EOS: eos, Limiter: sc.Limiter}
	saved := make([]Conserved, len(lv.U))
	copy(saved, lv.U)
	savedQ := make([]Primitive, len(lv.Q))
	copy(savedQ, lv.Q)
	restore := func() {
		lv.U = saved
		lv.Q = savedQ
	}

	// Stage 0 at time n.
	if err = fill0(lv); err != nil {
		restore()
		return nil, nil, err
	}
	F0 = faceFluxes(lv, eos, sc, nf, rec)
	u1, err := rk2Update(lv, eos, saved, F0, dt, active)
	if err != nil {
		restore()
		return nil, nil, err
	}
	lv.U = u1
	// Stage 1 at predictor time n+1.
	if err = fill1(lv); err != nil {
		restore()
		return nil, nil, err
	}
	F1 = faceFluxes(lv, eos, sc, nf, rec)
	u2, err := rk2Combo(lv, eos, saved, F0, F1, dt, active)
	if err != nil {
		restore()
		return nil, nil, err
	}
	lv.U = u2
	if err = fill1(lv); err != nil {
		restore()
		return nil, nil, err
	}
	return F0, F1, nil
}

// rk2Step is the standalone-grid convenience wrapper using physical ghosts.
func rk2Step(lv *Level, eos *EOS, sc StepConfig, dt float64) (F0, F1 []Conserved, err error) {
	fill := func(l *Level) error {
		return fillConservativeGhosts(l, eos, sc.Boundary)
	}
	return rk2StepWithGhost(lv, eos, sc, dt, fill, fill, nil)
}

// maxWaveSpeed returns max_i(|u_i|+c_i) over the interior cells of lv.
func maxWaveSpeed(lv *Level, eos *EOS) float64 {
	lo, hi := NGhost, NGhost+lv.NCells()
	s := 0.0
	for i := lo; i < hi; i++ {
		v := abs(lv.Q[i].U) + eos.SoundSpeed(lv.Q[i])
		if v > s {
			s = v
		}
	}
	return s
}

// StableDT returns the largest explicit time step satisfying the CFL condition
// dt <= cfl * dx / smax for every level of the mesh. Zero means no estimate.
func StableDT(m *Mesh, eos *EOS, cfl float64) float64 {
	dt := infFloat()
	for _, lv := range m.Levels {
		if lv.NCells() == 0 {
			continue
		}
		s := maxWaveSpeed(lv, eos)
		if s <= 0 {
			continue
		}
		d := cfl * lv.DX / s
		if d < dt {
			dt = d
		}
	}
	if dt == infFloat() {
		return 0
	}
	return dt
}

// advanceLevel performs one SSP-RK2 step on a single level.
func advanceLevel(lv *Level, eos *EOS, sc StepConfig, dt float64) error {
	_, _, err := rk2Step(lv, eos, sc, dt)
	return err
}
