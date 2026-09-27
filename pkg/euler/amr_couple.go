package euler

import (
	"fmt"
	"os"
)

// activeCoarseIntervals returns the coarse cells genuinely covered by the fine
// level (derived from fine.Active), in coarse global index scale.
func (s *Solver) activeCoarseIntervals() []Interval {
	if len(s.mesh.Levels) < 2 {
		return nil
	}
	fine := s.mesh.Levels[1]
	r := s.cfg.AMR.Ratio
	out := make([]Interval, 0, len(fine.Active))
	for _, iv := range fine.Active {
		out = append(out, Interval{Lo: iv.Lo / r, Hi: iv.Hi / r})
	}
	return out
}

// fineActiveMask returns a length-NCells bool mask marking genuinely refined
// fine interior cells (gap cells are false).
func fineActiveMask(fine *Level) []bool {
	mask := make([]bool, fine.NCells())
	for _, iv := range fine.Active {
		for fi := iv.Lo; fi < iv.Hi; fi++ {
			mask[fi-fine.IV.Lo] = true
		}
	}
	return mask
}

// refillFineGaps writes coarse-conservative states into every non-active span
// cell so gaps always equal the supplied coarse snapshot.
func (s *Solver) refillFineGaps(fine *Level, cu []Conserved) {
	for fi := fine.IV.Lo; fi < fine.IV.Hi; fi++ {
		active := false
		for _, iv := range fine.Active {
			if fi >= iv.Lo && fi < iv.Hi {
				active = true
				break
			}
		}
		if !active {
			fine.U[fine.at(fi)] = s.prolongFrom(cu, fi)
		}
	}
}

// restrictFineToCoarse averages r fine children into covered coarse cells for
// the given (fine-scale) intervals.
func (s *Solver) restrictFineToCoarse(fine *Level, fivs []Interval) {
	r := s.cfg.AMR.Ratio
	base := s.mesh.Levels[0]
	for _, iv := range fivs {
		for fi := iv.Lo; fi < iv.Hi; fi += r {
			kids := make([]Conserved, r)
			for k := 0; k < r; k++ {
				kids[k] = fine.U[fine.at(fi+k)]
			}
			base.U[base.at(fi/r)] = restrictChildren(kids)
		}
	}
}

// coupledAdvance advances base+fine by dt with SSP-RK2 and time-consistent
// interface coupling. Only fine.Active cells evolve; span gaps are continuously
// held at coarse-conservative values and serve as stencil data. Flux mismatch
// at each coarse/fine interface is refluxed so the composite update is exactly
// conservative.
func (s *Solver) coupledAdvance(dt float64) error {
	base := s.mesh.Levels[0]
	fine := s.mesh.Levels[1]
	sc := s.stepConfig()
	r := s.cfg.AMR.Ratio
	active := fineActiveMask(fine)

	// Snapshot coarse time-n states.
	cN := make([]Conserved, len(base.U))
	copy(cN, base.U)

	nf := NewFlux(sc.Flux)
	rec := &Reconstructor{EOS: s.eos, Limiter: sc.Limiter}

	// ---- Coarse RK2 (standalone physical ghosts). ----
	if err := fillConservativeGhosts(base, s.eos, sc.Boundary); err != nil {
		return err
	}
	cF0 := faceFluxes(base, s.eos, sc, nf, rec)
	cPred, err := rk2Update(base, s.eos, cN, cF0, dt, nil)
	if err != nil {
		return err
	}
	base.U = cPred
	if err := fillConservativeGhosts(base, s.eos, sc.Boundary); err != nil {
		return err
	}
	cF1 := faceFluxes(base, s.eos, sc, nf, rec)
	cNew, err := rk2Combo(base, s.eos, cN, cF0, cF1, dt, nil)
	if err != nil {
		return err
	}
	base.U = cNew

	// ---- Fine RK2 with stage-matched coarse ghosts; gaps held coarse. ----
	if err := s.fillFineGhostsFrom(fine, cN); err != nil {
		return err
	}
	s.refillFineGaps(fine, cN)
	if err := fine.syncQ(s.eos); err != nil {
		return err
	}
	fN := make([]Conserved, len(fine.U))
	copy(fN, fine.U)
	fF0 := faceFluxes(fine, s.eos, sc, nf, rec)
	fPred, err := rk2Update(fine, s.eos, fN, fF0, dt, active)
	if err != nil {
		return err
	}
	fine.U = fPred
	if err := s.fillFineGhostsFrom(fine, cPred); err != nil {
		return err
	}
	s.refillFineGaps(fine, cPred)
	if err := fine.syncQ(s.eos); err != nil {
		return err
	}
	fF1 := faceFluxes(fine, s.eos, sc, nf, rec)
	fNew, err := rk2Combo(fine, s.eos, fN, fF0, fF1, dt, active)
	if err != nil {
		return err
	}
	fine.U = fNew

	// ---- Reflux at every coarse/fine interface (per Active interval). ----
	civs := s.activeCoarseIntervals()
	d := func(cj, fj int) Conserved {
		return Conserved{
			R: 0.5 * (cF0[cj].R + cF1[cj].R - fF0[fj].R - fF1[fj].R),
			M: 0.5 * (cF0[cj].M + cF1[cj].M - fF0[fj].M - fF1[fj].M),
			E: 0.5 * (cF0[cj].E + cF1[cj].E - fF0[fj].E - fF1[fj].E),
		}
	}
	for _, iv := range civs {
		if iv.Lo-1 >= base.IV.Lo {
			dL := d(base.at(iv.Lo-1), fine.at(iv.Lo*r-1))
			s.addCorrection(base.at(iv.Lo-1), dL, +dt/base.DX)
		}
		if iv.Hi < base.IV.Hi {
			dR := d(base.at(iv.Hi-1), fine.at(iv.Hi*r-1))
			s.addCorrection(base.at(iv.Hi), dR, -dt/base.DX)
		}
	}

	// ---- Covered coarse cells <- conservative fine average. ----
	s.restrictFineToCoarse(fine, fine.Active)
	FillGhostCons(base.U, NGhost, NGhost+base.NCells(), NGhost, sc.Boundary)
	if err := base.syncQ(s.eos); err != nil {
		amrDebugLog(base, civs, err)
		return err
	}
	return s.fillFineGhosts(fine)
}

// addCorrection adds scale*dF to a coarse conserved cell.
func (s *Solver) addCorrection(packed int, dF Conserved, scale float64) {
	u := &s.mesh.Levels[0].U[packed]
	u.R += scale * dF.R
	u.M += scale * dF.M
	u.E += scale * dF.E
}

// uniqueCellsAMR counts every spatial point exactly once: only fine.Active
// cells replace their coarse parents.
func (s *Solver) uniqueCellsAMR() []cellPart {
	base := s.mesh.Levels[0]
	fine := s.mesh.Levels[1]
	r := s.cfg.AMR.Ratio
	covered := make([]bool, base.NCells())
	for _, iv := range fine.Active {
		for ci := iv.Lo / r; ci < iv.Hi/r; ci++ {
			covered[ci-base.IV.Lo] = true
		}
	}
	var coarseIdx []int
	for ci := base.IV.Lo; ci < base.IV.Hi; ci++ {
		if !covered[ci-base.IV.Lo] {
			coarseIdx = append(coarseIdx, base.at(ci))
		}
	}
	var fineIdx []int
	for _, iv := range fine.Active {
		for fi := iv.Lo; fi < iv.Hi; fi++ {
			fineIdx = append(fineIdx, fine.at(fi))
		}
	}
	return []cellPart{{lv: base, idx: coarseIdx}, {lv: fine, idx: fineIdx}}
}

// amrDebugLog prints reflux boundary states when AMR_POSDBG is set.
func amrDebugLog(base *Level, civs []Interval, err error) {
	if os.Getenv("AMR_POSDBG") == "" {
		return
	}
	for _, iv := range civs {
		for _, ci := range []int{iv.Lo - 1, iv.Hi} {
			if ci < 0 || ci >= base.NCells() {
				continue
			}
			fmt.Fprintf(os.Stderr, "  dbg ci=%d U=%+v\n", ci, base.U[base.at(ci)])
		}
	}
	fmt.Fprintf(os.Stderr, "  syncQ: %v\n", err)
}
