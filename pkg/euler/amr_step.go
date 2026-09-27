package euler

// floorDiv returns Euclidean floor(a/b) for b>0 (Go's / truncates to zero).
func floorDiv(a, b int) int {
	q := a / b
	if a%b < 0 {
		q--
	}
	return q
}

// dilate grows the true entries of mask by b cells on each side. Dilation is
// non-wrapping even on a periodic domain: the later interval construction
// peels seam-touching runs, so a circular wrap would incorrectly flag the
// entire domain whenever a flag landed near a seam.
func dilate(mask []bool, b int, periodic bool) []bool {
	n := len(mask)
	out := make([]bool, n)
	for i, v := range mask {
		if !v {
			continue
		}
		lo := i - b
		if lo < 0 {
			lo = 0
		}
		hi := i + b
		if hi >= n {
			hi = n - 1
		}
		for j := lo; j <= hi; j++ {
			out[j] = true
		}
	}
	return out
}

// maskIntervals converts a boolean mask into merged intervals.
func maskIntervals(mask []bool) []Interval {
	var ivs []Interval
	i := 0
	for i < len(mask) {
		if !mask[i] {
			i++
			continue
		}
		lo := i
		for i < len(mask) && mask[i] {
			i++
		}
		ivs = append(ivs, Interval{Lo: lo, Hi: i})
	}
	return ivs
}

// refineCoarseIntervals returns the coarse intervals selected for refinement.
// On a periodic domain intervals are kept strictly interior (at least
// buffer+1 cells away from either seam) unless essentially the whole domain is
// flagged, in which case the entire span is refined. Keeping every coarse/fine
// interface interior lets a single conservative reflux correction close the
// global balance exactly, including flux crossing the periodic seam.
// Non-seam components closer than 2*buffer+1 are bridged so the fine span has
// no internal coarse/fine gaps (which would be non-conservative at regrid).
func refineCoarseIntervals(mask []bool, periodic bool, buffer int) []Interval {
	n := len(mask)
	ivs := maskIntervals(mask)
	if periodic {
		// Keep the fine level strictly interior: peel any run touching a
		// physical seam. Then every coarse/fine interface is an ordinary
		// interior face and a single reflux correction closes the global
		// balance exactly; the seam stays coarse, where periodic telescoping
		// fluxes cancel by construction.
		var in []Interval
		for _, iv := range ivs {
			if iv.Lo == 0 {
				iv.Lo = buffer + 1
			}
			if iv.Hi == n {
				iv.Hi = n - buffer - 1
			}
			if iv.Hi > iv.Lo {
				in = append(in, iv)
			}
		}
		ivs = in
	}
	bridge := 2*buffer + 1
	if len(ivs) <= 1 {
		return ivs
	}
	merged := []Interval{ivs[0]}
	for _, iv := range ivs[1:] {
		last := &merged[len(merged)-1]
		if iv.Lo-last.Hi <= bridge {
			last.Hi = iv.Hi
		} else {
			merged = append(merged, iv)
		}
	}
	return merged
}

// adapt rebuilds the fine level(s) from current data. Existing fine cells that
// survive are copied bit-for-bit; newly refined cells are injected
// conservatively from the coarse parent; coarsened cells disappear and their
// content already lives in the coarse level via restriction. One Level object
// per AMR depth is allocated, spanning the union range; gap cells are
// coarse-injected but masked out of integrals.
func (s *Solver) adapt() error {
	cfg := s.cfg.AMR
	base := s.mesh.Levels[0]
	n := base.NCells()
	if cfg.MaxLevels <= 1 {
		s.mesh.Levels = s.mesh.Levels[:1]
		return nil
	}
	if len(s.mesh.Levels) > 1 {
		s.restrictFineToCoarse(s.mesh.Levels[1], s.fineIntervals())
	}

	// The regrid decision must depend only on serialisable parameters so a
	// restored run follows the identical mesh trajectory.
	crit := NewJumpCriterion(cfg.RefineThreshold, cfg.CoarsenThreshold, cfg.CriterionEps)
	mask := crit.Flag(base)
	mask = dilate(mask, cfg.MergeBufferCells, s.cfg.Boundary == BoundaryPeriodic)
	if s.cfg.Boundary == BoundaryReflecting && n > 2 {
		mask[0] = false
		mask[n-1] = false
	}
	civs := refineCoarseIntervals(mask, s.cfg.Boundary == BoundaryPeriodic, cfg.MergeBufferCells)
	if len(civs) == 0 {
		s.mesh.Levels = s.mesh.Levels[:1]
		return nil
	}
	r := cfg.Ratio
	lo, hi := civs[0].Lo, civs[0].Hi
	for _, iv := range civs[1:] {
		if iv.Lo < lo {
			lo = iv.Lo
		}
		if iv.Hi > hi {
			hi = iv.Hi
		}
	}
	flo, fhi := lo*r, hi*r
	fdx := s.mesh.BaseDX() / float64(r)
	fine := NewLevel(flo, fhi, fdx)
	for _, iv := range civs {
		fine.Active = append(fine.Active, Interval{Lo: iv.Lo * r, Hi: iv.Hi * r})
	}
	var old *Level
	if len(s.mesh.Levels) > 1 {
		old = s.mesh.Levels[1]
	}
	oldActive := func(fi int) bool {
		if old == nil {
			return false
		}
		for _, iv := range old.Active {
			if fi >= iv.Lo && fi < iv.Hi {
				return true
			}
		}
		return false
	}
	for fi := flo; fi < fhi; fi++ {
		p := fine.at(fi)
		if oldActive(fi) {
			fine.U[p] = old.U[old.at(fi)] // preserved bit-for-bit
			continue
		}
		ci := floorDiv(fi, r)
		fine.U[p] = s.coarseConsFor(ci) // conservative injection (incl. gaps)
	}
	if err := s.fillFineGhosts(fine); err != nil {
		return err
	}
	s.mesh.Levels = []*Level{base, fine}
	return nil
}

// coarseConsFor returns the base-level conserved state of coarse cell ci.
func (s *Solver) coarseConsFor(ci int) Conserved {
	base := s.mesh.Levels[0]
	return base.U[base.at(ci)]
}

// fineIntervals returns the genuinely refined intervals (fine global scale).
func (s *Solver) fineIntervals() []Interval {
	if len(s.mesh.Levels) < 2 {
		return nil
	}
	return s.mesh.Levels[1].Active
}

// fillFineGhosts populates fine ghosts from the current coarse state.
func (s *Solver) fillFineGhosts(fine *Level) error {
	return s.fillFineGhostsFrom(fine, s.mesh.Levels[0].U)
}

// prolongFrom returns the r-way child of coarse parent ci (in coarse array cu)
// containing fine global index fi.
func (s *Solver) prolongFrom(cu []Conserved, fi int) Conserved {
	r := s.cfg.AMR.Ratio
	base := s.mesh.Levels[0]
	ci := floorDiv(fi, r)
	parent := cu[base.at(ci)]
	return prolongCell(parent, r)[fi-ci*r]
}

// fillFineGhostsFrom populates fine ghosts using the supplied coarse conserved
// array (which may be a time-n or predictor snapshot). Interior ghosts are
// prolongations of coarse parents; physical-boundary ghosts apply the
// configured treatment (mirror at a wall, wrap on a periodic seam).
func (s *Solver) fillFineGhostsFrom(fine *Level, cu []Conserved) error {
	r := s.cfg.AMR.Ratio
	nf := s.mesh.NBase * r
	bc := s.cfg.Boundary

	// Gaps inside the allocated span always track the supplied coarse state.
	s.refillFineGaps(fine, cu)

	for pos := fine.IV.Lo - NGhost; pos < fine.IV.Lo; pos++ {
		fine.U[fine.at(pos)] = s.prolongFrom(cu, pos)
	}
	for pos := fine.IV.Hi; pos < fine.IV.Hi+NGhost; pos++ {
		fine.U[fine.at(pos)] = s.prolongFrom(cu, pos)
	}

	seam := func(fi int) Conserved {
		if fi >= fine.IV.Lo && fi < fine.IV.Hi {
			return fine.U[fine.at(fi)]
		}
		return s.prolongFrom(cu, fi)
	}
	if fine.IV.Lo == 0 {
		for k := 1; k <= NGhost; k++ {
			if bc == BoundaryPeriodic {
				fine.U[fine.at(-k)] = seam(nf - k)
			} else {
				v := fine.U[fine.at(k-1)]
				fine.U[fine.at(-k)] = Conserved{R: v.R, M: -v.M, E: v.E}
			}
		}
	}
	if fine.IV.Hi == nf {
		for k := 0; k < NGhost; k++ {
			if bc == BoundaryPeriodic {
				fine.U[fine.at(nf+k)] = seam(k)
			} else {
				v := fine.U[fine.at(nf-1-k)]
				fine.U[fine.at(nf+k)] = Conserved{R: v.R, M: -v.M, E: v.E}
			}
		}
	}
	return fine.syncQ(s.eos)
}
