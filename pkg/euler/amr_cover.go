package euler

import "sort"

// Interval on a level is a contiguous run of global coarse (level-0) indices
// [Lo,Hi). Fine intervals are stored at their own level's global index scale.
// (The Interval type in grid.go is reused.)

// unionIntervals merges possibly overlapping/touching sorted intervals.
func unionIntervals(ivs []Interval) []Interval {
	if len(ivs) == 0 {
		return nil
	}
	cp := make([]Interval, len(ivs))
	copy(cp, ivs)
	sort.Slice(cp, func(i, j int) bool {
		if cp[i].Lo != cp[j].Lo {
			return cp[i].Lo < cp[j].Lo
		}
		return cp[i].Hi < cp[j].Hi
	})
	out := []Interval{cp[0]}
	for _, iv := range cp[1:] {
		last := &out[len(out)-1]
		if iv.Lo <= last.Hi {
			if iv.Hi > last.Hi {
				last.Hi = iv.Hi
			}
			continue
		}
		out = append(out, iv)
	}
	return out
}

// expandIntervals grows every interval by b cells on each side, clips to
// [0,n), and re-merges.
func expandIntervals(ivs []Interval, b, n int) []Interval {
	var out []Interval
	for _, iv := range ivs {
		lo := iv.Lo - b
		if lo < 0 {
			lo = 0
		}
		hi := iv.Hi + b
		if hi > n {
			hi = n
		}
		out = append(out, Interval{Lo: lo, Hi: hi})
	}
	return unionIntervals(out)
}

// partitionCells returns level l's packed indices that are NOT covered by the
// next finer level, i.e. the cells that count uniquely for this level.
// fineIV is given on level l+1's global index scale; ratio aligned fine runs
// [r*i, r*(i+1)) cover coarse cell i.
func partitionCells(lv *Level, l int, m *Mesh, fineIV []Interval) []int {
	r := m.Ratio
	covered := func(coarseI int) bool {
		fl := coarseI * r
		fr := fl + r
		for _, iv := range fineIV {
			if iv.Lo < fr && fl < iv.Hi {
				return true
			}
		}
		return false
	}
	var idx []int
	for i := lv.IV.Lo; i < lv.IV.Hi; i++ {
		if !covered(i) {
			idx = append(idx, lv.at(i))
		}
	}
	return idx
}
