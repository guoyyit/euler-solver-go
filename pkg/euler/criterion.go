package euler

import "math"

// JumpCriterion flags a cell when the normalised first-order undivided
// difference of any component exceeds RefineThreshold, and clears the flag
// only once the same measure drops below CoarsenThreshold (hysteresis avoids
// rapid refine/coarsen chatter). Both thresholds are user configurable.
//
// The measure is scale-free:
//
//	|w[i+1]-w[i-1]| / (|w[i+1]|+|w[i-1]|+eps)
//
// evaluated on density and pressure (velocity uses a Mach-normalised form).
type JumpCriterion struct {
	RefineThreshold  float64
	CoarsenThreshold float64
	Eps              float64
}

// NewJumpCriterion builds a criterion with the given thresholds.
func NewJumpCriterion(refine, coarsen, eps float64) *JumpCriterion {
	return &JumpCriterion{RefineThreshold: refine, CoarsenThreshold: coarsen, Eps: eps}
}

// Flag implements RefineCriterion.
func (c *JumpCriterion) Flag(lv *Level) []bool {
	n := lv.NCells()
	flag := make([]bool, n)
	lo := NGhost
	for i := 0; i < n; i++ {
		p := lo + i
		qL, qC, qR := lv.Q[p-1], lv.Q[p], lv.Q[p+1]
		meas := math.Abs(qR.Rho-qL.Rho) / (math.Abs(qR.Rho) + math.Abs(qL.Rho) + c.Eps)
		mp := math.Abs(qR.P-qL.P) / (math.Abs(qR.P) + math.Abs(qL.P) + c.Eps)
		if mp > meas {
			meas = mp
		}
		scale := math.Abs(qL.U) + math.Abs(qR.U) + c.Eps
		mu := math.Abs(qR.U-qL.U) / scale
		if mu > meas {
			meas = mu
		}
		_ = qC
		if meas >= c.RefineThreshold {
			flag[i] = true
		} else if meas <= c.CoarsenThreshold {
			flag[i] = false
		}
	}
	return flag
}
