package euler

import "math"

func sqrtMath(x float64) float64 { return math.Sqrt(x) }

func mathFloat64bits(f float64) uint64     { return math.Float64bits(f) }
func mathFloat64frombits(u uint64) float64 { return math.Float64frombits(u) }

// minmod returns the minmod of a,b (limited slope).
func minmod(a, b float64) float64 {
	if a*b <= 0 {
		return 0
	}
	if abs(a) < abs(b) {
		return a
	}
	return b
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func infFloat() float64 { return math.Inf(1) }
