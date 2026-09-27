package euler

// Limiter selects the slope limiter used in MUSCL reconstruction.
// Both choices are total-variation diminishing, so no spurious extrema
// (and hence no non-physical overshoots) are generated across shocks or
// rarefactions.
type Limiter string

const (
	LimiterMinmod  Limiter = "minmod" // most robust; default
	LimiterVanLeer Limiter = "vanleer"
)

// minmod returns the monotone slope: 0 if arguments differ in sign,
// otherwise the smaller-magnitude argument.
func minmod(a, b float64) float64 {
	if a*b <= 0 {
		return 0
	}
	if absF(a) < absF(b) {
		return a
	}
	return b
}

func absF(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

// vanLeerSlope is the symmetric van Leer limiter slope.
func vanLeerSlope(a, b float64) float64 {
	if a*b <= 0 {
		return 0
	}
	return 2 * a * b / (a + b)
}

// Slope computes the limited slope given the one-sided differences
// a = c-l and b = r-c (equal cell widths are assumed by the caller).
func Slope(lim Limiter, l, c, r float64) float64 {
	switch lim {
	case LimiterVanLeer:
		return vanLeerSlope(c-l, r-c)
	default:
		return minmod(c-l, r-c)
	}
}

// MUSCLFace reconstructs primitive traces at the shared face of the center
// and right cells from three consecutive equal-width cell averages:
//
//	cR = c + 0.5*slope(l,c,r)   (right-going trace from center)
//	rL = r - 0.5*slope(c,r,rR)  needs rR -> caller passes rR as fourth...
//
// To keep a compact signature, this helper reconstructs the center cell's
// right trace and the right cell's left trace using slopes evaluated on the
// stencils (l,c,r) and (c,r,rR). It therefore takes FOUR values.
func MUSCLFace(lim Limiter, l, c, r, rR Primitive) (cR, rL Primitive) {
	cR = Primitive{
		Rho: c.Rho + 0.5*Slope(lim, l.Rho, c.Rho, r.Rho),
		U:   c.U + 0.5*Slope(lim, l.U, c.U, r.U),
		P:   c.P + 0.5*Slope(lim, l.P, c.P, r.P),
	}
	rL = Primitive{
		Rho: r.Rho - 0.5*Slope(lim, c.Rho, r.Rho, rR.Rho),
		U:   r.U - 0.5*Slope(lim, c.U, r.U, rR.U),
		P:   r.P - 0.5*Slope(lim, c.P, r.P, rR.P),
	}
	return cR, rL
}
