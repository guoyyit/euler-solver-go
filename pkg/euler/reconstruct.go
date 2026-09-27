package euler

// Limiter identifies a slope-limiter for piecewise-linear reconstruction.
type Limiter int

const (
	// LimitNone gives an unlimited centred slope (2nd order, not TVD).
	LimitNone Limiter = iota
	// LimitMinmod is the robust minmod limiter (diffusive, always monotone).
	LimitMinmod
	// LimitMC is the monotonized-centred limiter (sharper, still TVD).
	LimitMC
)

// Reconstructor turns cell averages into left/right face states using MUSCL
// reconstruction applied component by component in primitive variables.
//
// Both minmod and MC limiters are total variation diminishing: every
// extrapolated face value lies inside the three-cell stencil, so positive
// stencil states produce positive face states.
type Reconstructor struct {
	EOS     *EOS
	Limiter Limiter
}

// limitSlope computes the limited slope given left (b=a-prev) and right
// (c=next-a) differences.
func (r *Reconstructor) limitSlope(b, c float64) float64 {
	switch r.Limiter {
	case LimitNone:
		return 0.5 * (b + c)
	case LimitMC:
		s := 0.5 * (b + c)
		return minmod(s, minmod(2.0*b, 2.0*c))
	default:
		return minmod(b, c)
	}
}

// Reconstruct fills face states for every face j in [lo-1, hi). qL[j] is the
// state on the left side of face j (extrapolated from cell j) and qR[j] the
// state on its right (extrapolated from cell j+1). The arrays are indexed like
// w and only faces [lo-1,hi) are written; ghost cells must be populated.
func (r *Reconstructor) Reconstruct(w []Primitive, lo, hi int) (qL, qR []Primitive) {
	n := len(w)
	qL = make([]Primitive, n)
	qR = make([]Primitive, n)
	for j := lo - 1; j < hi; j++ {
		iL := j     // cell left of face j (face between j and j+1)
		iR := j + 1 // cell right of face
		sL := r.slope(w, iL)
		sR := r.slope(w, iR)
		qL[j] = Primitive{
			Rho: w[iL].Rho + 0.5*sL.rho,
			U:   w[iL].U + 0.5*sL.u,
			P:   w[iL].P + 0.5*sL.p,
		}
		qR[j] = Primitive{
			Rho: w[iR].Rho - 0.5*sR.rho,
			U:   w[iR].U - 0.5*sR.u,
			P:   w[iR].P - 0.5*sR.p,
		}
	}
	return qL, qR
}

type primSlope struct{ rho, u, p float64 }

func (r *Reconstructor) slope(w []Primitive, i int) primSlope {
	return primSlope{
		rho: r.limitSlope(w[i].Rho-w[i-1].Rho, w[i+1].Rho-w[i].Rho),
		u:   r.limitSlope(w[i].U-w[i-1].U, w[i+1].U-w[i].U),
		p:   r.limitSlope(w[i].P-w[i-1].P, w[i+1].P-w[i].P),
	}
}
