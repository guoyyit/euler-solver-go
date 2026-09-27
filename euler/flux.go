package euler

import "math"

// FluxScheme selects the numerical Riemann flux.
type FluxScheme string

const (
	FluxRusanov FluxScheme = "rusanov" // robust local Lax-Friedrichs flux (default)
	FluxHLLC    FluxScheme = "hllc"    // full-wave HLLC contact-resolving flux
)

// PhysicalFlux returns f(U) = [rho u, rho u^2 + p, u(E+p)].
func PhysicalFlux(g Gas, c Conserved, p Primitive) Conserved {
	u := p.U
	return Conserved{
		Rho:  c.Rho * u,
		RhoU: c.Rho*u*u + p.P,
		E:    u * (c.E + p.P),
	}
}

// NumericalFlux evaluates the upwind flux at an interface given reconstructed
// primitive traces from the left and right cells.
func NumericalFlux(g Gas, scheme FluxScheme, left, right Primitive) Conserved {
	switch scheme {
	case FluxHLLC:
		return hllcFlux(g, left, right)
	default:
		return rusanovFlux(g, left, right)
	}
}

// maxWaveSpeed returns max(|u|+c) across two states.
func maxWaveSpeed(g Gas, l, r Primitive) float64 {
	s := math.Abs(l.U) + g.SoundSpeed(l)
	if sr := math.Abs(r.U) + g.SoundSpeed(r); sr > s {
		return sr
	}
	return s
}

// rusanovFlux is F = 0.5(FL+FR) - 0.5 smax (UR-UL), always positivity-friendly.
func rusanovFlux(g Gas, left, right Primitive) Conserved {
	cL, _ := g.PrimitiveToConservative(left)
	cR, _ := g.PrimitiveToConservative(right)
	fL := PhysicalFlux(g, cL, left)
	fR := PhysicalFlux(g, cR, right)
	s := maxWaveSpeed(g, left, right)
	return Conserved{
		Rho:  0.5*(fL.Rho+fR.Rho) - 0.5*s*(cR.Rho-cL.Rho),
		RhoU: 0.5*(fL.RhoU+fR.RhoU) - 0.5*s*(cR.RhoU-cL.RhoU),
		E:    0.5*(fL.E+fR.E) - 0.5*s*(cR.E-cL.E),
	}
}

// hllcFlux implements Toro's HLLC with a primitive-variable (PVRS) wave estimate.
func hllcFlux(g Gas, L, R Primitive) Conserved {
	gm := g.Gamma
	cL, _ := g.PrimitiveToConservative(L)
	cR, _ := g.PrimitiveToConservative(R)
	fL := PhysicalFlux(g, cL, L)
	fR := PhysicalFlux(g, cR, R)

	aL := g.SoundSpeed(L)
	aR := g.SoundSpeed(R)

	// PVRS pressure estimate (Toro eq. 10.60), clamped positive.
	pPV := 0.5*(L.P+R.P) - 0.25*(R.U-L.U)*(L.Rho+R.Rho)*(aL+aR)
	if pPV < 1e-9 {
		pPV = 1e-9
	}

	// Wave speed estimates (Toro eqs. 10.52-10.53).
	var qL, qR float64
	if pPV <= L.P {
		qL = 1
	} else {
		qL = math.Sqrt(1 + 0.5*(gm+1)/gm*(pPV/L.P-1))
	}
	if pPV <= R.P {
		qR = 1
	} else {
		qR = math.Sqrt(1 + 0.5*(gm+1)/gm*(pPV/R.P-1))
	}
	sL := L.U - aL*qL
	sR := R.U + aR*qR

	// Contact wave speed (Toro eq. 10.58).
	den := L.Rho*(sL-L.U) - R.Rho*(sR-R.U)
	if den == 0 {
		return rusanovFlux(g, L, R)
	}
	sS := (R.P - L.P + L.Rho*L.U*(sL-L.U) - R.Rho*R.U*(sR-R.U)) / den

	switch {
	case 0 <= sL:
		return fL
	case sL < 0 && 0 <= sS:
		// Left star state.
		fac := L.Rho * (sL - L.U) / (sL - sS)
		uStar := Conserved{
			Rho:  fac,
			RhoU: fac * sS,
			E:    fac * (cL.E/L.Rho + (sS-L.U)*(sS+L.P/(L.Rho*(sL-L.U)))),
		}
		return Conserved{
			Rho:  fL.Rho + sL*(uStar.Rho-cL.Rho),
			RhoU: fL.RhoU + sL*(uStar.RhoU-cL.RhoU),
			E:    fL.E + sL*(uStar.E-cL.E),
		}
	case sS < 0 && 0 <= sR:
		// Right star state.
		fac := R.Rho * (sR - R.U) / (sR - sS)
		uStar := Conserved{
			Rho:  fac,
			RhoU: fac * sS,
			E:    fac * (cR.E/R.Rho + (sS-R.U)*(sS+R.P/(R.Rho*(sR-R.U)))),
		}
		return Conserved{
			Rho:  fR.Rho + sR*(uStar.Rho-cR.Rho),
			RhoU: fR.RhoU + sR*(uStar.RhoU-cR.RhoU),
			E:    fR.E + sR*(uStar.E-cR.E),
		}
	default:
		return fR
	}
}
