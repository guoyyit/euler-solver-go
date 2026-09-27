package euler

// FluxScheme selects a numerical flux function.
type FluxScheme int

const (
	// FluxRusanov is the local Lax-Friedrichs (Rusanov) flux.
	FluxRusanov FluxScheme = iota
	// FluxHLL is the Harten-Lax-van Leer flux.
	FluxHLL
)

// NumericalFlux evaluates the convective flux across one face.
type NumericalFlux interface {
	Flux(eos *EOS, qL, qR Primitive) Conserved
}

// rusanovFlux is F = 0.5*(FL+FR) - 0.5*smax*(UR-UL) with
// smax = max(|uL|+cL, |uR|+cR). It is monotone and introduces no
// non-physical oscillations, at the cost of extra numerical diffusion.
type rusanovFlux struct{}

func (rusanovFlux) Flux(eos *EOS, qL, qR Primitive) Conserved {
	uL := eos.PrimToCons(qL)
	uR := eos.PrimToCons(qR)
	fL := eos.PrimFlux(qL)
	fR := eos.PrimFlux(qR)
	sL := abs(qL.U) + eos.SoundSpeed(qL)
	sR := abs(qR.U) + eos.SoundSpeed(qR)
	s := sL
	if sR > s {
		s = sR
	}
	return Conserved{
		R: 0.5*(fL.R+fR.R) - 0.5*s*(uR.R-uL.R),
		M: 0.5*(fL.M+fR.M) - 0.5*s*(uR.M-uL.M),
		E: 0.5*(fL.E+fR.E) - 0.5*s*(uR.E-uL.E),
	}
}

// hllFlux uses the Davis wave-speed estimates:
//
//	sL=min(uL-cL,uR-cR), sR=max(uL+cL,uR+cR).
type hllFlux struct{}

func (hllFlux) Flux(eos *EOS, qL, qR Primitive) Conserved {
	fL := eos.PrimFlux(qL)
	fR := eos.PrimFlux(qR)
	uL := eos.PrimToCons(qL)
	uR := eos.PrimToCons(qR)
	slL := qL.U - eos.SoundSpeed(qL)
	slR := qR.U - eos.SoundSpeed(qR)
	srL := qL.U + eos.SoundSpeed(qL)
	srR := qR.U + eos.SoundSpeed(qR)
	sMin := slL
	if slR < sMin {
		sMin = slR
	}
	sMax := srL
	if srR > sMax {
		sMax = srR
	}
	switch {
	case sMin >= 0:
		return fL
	case sMax <= 0:
		return fR
	default:
		// (sR*FL - sL*FR + sL*sR*(UR-UL)) / (sR - sL)
		d := sMax - sMin
		return Conserved{
			R: (sMax*fL.R - sMin*fR.R + sMin*sMax*(uR.R-uL.R)) / d,
			M: (sMax*fL.M - sMin*fR.M + sMin*sMax*(uR.M-uL.M)) / d,
			E: (sMax*fL.E - sMin*fR.E + sMin*sMax*(uR.E-uL.E)) / d,
		}
	}
}

// NewFlux builds the numerical flux selected by scheme.
func NewFlux(scheme FluxScheme) NumericalFlux {
	switch scheme {
	case FluxHLL:
		return hllFlux{}
	default:
		return rusanovFlux{}
	}
}
