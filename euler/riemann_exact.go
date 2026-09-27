package euler

import "math"

// RiemannData describes an ideal-gas Riemann problem with left/right
// primitive states, heat ratio and initial discontinuity location x0.
type RiemannData struct {
	PL, PR Primitive
	Gamma  float64
	X0     float64
}

// ExactState evaluates the exact solution at position x and time t>0
// (Toro, "Riemann Solvers and Numerical Methods", Chap. 4).
func (r *RiemannData) ExactState(x, t float64) Primitive {
	if t <= 0 {
		if x < r.X0 {
			return r.PL
		}
		return r.PR
	}
	g := r.Gamma
	g1, g2, g3, g4 := 0.5*(g-1)/g, 0.5*(g+1)/g, 2*g/(g-1), 2/(g-1)
	g5 := 2 / (g + 1)
	g6 := (g - 1) / (g + 1)
	g7 := 0.5 * (g - 1)
	g9 := 1 / g

	pStar, uStar := r.starState(g1, g2, g3, g4, g5, g6)

	cL := math.Sqrt(g * r.PL.P / r.PL.Rho)
	cR := math.Sqrt(g * r.PR.P / r.PR.Rho)
	xi := (x - r.X0) / t

	// Left wave.
	if xi < uStar {
		if pStar <= r.PL.P { // rarefaction
			sHL := r.PL.U - cL
			if xi < sHL {
				return r.PL
			}
			cML := cL * math.Pow(pStar/r.PL.P, g1)
			sTL := uStar - cML
			if xi < sTL { // fan
				fac := g5 + g6/cL*(r.PL.U-xi)
				ro := r.PL.Rho * math.Pow(fac, g4)
				u := g5 * (cL + g7*r.PL.U + xi)
				p := r.PL.P * math.Pow(fac, g3)
				return Primitive{ro, u, p}
			}
		} else { // shock
			pML := pStar / r.PL.P
			sL := r.PL.U - cL*math.Sqrt(g2*pML+g1)
			if xi < sL {
				return r.PL
			}
		}
		// left star
		var ro float64
		if pStar <= r.PL.P {
			ro = r.PL.Rho * math.Pow(pStar/r.PL.P, g9)
		} else {
			ro = r.PL.Rho * (pMLShock(pStar, r.PL.P) + g6) / (pMLShock(pStar, r.PL.P)*g6 + 1)
		}
		return Primitive{ro, uStar, pStar}
	}

	// Right wave.
	if pStar >= r.PR.P { // shock
		pMR := pStar / r.PR.P
		sR := r.PR.U + cR*math.Sqrt(g2*pMR+g1)
		if xi > sR {
			return r.PR
		}
		ro := r.PR.Rho * (pMR + g6) / (pMR*g6 + 1)
		return Primitive{ro, uStar, pStar}
	}
	// rarefaction
	sHR := r.PR.U + cR
	if xi > sHR {
		return r.PR
	}
	cMR := cR * math.Pow(pStar/r.PR.P, g1)
	sTR := uStar + cMR
	if xi > sTR { // fan
		fac := g5 - g6/cR*(r.PR.U-xi)
		ro := r.PR.Rho * math.Pow(fac, g4)
		u := g5 * (-cR + g7*r.PR.U + xi)
		p := r.PR.P * math.Pow(fac, g3)
		return Primitive{ro, u, p}
	}
	ro := r.PR.Rho * math.Pow(pStar/r.PR.P, g9)
	return Primitive{ro, uStar, pStar}
}

func pMLShock(pStar, p float64) float64 { return pStar / p }

// starState solves the pressure function f_L(p)+f_R(p)+(uR-uL)=0 for (p*,u*).
func (r *RiemannData) starState(g1, g2, g3, g4, g5, g6 float64) (p, u float64) {
	cL := math.Sqrt(r.Gamma * r.PL.P / r.PL.Rho)
	cR := math.Sqrt(r.Gamma * r.PR.P / r.PR.Rho)

	// PVRS initial guess.
	pPV := 0.5*(r.PL.P+r.PR.P) - 0.125*(r.PR.U-r.PL.U)*(r.PL.Rho+r.PR.Rho)*(cL+cR)
	if pPV < 1e-9 {
		pPV = 1e-9
	}
	p = pPV
	var fL, fR, dfL, dfR float64
	for iter := 0; iter < 30; iter++ {
		fL, dfL = preFun(p, r.PL.Rho, cL, r.PL.P, r.PL.U, g1, g2, g4, g5)
		fR, dfR = preFun(p, r.PR.Rho, cR, r.PR.P, r.PR.U, g1, g2, g4, g5)
		pn := p - (fL+fR+r.PR.U-r.PL.U)/(dfL+dfR)
		if pn < 0 {
			pn = 1e-9
		}
		if 2*math.Abs(pn-p)/(pn+p) < 1e-12 {
			p = pn
			break
		}
		p = pn
	}
	u = 0.5*(r.PL.U+r.PR.U) + 0.5*(fR-fL)
	return p, u
}

// preFun evaluates f_K(p) and its derivative for K=L or R.
func preFun(p, dK, cK, pK, uK, g1, g2, g4, g5 float64) (f, fp float64) {
	if p <= pK { // rarefaction
		pr := p / pK
		f = g4 * cK * (math.Pow(pr, g1) - 1)
		fp = math.Pow(pr, -g2) / (dK * cK)
		return f, fp
	}
	// shock
	aK := g5 / dK
	bK := g1 * pK
	sqrt := math.Sqrt(aK / (p + bK))
	f = (p - pK) * sqrt
	fp = sqrt * (1 - 0.5*(p-pK)/(p+bK))
	return f, fp
}
