package euler

import (
	"math"
	"testing"
)

func approxCons(a, b Conserved, tol float64) bool {
	d := func(x, y float64) float64 { return math.Abs(x - y) }
	return d(a.R, b.R) <= tol && d(a.M, b.M) <= tol && d(a.E, b.E) <= tol
}

func TestFluxConsistency(t *testing.T) {
	eos := NewEOS(1.4)
	q := Primitive{Rho: 1.2, U: 0.7, P: 2.1}
	phys := eos.PrimFlux(q)
	for _, sc := range []FluxScheme{FluxRusanov, FluxHLL} {
		f := NewFlux(sc).Flux(eos, q, q)
		if !approxCons(f, phys, 1e-14) {
			t.Fatalf("%v: numerical flux %v != physical flux %v", sc, f, phys)
		}
	}
}

func TestFluxMonotonicityAtSupersonic(t *testing.T) {
	eos := NewEOS(1.4)
	qL := Primitive{Rho: 1, U: 10, P: 1}
	qR := Primitive{Rho: 0.5, U: 9, P: 0.8}
	f := hllFlux{}.Flux(eos, qL, qR)
	if f != eos.PrimFlux(qL) {
		t.Fatalf("supersonic HLL should upwind: %v", f)
	}
}
