package euler

const testICName = "__test_sod__"

func init() {
	RegisterProblem(Problem{Name: testICName, IC: func(x float64) Primitive {
		if x < 0.5 {
			return Primitive{Rho: 1.0, U: 0, P: 1.0}
		}
		return Primitive{Rho: 0.125, U: 0, P: 0.1}
	}})
}

func testConfig(n int, boundary BoundaryType) Config {
	return Config{
		ICName:   testICName,
		Xmin:     0,
		Xmax:     1,
		N0:       n,
		Gamma:    1.4,
		Boundary: boundary,
		Flux:     FluxRusanov,
		Limiter:  LimiterMinmod,
		CFL:      0.8,
		TEnd:     0.2,
	}
}

func amrConfig(n int, boundary BoundaryType) Config {
	c := testConfig(n, boundary)
	c.Refine = RefineConfig{
		Enabled: true, MinLevel: 0, MaxLevel: 3,
		RefineAt: 0.05, CoarsenAt: 0.01,
		Criterion: RefineDensityGradient,
	}
	c.RegridEvery = 4
	return c
}
