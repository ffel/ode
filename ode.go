// Package ode provides numerical methods for ordinary differential equations
package ode

// Num as short for float64
type Num float64

// Result is a row in a table of calculated results
type Result struct {
	t  Num
	xx []Num
}

// Ode is a first order differential equation
type Ode func([]Num, Num) Num

// Integrator is the integration method, like mid point
type Integrator func([]Num, Num, Num, []Ode) []Num

// FixedStep iterates over a set of ode's with fixed step h
func FixedStep(method Integrator, dxdt []Ode, xx []Num, t0, tmax, h Num) []Result {
	T := t0

	r := make([]Result, 0, 200)

	for T <= tmax {
		kk := method(xx, T, h, dxdt)

		// store a copy of xx in r, not x itself for that will change
		x := make([]Num, len(xx))
		copy(x, xx)

		r = append(r, Result{T, x})

		for i, k := range kk {
			xx[i] += k
		}

		T += h
	}

	return r
}

// AdaptiveStep iterates over a set of ode's with adaptive h
// starts with h and minimum hmin
func AdaptiveStep(method Integrator, dxdt []Ode, xx []Num, t0, tmax, hmin, h Num) []Result {
	T := t0

	var kkFull []Num

	var H Num

	r := make([]Result, 0, 200)

	for T <= tmax {

		// max 5 decrements
		for a := 0; a < 5; a++ {
			xFullTmp := make([]Num, len(xx))
			xHalfTmp := make([]Num, len(xx))

			copy(xHalfTmp, xx)

			kkFull = method(xx, T, h, dxdt)

			for i, k := range kkFull {
				xFullTmp[i] = xx[i] + k
			}

			var kkHalf []Num

			for halfs := 0; halfs <= 1; halfs++ {
				kkHalf = method(xHalfTmp, T, h/2, dxdt)

				for i, k := range kkHalf {
					xHalfTmp[i] += k
				}
			}

			q := quality(xFullTmp, xHalfTmp, h)

			// store h as the used value
			H = h

			if h < hmin {
				break
			} else if q > 0.005 {
				h /= 2
			} else if q < 0.0005 {
				h *= 2
				break
			} else {
				break
			}
		}

		x := make([]Num, len(xx))
		copy(x, xx)

		r = append(r, Result{T, x})

		T += H

		for i, k := range kkFull {
			xx[i] += k
		}
	}

	return r
}

// quality compares results with h and h/2
func quality(xFull []Num, xHalf []Num, h Num) Num {
	var q Num

	for i, full := range xFull {
		var c Num

		if diff := full - xHalf[i]; diff >= 0 {
			c = diff / h
		} else {
			c = -diff / h
		}

		if c > q {
			q = c
		}
	}

	return q
}

// Euler integration method (deriv a start of interval)
func Euler(xx []Num, t, h Num, dxdt []Ode) (kk []Num) {
	dd := make([]Num, len(xx))
	for i, f := range dxdt {
		dd[i] = f(xx, t)
	}

	kk = make([]Num, len(xx))
	for i, d := range dd {
		kk[i] = h * d
	}

	return kk
}

// MidPoint integration method (derive half way interval)
func MidPoint(xx []Num, t, h Num, dxdt []Ode) (kk []Num) {
	dd := make([]Num, len(xx))

	for i, f := range dxdt {
		dd[i] = f(xx, t)
	}

	xx2 := make([]Num, len(xx))

	for i, x := range xx {
		xx2[i] = x + dd[i]*h/2
	}

	for i, f := range dxdt {
		dd[i] = f(xx2, t+h/2)
	}

	kk = make([]Num, len(xx))
	for i, d := range dd {
		kk[i] = h * d
	}

	return kk
}

// Rk4 Runge Kutta integration method (weighted avarage deriv)
func Rk4(xx0 []Num, t, h Num, dxdt []Ode) (kk []Num) {
	dd0 := make([]Num, len(xx0))
	dd1 := make([]Num, len(xx0))
	dd2 := make([]Num, len(xx0))
	dd3 := make([]Num, len(xx0))

	xxNext := make([]Num, len(xx0))

	for i, f := range dxdt {
		dd0[i] = f(xx0, t)
		xxNext[i] = xx0[i] + h/2*dd0[i]
	}

	for i, f := range dxdt {
		dd1[i] = f(xxNext, t+h/2)
		xxNext[i] = xx0[i] + h/2*dd1[i]
	}

	for i, f := range dxdt {
		dd2[i] = f(xxNext, t+h/2)
		xxNext[i] = xx0[i] + h*dd2[i]
	}

	for i, f := range dxdt {
		dd3[i] = f(xxNext, t+h)
	}

	kk = make([]Num, len(xxNext))

	for i := range xx0 {
		kk[i] = h / 6 * (dd0[i] + 2*dd1[i] + 2*dd2[i] + dd3[i])
	}

	return kk
}
