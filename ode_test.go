package ode

import "fmt"

func Example_Rk4() {
	var k, m, b Num

	k = 1
	m = 1
	b = 0.4

	dxdt := func(xx []Num, t Num) Num { v := xx[1]; return v }
	dvdt := func(xx []Num, t Num) Num { x, v := xx[0], xx[1]; return -k*x/m - b*v/m }

	odes := []Ode{dxdt, dvdt}

	fmt.Printf("%9s %9s %9s %9s %9s\n", "t", "x", "v", "x'", "v'")

	// select one of the steppers and one of the integration methods
	result := AdaptiveStep(Rk4, odes, []Num{-0.5, 0}, 0, 15, 0.01, 0.5)

	d := result[len(result)-1]

	fmt.Printf("%9.3f %9.3f %9.3f %9.3f %9.3f\n",
		d.t, d.xx[0], d.xx[1], odes[0](d.xx, d.t), odes[1](d.xx, d.t))

	// output:
	//         t         x         v        x'        v'
	//    14.875    -0.003     0.016     0.016    -0.004
}
