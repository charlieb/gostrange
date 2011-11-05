package main

import (
	"sdl"
	"fmt"
	"rand"
	"math"
	"time"
	"gl"
)

func feq(x float64, y float64, within int) bool {
	z := x
	for i := 0; i < within; i++ {
		z = math.Nextafter(z, y)
		if y == z { return true }
	}
	return false
}

func fgteq(x float64, y float64, within int) bool {
	return feq(x, y, within) || x > y
}

func flteq(x float64, y float64, within int) bool {
	return feq(x, y, within) || x < y
}


type F2D func(float64, float64) (float64, float64)

func lyapunov_exponent2d(fn F2D, startx, starty float64) float64 {
	const (
		ignore int = 500
	)
	var (
		x1, y1, x2, y2, x3, y3, d0, d1, ltot float64
		i, its int
	)
	// Start with any initial condition in the basin of attraction.
	x1 = startx
	y1 = starty

	// Iterate until the orbit is on the attractor.
	for i = 0; i < ignore; i++ { x1, y1 = fn(x1, y1) }

	// Select (almost any) nearby point (separated by d0).
	// it's about 1000 MIN_FLOATs away from d
	// sqrt 1000 = 31
	x2, y2 = x1, y1
	for i = 0; i < 31; i++ {
		x2 = math.Nextafter(x2, 1000)
		y2 = math.Nextafter(y2, 1000)
	}

	// Get a d0 of the right scale
	d0 = x1 - x2

	//fmt.Println(x1, y1, " - ", x2, y2, x1 - x2, y1 - y2)

	for its = 0; its < 1000; its++ {
		if math.Fabs(x1) > 10000 || math.Fabs(y1) > 10000 ||
			math.IsNaN(x1) || math.IsNaN(x2) {	
			return 0 
		}
		// Advance both orbits one iteration and calculate the new separation d1.
		//fmt.Println("L1", its, x1, y1)
		//fmt.Println("L2", its, x2, y2)
		x1, y1 = fn(x1, y1)
		x3, y3 = fn(x2, y2)

		if feq(x1, x3, 2) && feq(y1, y3, 2) {
//			its--
			continue
		}

		d1 = math.Sqrt(math.Pow(x3 - x1, 2) + math.Pow(y3 - y1, 2))

		// Evaluate log |d1/d0| in any convenient base.
		ltot += math.Log(math.Fabs(d1 / d0))

		// Readjust one orbit so its separation is d0 in the same direction as d1.
		x2 = x1 + d0 * (x3 - x1) / d1
		y2 = y1 + d0 * (y3 - y1) / d1
	}
	// Repeat many times and calculate the average 
	return ltot / float64(its)
}

var (
	//rot float64 = 45.0
	xoff, yoff,	zoff float64 = 0,0,0
	scale float64 = 0.5
	draw_reticule bool = true
)

func plot_vertices(fn func()) {
	gl.Enable(gl.BLEND)
	gl.Enable(gl.POINT_SMOOTH)
	gl.BlendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)

	gl.PointSize(1.0)
	gl.LoadIdentity()
/*
	gl.Rotated(rot, 1.0, 0.0, 1.0)
  gl.Rotated(rot, 1.0, 0.0, 1.0)
*/
	//rot += 0.2
	gl.Translated(xoff, yoff, zoff)
  gl.Scaled(scale, scale, scale)
	gl.Begin(gl.POINTS)
	
	if draw_reticule {
		gl.Color4d(1,0,0, 1.0)
		gl.Vertex2d(0.1,0.1)
		gl.Vertex2d(0.1,-0.1)
		gl.Vertex2d(-0.1,-0.1)
		gl.Vertex2d(-0.1,0.1)
	}

	gl.Color4d(1,1,1, 0.25)

	fn()

	gl.End()

}

func iterate_plot2d(fn F2D, startx, starty float64, disguard, plot int) {
	var (
		plot_fn = func() {
			var (
				x, y float64 = startx, starty
				i int
			)
			
			for i = 0; i < disguard; i++ { x, y = fn(x, y) }
			for i = 0; i < plot; i++ {
				x, y = fn(x, y)
				gl.Vertex2d(x, y)
			}
		}
	)
	plot_vertices(plot_fn)
}


func plot2d(fn func(float64, float64, float64) (float64, float64), 
	disguard, plot int, 
	from, to, step float64) {
	
	var (
		plot_fn = func() {
			var (
				r, x, y float64
				i int
			)
			//xoff += 0.01
			//yoff += 0.01
			
			for r = from; r < to; r += step {
				y = rand.Float64()
				x = rand.Float64()
				//y = 0.1
				for i = 0; i < disguard; i++ { x, y = fn(x, y, r) }
				for i = 0; i < plot; i++ {
					x, y = fn(x, y, r)
					gl.Vertex2d(x, y)
				}
			}
		}
	)

	plot_vertices(plot_fn)
}


func iterate_plot1d(fn func(float64, float64) (float64), 
	disguard, plot int, 
	from, to, step float64) {

	var (
		plot_fn = func() {
			var (
				x, y float64
				i int
			)
			
			//xoff += 0.01
			//yoff += 0.01
			for x = from; x < to; x += step {
				y = rand.Float64()
				//y = 0.1
				for i = 0; i < disguard; i++ { y = fn(y, x) }
				for i = 0; i < plot; i++ {
					y = fn(y, x)
					gl.Vertex2d(x-from,y)
				}
			}
		}
	)

	plot_vertices(plot_fn)
}

type Coeffs2DQuad [12]float64

func make_quadric_map_fn(coeffs *Coeffs2DQuad) F2D {

	return func(x, y float64) (float64, float64) {
		/*
		 Xn+1 = a1 + a2Xn + a3Xn2 + a4XnYn + a5Yn + a6Yn2
		 Yn+1 = a7 + a8Xn + a9Xn2 + a10XnYn + a11Yn + a12Yn2
		 */
		return coeffs[0] + 
			x * coeffs[1] + 
			x * x * coeffs[2] + 
			x * y * coeffs[3] + 
			y * coeffs[4] + 
			y * y * coeffs[5],

			coeffs[6] + 
			x * coeffs[7] + 
			x * x * coeffs[8] + 
			x * y * coeffs[9] + 
			y * coeffs[10] + 
			y * y * coeffs[11]
	}
}

func find_quadric_map_with_L(min, max float64) (cp *Coeffs2DQuad, startx, starty float64) {
	var (
		L float64 = 0
		rejected int = -1
		coeffs Coeffs2DQuad
	)
	startx, starty = rand.Float64(), rand.Float64()

	for L < min || L > max {
		rejected++
		for i := range coeffs {
			coeffs[i] = 1.2 - 2.4 * rand.Float64()
		}
		L = lyapunov_exponent2d(make_quadric_map_fn(&coeffs), startx, starty)
	}
	fmt.Println("Rejected:", rejected)
	fmt.Println("L:", L)
	fmt.Println(coeffs)

	cp = &coeffs

	return
}

func make_henon_fn() F2D {
	//  a1 = 1 a3 = -1.4 a5 = 0.3 a8 = 1 
	coeffs := Coeffs2DQuad{1, 0, -1.4, 0, 0.3, 0, 0, 1, 0,0,0,0}
	return make_quadric_map_fn(&coeffs)
}

func testPlot() {
	var (
		iterations int64 = 0
		start = time.Nanoseconds()
		split int64 = start
		total int64

		new_attractor bool = true

		coeffs *Coeffs2DQuad
		offsets, offset_coeffs Coeffs2DQuad
		startx, starty float64
	)

	for i := range offsets { offsets[i] += rand.Float64() }

	for handleEvents(&new_attractor) {

		if new_attractor { 
			coeffs, startx, starty = find_quadric_map_with_L(0.01, 0.2) 
			new_attractor = false
		}
	
		//iterate_plot1d(itfn1d, 500, 100, 2, 4, 0.001)
		//plot2d(itfn2d, 500, 100, 2, 4, 0.001)
		iterate_plot2d(make_quadric_map_fn(&offset_coeffs), startx, starty, 500, 100000)
		for i := range offsets { offsets[i] += 0.05 + 0.001 * rand.Float64() }
		offset_coeffs = *coeffs
		for i := range offsets { 
			offset_coeffs[i] = coeffs[i] + 0.05 * math.Sin(offsets[i]) 
		}
		sdl.GL_SwapBuffers()
		gl.Clear(gl.COLOR_BUFFER_BIT)

		
		iterations++
		if time.Nanoseconds() - split > 10e9 {
			total = time.Nanoseconds() - start
			fmt.Println(iterations, "iterations in", total / 1e9, "gives", 
				iterations * 1e9 / total, "fps")
			split = time.Nanoseconds()
		}
	}
}

/**************************************************/

func initScreen() {
	sdl.Init(sdl.INIT_VIDEO)
	const (
		resx int = 640
		resy int = 480
	)
	var (
		screen = sdl.SetVideoMode(resx, resy, 16, sdl.OPENGL)
	)

	if screen == nil {
		sdl.Quit()
		panic("Couldn't set GL video mode: " + sdl.GetError() + "\n")
	}

	if gl.Init() != 0 {
		panic("gl error")	
	}

	gl.MatrixMode(gl.PROJECTION)

	gl.Viewport(0, 0, int(screen.W), int(screen.H))
	gl.LoadIdentity()
	gl.Ortho(0, float64(screen.W), float64(screen.H), 0, -1.0, 1.0)

	gl.ClearColor(0, 0, 0, 0)
	gl.Clear(gl.COLOR_BUFFER_BIT)
}

// sdl.GetModState() doesn't work properly so we store state here :(
var (
	mod = sdl.KMOD_NONE
)
func handleEvents(new_attractor *bool) bool {
	for ev := sdl.PollEvent(); ev != nil; ev = sdl.PollEvent() {
		switch e := ev.(type) {
		case *sdl.QuitEvent:
			return false
		case *sdl.KeyboardEvent:
			//fmt.Println(e.Keysym.Sym, ": ", sdl.GetKeyName(sdl.Key(e.Keysym.Sym)))
			switch e.Type {
			case sdl.KEYUP:
				switch sdl.Key(e.Keysym.Sym) {
				case sdl.K_LCTRL: mod = sdl.KMOD_NONE
				}
			case sdl.KEYDOWN:
				switch mod {
				case sdl.KMOD_NONE:
					switch sdl.Key(e.Keysym.Sym) {
					case sdl.K_LCTRL:  mod = sdl.KMOD_LCTRL
					case sdl.K_ESCAPE: return false
					case sdl.K_UP:        yoff += 0.05
					case sdl.K_DOWN:      yoff -= 0.05
					case sdl.K_LEFT:      xoff -= 0.05
					case sdl.K_RIGHT:     xoff += 0.05
					case sdl.K_n:         *new_attractor = true
					default:
					}
				case sdl.KMOD_LCTRL:
					fmt.Println("CTRL")
					switch sdl.Key(e.Keysym.Sym) {
					case sdl.K_LEFT:  scale -= 0.05
					case sdl.K_RIGHT: scale += 0.05
					default:
					}
				}
			}
		default:
		}
	}
	return true
}


/************************/

func main() {
	//test()

	initScreen()
	//testCircleIntersect()
	testPlot()

}