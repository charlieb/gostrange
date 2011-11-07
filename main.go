package main

import (
	"sdl"
  "fmt"
  "rand"
  "math"
  "time"
  "gl"

	"runtime/pprof"
	"flag"
	"os"
	"log"
	"strconv"
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

  for its = 0; its < 10000; its++ {
    if math.Fabs(x1) > 10000 || math.Fabs(y1) > 10000 ||
			math.Fabs(x2) > 10000 || math.Fabs(y2) > 10000 ||
      math.IsNaN(x1) || math.IsNaN(y1) ||
			math.IsNaN(x1) || math.IsNaN(y2) {  
      return 0 
    }
    // Advance both orbits one iteration and calculate the new separation d1.
    //fmt.Println("L1", its, x1, y1)
    //fmt.Println("L2", its, x2, y2)
    x1, y1 = fn(x1, y1)
    x3, y3 = fn(x2, y2)

    if feq(x1, x3, 2) && feq(y1, y3, 2) {
//      its--
      continue
    }

    d1 = math.Sqrt(math.Pow(x3 - x1, 2) + math.Pow(y3 - y1, 2))

    // Evaluate log |d1/d0| in any convenient base.
    ltot += math.Log2(math.Fabs(d1 / d0))

    // Readjust one orbit so its separation is d0 in the same direction as d1.
    x2 = x1 + d0 * (x3 - x1) / d1
    y2 = y1 + d0 * (y3 - y1) / d1
  }
  // Repeat many times and calculate the average 
  return ltot / float64(its)
}

var (
  //rot float64 = 45.0
  xoff, yoff, zoff float64 = 0,0,0
  scale float64 = 0.5
)

func make_2D_plot_fn(fn F2D, startx, starty float64, disguard, plot int) func() {
	return func() {
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
}

func generate_list(fn func()) {
  gl.Begin(gl.POINTS)
  gl.Color4d(1,1,1, 0.25)
  fn()
  gl.End()
}

func plot_list(list uint) {
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
  gl.Scaled(scale, scale, scale)
  gl.Translated(xoff, yoff, zoff)
	gl.CallList(list)
	gl.Flush()
}

func ncoeffs(order int) int { return (order + 1)*(order + 2) }

func multiply_out_terms(order int, x, y float64, mults []float64) {
		
	var coeff int = -1

  for i := 0; i < ncoeffs(order) / 2; i++ { mults[i] = 1 }
	for i := 0; i <= order; i++ {
		for j := 0; j <= order - i; j++ {
			coeff++
			//fmt.Print(coeff, ": ")
			for xm := 0; xm < i; xm++ {
				//fmt.Print("x")
				mults[coeff] *= x
			}
			for ym := 0; ym < j; ym++ {
				//fmt.Print("y")
				mults[coeff] *= y
			}
			//fmt.Print("\n")
		}
	}    
}

func apply_quadric_coeffs(x, y float64, coeffs []float64) (float64, float64) {
  return coeffs[0] + 
		coeffs[1] * y + 
		coeffs[2] * y * y +
		coeffs[3] * x +
    coeffs[4] * x * y + 
    coeffs[5] * x * x,
	coeffs[6] + 
		coeffs[7] * y + 
		coeffs[8] * y * y +
		coeffs[9] * x +
    coeffs[10] * x * y + 
    coeffs[11] * x * x
}


func apply_cubic_coeffs(x, y float64, coeffs []float64) (float64, float64) {
  return coeffs[0] + 
		coeffs[1] * y + 
		coeffs[2] * y * y +
		coeffs[3] * y * y * y +
		coeffs[4] * x +
    coeffs[5] * x * y + 
    coeffs[6] * x * y * y + 
    coeffs[7] * x * x +
    coeffs[8] * x * x * y +
    coeffs[9] * x * x * x,
	coeffs[10] + 
		coeffs[11] * y + 
		coeffs[12] * y * y +
		coeffs[13] * y * y * y +
		coeffs[14] * x +
    coeffs[15] * x * y + 
    coeffs[16] * x * y * y + 
    coeffs[17] * x * x +
    coeffs[18] * x * x * y +
    coeffs[19] * x * x * x
}

func apply_quartic_coeffs(x, y float64, coeffs []float64) (float64, float64) {
  return coeffs[0] + 
		coeffs[1] * y + 
		coeffs[2] * y * y +
		coeffs[3] * y * y * y +
		coeffs[4] * y * y * y + y +
		coeffs[5] * x +
    coeffs[6] * x * y + 
    coeffs[7] * x * y * y + 
    coeffs[8] * x * y * y * y + 
    coeffs[9] * x * x +
    coeffs[10] * x * x * y +
    coeffs[11] * x * x * y * y +
    coeffs[12] * x * x * x + 
    coeffs[13] * x * x * x * y +
    coeffs[14] * x * x * x * x,
	coeffs[15] + 
		coeffs[16] * y + 
		coeffs[17] * y * y +
		coeffs[18] * y * y * y +
		coeffs[19] * y * y * y + y +
		coeffs[20] * x +
    coeffs[21] * x * y + 
    coeffs[22] * x * y * y + 
    coeffs[23] * x * y * y * y + 
    coeffs[24] * x * x +
    coeffs[25] * x * x * y +
    coeffs[26] * x * x * y * y +
    coeffs[27] * x * x * x + 
    coeffs[28] * x * x * x * y +
    coeffs[29] * x * x * x * x
}

func apply_quintic_coeffs(x, y float64, coeffs []float64) (float64, float64) {
  return coeffs[0] + 
		coeffs[1] * y + 
		coeffs[2] * y * y +
		coeffs[3] * y * y * y +
		coeffs[4] * y * y * y + y +
		coeffs[5] * y * y * y + y * y +
		coeffs[6] * x +
    coeffs[7] * x * y + 
    coeffs[8] * x * y * y + 
    coeffs[9] * x * y * y * y + 
    coeffs[10] * x * y * y * y * y + 
    coeffs[11] * x * x +
    coeffs[12] * x * x * y +
    coeffs[13] * x * x * y * y +
    coeffs[14] * x * x * y * y * y +
    coeffs[15] * x * x * x + 
    coeffs[16] * x * x * x * y +
    coeffs[17] * x * x * x * y * y +
    coeffs[18] * x * x * x * x +
    coeffs[19] * x * x * x * x * y +
    coeffs[20] * x * x * x * x * x,
	coeffs[21] + 
		coeffs[22] * y + 
		coeffs[23] * y * y +
		coeffs[24] * y * y * y +
		coeffs[25] * y * y * y + y +
		coeffs[26] * y * y * y + y * y +
		coeffs[27] * x +
    coeffs[28] * x * y + 
    coeffs[29] * x * y * y + 
    coeffs[30] * x * y * y * y + 
    coeffs[31] * x * y * y * y * y + 
    coeffs[32] * x * x +
    coeffs[33] * x * x * y +
    coeffs[34] * x * x * y * y +
    coeffs[35] * x * x * y * y * y +
    coeffs[36] * x * x * x + 
    coeffs[37] * x * x * x * y +
    coeffs[38] * x * x * x * y * y +
    coeffs[39] * x * x * x * x +
    coeffs[40] * x * x * x * x * y +
    coeffs[41] * x * x * x * x * x 
}

func make_map_fn(order int, coeffs []float64) F2D {

	switch order {
	case 2:	return func(x, y float64) (float64, float64) {
			return apply_quadric_coeffs(x ,y, coeffs) 
		}
	case 3:	return func(x, y float64) (float64, float64) {
			return apply_cubic_coeffs(x ,y, coeffs) 
		}
	case 4:	return func(x, y float64) (float64, float64) {
			return apply_quartic_coeffs(x ,y, coeffs) 
		}
	case 5:	return func(x, y float64) (float64, float64) {
			return apply_quintic_coeffs(x ,y, coeffs) 
		}
	}

	var	(
		nterms = ncoeffs(order) / 2
		term_mults []float64 = make([]float64, nterms)
	)
	
	return func(x, y float64) (float64, float64) {
		var (
			xnext float64 = 0
			ynext float64 = 0
		)
		multiply_out_terms(order, x, y, term_mults)
		for i := 0; i < nterms; i++ {
			xnext += term_mults[i] * coeffs[i]
			ynext += term_mults[i] * coeffs[i + nterms]
		}
		return xnext, ynext
	}
}

func find_map_with_L(order int, min, max float64) (coeffs []float64, startx, starty float64) {
  var (
    L float64 = 0
    rejected int = -1
	)
	coeffs = make([]float64, ncoeffs(order))
  startx, starty = rand.Float64(), rand.Float64()

  for L < min || L > max {
    rejected++
    for i := range coeffs {
      coeffs[i] = 1.2 - 2.4 * rand.Float64()
    }
    L = lyapunov_exponent2d(make_map_fn(order, coeffs), startx, starty)
  }
  fmt.Println("Rejected:", rejected)
  fmt.Println("L:", L)
  fmt.Println(coeffs)

  return
}

func make_henon_fn() F2D {
  //  a1 = 1 a3 = -1.4 a5 = 0.3 a8 = 1 
  coeffs := []float64{1, 0, -1.4, 0, 0.3, 0, 0, 1, 0,0,0,0}
  return make_map_fn(2, coeffs)
}

func testPlot(order int) {
  var (
    iterations int64 = 0
    start = time.Nanoseconds()
    split int64 = start
    total int64

    new_attractor bool = true

    coeffs []float64
    //offsets, offset_coeffs []float64
    startx, starty float64

		attractor = gl.GenLists(1);
  )

	//offsets = make([]float64, ncoeffs(order))
	//offset_coeffs = make([]float64, ncoeffs(order))
	coeffs = make([]float64, ncoeffs(order))

  //for i := range offsets { offsets[i] += rand.Float64() }

  for handleEvents(&new_attractor) {

    if new_attractor { 
      coeffs, startx, starty = find_map_with_L(order, 0.1, 0.4) 
			gl.NewList(attractor, gl.COMPILE);
			generate_list(make_2D_plot_fn(
				make_map_fn(order, coeffs), 
				startx, starty, 500, 1000000))
			gl.EndList();
      new_attractor = false
    }
  
		plot_list(attractor)
    //for i := range offsets { offsets[i] += 0.05 + 0.001 * rand.Float64() }
    //offset_coeffs = coeffs
    //for i := range offsets { 
    //  offset_coeffs[i] = coeffs[i] + 0.05 * math.Sin(offsets[i]) 
    //}
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
	xvel, yvel, zvel, svel float64 = 0,0,0,0
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
        case sdl.K_UP, sdl.K_DOWN:     yvel = 0
        case sdl.K_LEFT, sdl.K_RIGHT:  xvel = 0; svel = 0
        }
      case sdl.KEYDOWN:
        switch mod {
        case sdl.KMOD_NONE:
          switch sdl.Key(e.Keysym.Sym) {
          case sdl.K_LCTRL:  mod = sdl.KMOD_LCTRL
          case sdl.K_ESCAPE: return false
          case sdl.K_UP:        yvel = -0.05 / scale
          case sdl.K_DOWN:      yvel = 0.05 / scale
          case sdl.K_LEFT:      xvel = 0.05 / scale
          case sdl.K_RIGHT:     xvel = -0.05 / scale
          case sdl.K_n:         *new_attractor = true
						case sdl.K_z:        
						xoff, yoff, zoff = 0,0,0
						scale = 0.5
					default:
          }
        case sdl.KMOD_LCTRL:
          switch sdl.Key(e.Keysym.Sym) {
          case sdl.K_LEFT:  svel -= 0.05
          case sdl.K_RIGHT: svel += 0.05
          default:
          }
        }
      }
    default:
    }
  }
	xoff += xvel
	yoff += yvel
	zoff += zvel
	scale += svel
  return true
}


/************************/
var (
	cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")
	orderflag = flag.String("order", "5", "Order of attractors to generate")
)

func main() {
	
  flag.Parse()
  if *cpuprofile != "" {
    f, err := os.Create(*cpuprofile)
    if err != nil {
      log.Fatal(err)
    }
    pprof.StartCPUProfile(f)
    defer pprof.StopCPUProfile()
  }

	initScreen()
  if *orderflag != "" {
		var (
			order int
			err os.Error
		)
		order, err = strconv.Atoi(*orderflag)
		if err != nil {
			fmt.Println(err)
		} else {
			testPlot(order)
		}
	}
}