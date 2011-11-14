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

type F func(Matrix) Matrix

func lyapunov_exponent(fn F, start Matrix) float64 {
  const (
    ignore int = 500
  )
  var (
    p1 Matrix = Vector(len(start))
		p2 Matrix = Vector(len(start))
		p3 Matrix = Vector(len(start))
		d0, d1, ltot, dim_sum float64
    i, its, dim int
  )
  // Start with any initial condition in the basin of attraction.
	for i = 0; i < len(start); i++ { p1[i][0] = start[i][0] }

  // Iterate until the orbit is on the attractor.
  for i = 0; i < ignore; i++ { p1 = fn(p1) }

  // Select (almost any) nearby point (separated by d0).
  // it's about 1000 MIN_FLOATs away from d
  // sqrt 1000 = 31
	for i = 0; i < len(start); i++ { p2[i][0] = p1[i][0] }
  for i = 0; i < 31; i++ {
		for dim = 0; dim < len(p2); dim++ {
			p2[dim][0] = math.Nextafter(p2[dim][0], 1000)
		}
  }

  // Get a d0 of the right scale
  d0 = p1[0][0] - p2[0][0]

  //fmt.Println(x1, y1, " - ", x2, y2, x1 - x2, y1 - y2)

  for its = 0; its < 10000; its++ {
		for dim = 0; dim < len(p2); dim++ {
			if math.Fabs(p1[dim][0]) > 10000 || math.Fabs(p2[dim][0]) > 10000 ||
				math.IsNaN(p1[dim][0]) || math.IsNaN(p2[dim][0]) {
				return 0 
			}
		}
  
		// Advance both orbits one iteration and calculate the new separation d1.
    //fmt.Println("L1", its, p1[dim], y1)
    //fmt.Println("L2", its, x2, y2)
    p1 = fn(p1)
    p3 = fn(p2)

		for dim = 0; dim < len(p2); dim++ {
			if feq(p1[dim][0], p3[dim][0], 2) { continue }
		}

		dim_sum = 0
		for dim = 0; dim < len(p2); dim++ {
			dim_sum += math.Pow(p3[dim][0] - p1[dim][0], 2)
		}
    d1 = math.Sqrt(dim_sum)

    // Evaluate log |d1/d0| in any convenient base.
    ltot += math.Log2(math.Fabs(d1 / d0))

    // Readjust one orbit so its separation is d0 in the same direction as d1.
		for dim = 0; dim < len(p2); dim++ {
			p2[dim][0] = p1[dim][0] + d0 * (p3[dim][0] - p1[dim][0]) / d1
		}
	}
  // Repeat many times and calculate the average 
  return ltot / float64(its)
}

var (
  //rot float64 = 45.0
  xoff, yoff, zoff float64 = 0,0,0
  xrot, yrot, zrot float64 = 0,0,0
  scale float64 = 0.5
)

func make_plot_fn(fn F, start Matrix, disguard, plot int) func() {
	return func() {
    var (
      p Matrix = Vector(len(start))
      i int
    )
		copy(p, start)

    for i = 0; i < disguard; i++ { p = fn(p) }
    for i = 0; i < plot; i++ {
      p = fn(p)
      gl.Vertex3d(p[0][0], p[1][0], p[2][0])
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
  gl.Rotated(xrot, 1,0,0)
  gl.Rotated(yrot, 0,1,0)
  gl.Rotated(zrot, 0,0,1)
  gl.Scaled(scale, scale, scale)
  gl.Translated(xoff, yoff, zoff)
	gl.CallList(list)
	gl.Flush()
}

func nCoeffs(dimension, order int) int {
	switch {
	case dimension == 0 || order == 0: return 0
	case dimension == 1: return order
	case order == 1: return dimension
	default: return 2 * nCoeffs(dimension - 1, order) -
			nCoeffs(dimension - 2, order) +
			nCoeffs(dimension, order - 1) - 
			nCoeffs(dimension - 1, order - 1)
	}
	return 0
}

func coeffCombinations(dimension, order int) [][]int {
	var (
		vals []int = make([]int, order)	
		vals_lt_dimension bool = true
		results [][]int = make([][]int, 0)
		add_result func([]int) = func(new_result []int) {
			new_copy := make([]int, order)
			copy(new_copy, new_result)
			results = append(results, new_copy)
		}
	)
	for vals_lt_dimension {
		vals_lt_dimension = false
		for i := 0; i < order; i++ {
			if vals[i] < dimension {
				vals[i]++
				vals_lt_dimension = true
				add_result(vals)
				break
			} else {
				for j := i + 1; j < order; j++ {
					if vals[j] < dimension {
						vals[j]++
						for k := i; k < j; k++ {
							vals[k] = vals[j]
						}
						add_result(vals)
						i = -1
						break
					}
				}
			}
		}
	}
	return results
}

func makeMapFn(dimension, order int, coeffs []float64) F {
	var (
		coeffCombs [][]int = coeffCombinations(dimension, order)
		nCoeffCombs int = nCoeffs(dimension, order)
	)
	return func(p Matrix) Matrix {
		var ( 
			product float64
			coeff int = 0
			next Matrix = Vector(int(math.Fmax(3, float64(dimension))))
		)
		for d := 0; d < dimension; d++ {
			next[d][0] = coeffs[coeff]; coeff++
			for c := 0; c < nCoeffCombs; c++ {
				product = coeffs[coeff]; coeff++
				for ord := 0; ord < order; ord++ {
					if coeffCombs[c][ord] > 0 {
						product *= p[coeffCombs[c][ord] - 1][0]
					}
				}
				next[d][0] += product
			}
		}
		return next
	}
}


func find_map_with_L(dimension, order int, min, max float64) (coeffs []float64, start Matrix) {
  var (
    L float64 = 0
    rejected int = -1
	)
	coeffs = make([]float64, dimension + nCoeffs(order, dimension) * dimension)
	start = Vector(int(math.Fmax(3, float64(dimension))))
	for i := 0; i < dimension; i++ { start[i][0] = rand.Float64() }

  for L < min || L > max {
    rejected++
    for i := range coeffs {
      coeffs[i] = 1.2 - 2.4 * rand.Float64()
    }
    L = lyapunov_exponent(makeMapFn(dimension, order, coeffs), start)
  }
  fmt.Println("Rejected:", rejected)
  fmt.Println("L:", L)
  fmt.Println(coeffs)

  return
}

func testPlot(dimension, order int) {
  var (
    iterations int64 = 0
    start_t = time.Nanoseconds()
    split int64 = start_t
    total int64

    new_attractor, redraw bool = true, true
		npoints int = 1e6

    coeffs []float64
    //offsets, offset_coeffs []float64
		start Matrix = Vector(int(math.Fmax(3, float64(dimension))))

		attractor = gl.GenLists(1);
  )

	//offsets = make([]float64, ncoeffs(order))
	//offset_coeffs = make([]float64, ncoeffs(order))
	//coeffs = make([]float64, dimension + nCoeffs(order, dimension) * dimension)

  //for i := range offsets { offsets[i] += rand.Float64() }

  for handleEvents(&new_attractor, &redraw, &npoints) {

    if new_attractor { 
      coeffs, start = find_map_with_L(dimension, order, 0.1, 0.4)
			redraw = true
      new_attractor = false
		}
		if redraw {
			gl.NewList(attractor, gl.COMPILE);
			generate_list(
				make_plot_fn(makeMapFn(dimension, order, coeffs), start, 500, npoints))
			gl.EndList();
			redraw = false
			
			fmt.Println("Redraw", npoints, "points")
    }
  
		plot_list(attractor)
    //for i := range offsets { offsets[i] += 0.05 + 0.001 * rand.Float64() }
    //offset_coeffs = coeffs
    //for i := range offsets { 
    //  offset_coeffs[i] = coeffs[i] + 0.05 * math.Sin(offsets[i]) 
    //}
    sdl.GL_SwapBuffers()
    gl.Clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT)

    
    iterations++
    if time.Nanoseconds() - split > 10e9 {
      total = time.Nanoseconds() - start_t
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
	//gl.MatrixMode(gl.MODELVIEW)
	
  gl.Viewport(0, 0, int(screen.W), int(screen.H))
  gl.LoadIdentity()
  gl.Ortho(0, float64(screen.W), float64(screen.H), 0, -1.0, 1.0)
	//gl.DepthRange(-1, 1)

  gl.ClearColor(0, 0, 0, 0)
  gl.Clear(gl.COLOR_BUFFER_BIT)
}

// sdl.GetModState() doesn't work properly so we store state here :(
var (
  mod = sdl.KMOD_NONE
	xvel, yvel, zvel, xrotvel, yrotvel, zrotvel, svel float64 = 0,0,0,0,0,0,0
)
func handleEvents(new_attractor, redraw *bool, npoints *int) bool {
  for ev := sdl.PollEvent(); ev != nil; ev = sdl.PollEvent() {
    switch e := ev.(type) {
    case *sdl.QuitEvent:
      return false
    case *sdl.KeyboardEvent:
      //fmt.Println(sdl.Key(e.Keysym.Sym), ":", sdl.GetKeyName(sdl.Key(e.Keysym.Sym)))
      switch e.Type {
      case sdl.KEYUP:
        switch sdl.Key(e.Keysym.Sym) {
        case sdl.K_LCTRL: mod = sdl.KMOD_NONE
        case sdl.K_UP, sdl.K_DOWN:     yvel = 0
        case sdl.K_LEFT, sdl.K_RIGHT:  xvel = 0; svel = 0
        case sdl.K_KP4, sdl.K_KP6:     yrotvel = 0
        case sdl.K_KP8, sdl.K_KP2:     xrotvel = 0
        case sdl.K_KP7, sdl.K_KP9:     zrotvel = 0
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
          case sdl.K_KP4:       yrotvel = -0.5
          case sdl.K_KP6:       yrotvel = +0.5
          case sdl.K_KP8:       xrotvel = -0.5
          case sdl.K_KP2:       xrotvel = +0.5
          case sdl.K_KP7:       zrotvel = -0.5
          case sdl.K_KP9:       zrotvel = +0.5
          case sdl.K_n:         *new_attractor = true
						case sdl.K_z:        
						xoff, yoff, zoff, xrot, yrot, zrot = 0,0,0,0,0,0
						scale = 0.5
					case sdl.K_COMMA:
						*npoints -= 1e5
						*redraw = true
					case sdl.K_PERIOD:
						*npoints += 1e5
						*redraw = true
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
	xrot += xrotvel
	yrot += yrotvel
	zrot += zrotvel
	scale += svel
  return true
}


/************************/
var (
	cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")
	orderflag = flag.String("order", "5", "Order of attractors to generate")
	dimensionflag = flag.String("dimension", "2", "Dimension of attractors to generate")
)

func main() {

	var (
		order, dimension int
		err os.Error
	)

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
		order, err = strconv.Atoi(*orderflag)
		if err != nil {
			fmt.Println(err)
			return
		}
	}

	if *dimensionflag != "" {
		dimension, err = strconv.Atoi(*dimensionflag)
		if err != nil {
			fmt.Println(err)
			return
		}
	}

	testPlot(dimension, order)
}