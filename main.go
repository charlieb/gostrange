package main

import (
	"github.com/banthar/Go-SDL/sdl"
  "fmt"
  "math/rand"
  "math"
  "time"
  "github.com/banthar/gl"

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
    p1 Matrix = MakeMatrix(1, start.Height())
		p2 Matrix = MakeMatrix(1, start.Height())
		p3 Matrix = MakeMatrix(1, start.Height())
		d0, d1, ltot, dim_sum float64
    i, its, dim int
  )

  // Start with any initial condition in the basin of attraction.
	copy(p1.Cols()[0], start.Cols()[0])

  // Iterate until the orbit is on the attractor.
  for i = 0; i < ignore; i++ { p1 = fn(p1) }

  // Select (almost any) nearby point (separated by d0).
  // it's about 1000 MIN_FLOATs away from d
  // sqrt 1000 = 31
	copy(p2.Cols()[0], p1.Cols()[0])
  for i = 0; i < 31; i++ {
		for dim = 0; dim < p2.Height(); dim++ {
			p2.Cols()[0][dim] = math.Nextafter(p2.Cols()[0][dim], 1000)
		}
  }

  // Get a d0 of the right scale
  d0 = p1.Cols()[0][0] - p2.Cols()[0][0]

  //fmt.Println(x1, y1, " - ", x2, y2, x1 - x2, y1 - y2)

  for its = 0; its < 10000; its++ {
		for dim = 0; dim < p2.Height(); dim++ {
			if math.Abs(p1.Cols()[0][dim]) > 10000 || 
				math.Abs(p2.Cols()[0][dim]) > 10000 ||
				math.IsNaN(p1.Cols()[0][dim]) || 
				math.IsNaN(p2.Cols()[0][dim]) {
				return 0 
			}
		}
  
		// Advance both orbits one iteration and calculate the new separation d1.
    //fmt.Println("L1", its, p1[dim], y1)
    //fmt.Println("L2", its, x2, y2)
    p1 = fn(p1)
    p3 = fn(p2)

		for dim = 0; dim < p2.Height(); dim++ {
			if feq(p1.Cols()[0][dim], p3.Cols()[0][dim], 2) { continue }
		}

		dim_sum = 0
		for dim = 0; dim < p2.Height(); dim++ {
			dim_sum += math.Pow(p3.Cols()[0][dim] - p1.Cols()[0][dim], 2)
		}
    d1 = math.Sqrt(dim_sum)

    // Evaluate log |d1/d0| in any convenient base.
    ltot += math.Log2(math.Abs(d1 / d0))

    // Readjust one orbit so its separation is d0 in the same direction as d1.
		for dim = 0; dim < p2.Height(); dim++ {
			p2.Cols()[0][dim] = p1.Cols()[0][dim] + 
				d0 * (p3.Cols()[0][dim] - p1.Cols()[0][dim]) / d1
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
      p Matrix = MakeMatrix(1, start.Height())
      i int
    )
		for i = 0; i < start.Width(); i++ { p.Cols()[i][0] = start.Cols()[i][0] }
    for i = 0; i < disguard; i++ { p = fn(p) }
    for i = 0; i < plot; i++ {
      p = fn(p)
      gl.Vertex3d(p.Cols()[0][0], p.Cols()[1][0], p.Cols()[2][0])
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
			next Matrix = MakeMatrix(1, int(math.Max(3, float64(dimension))))
		)
		for d := 0; d < dimension; d++ {
			next.Cols()[0][d] = coeffs[coeff]; coeff++
			for c := 0; c < nCoeffCombs; c++ {
				product = coeffs[coeff]; coeff++
				for ord := 0; ord < order; ord++ {
					if coeffCombs[c][ord] > 0 {
						product *= p.Cols()[0][coeffCombs[c][ord] - 1]
					}
				}
				next.Cols()[0][d] += product
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
	start = MakeMatrix(1, int(math.Max(3, float64(dimension))))
	for i := 0; i < dimension; i++ { start.Cols()[0][i] = rand.Float64() }

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

func MakePointMatrix(fn F, start Matrix, disguard, plot int) Matrix {
	var (
		res Matrix = MakeMatrix(plot, start.Height())
		p Matrix = MakeMatrix(1, start.Height())
		i int
	)
	copy(p.Cols()[0], start.Cols()[0])
  for i = 0; i < disguard; i++ { p = fn(p) }
  for i = 0; i < plot; i++ {
    p = fn(p)
		copy(res.Cols()[i], p.Cols()[0])
  }
	return res
}

var (
	xwrot float64 = 0
)

func testPlot(dimension, order int) {
  var (
    iterations int64 = 0
    start_t = time.Now()
    split time.Time = start_t
    total time.Duration

    new_attractor, redraw bool = true, true
		npoints int = 1e5

    coeffs []float64
    //offsets, offset_coeffs []float64
		start Matrix = MakeMatrix(1, int(math.Max(3, float64(dimension))))
		points, points2 /*, points3*/ Matrix

		//attractor = gl.GenLists(1);
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
			points = MakePointMatrix(makeMapFn(dimension, order, coeffs), start, 500, npoints)
			points2 = MakeMatrix(npoints, points.Height())
			//points3 = MakeMatrix(npoints, points.Height())
			redraw = false
			fmt.Println("Redraw", npoints, "points")
    }

		//points2 = points
		//RotationXW(xwrot).Apply(points, points2)
		//RotationYW(xwrot + 0.25).Apply(points2, points3)
		//RotationZW(xwrot + 0.50).Apply(points3, points2)

		if dimension == 4 {
			RotationZW(xwrot).Apply(points, points2)
			xwrot += 0.05
			ApplyW(points2, points2)
		} else {
			points2 = points
		}

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

		gl.Color4d(1,1,1, 0.25)

		gl.EnableClientState(gl.VERTEX_ARRAY)
		if dimension > 3 {
			gl.VertexPointer(3, (dimension - 3) * 32, points2.FlatCols())
		} else {
			gl.VertexPointer(3, 0, points2.FlatCols())
		}
		gl.DrawArrays(gl.POINTS, 0, points2.Width())
		gl.DisableClientState(gl.VERTEX_ARRAY)

    sdl.GL_SwapBuffers()
    gl.Clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT)

    iterations++
    if time.Since(split).Seconds() > 1.0 {
      total = time.Since(start_t)
      fmt.Println(iterations, "iterations in", total.Seconds(), "gives",
        iterations * 1e9 / total.Nanoseconds(), "fps")
      split = time.Now()
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
		err error
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
