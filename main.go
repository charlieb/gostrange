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
	"sort"
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

type F func([]float64) []float64

func lyapunov_exponent(fn F, start []float64) float64 {
  const (
    ignore int = 500
  )
  var (
    p1 []float64 = make([]float64, len(start))
		p2 []float64 = make([]float64, len(start))
		p3 []float64 = make([]float64, len(start))
		d0, d1, ltot, dim_sum float64
    i, its, dim int
  )
  // Start with any initial condition in the basin of attraction.
  copy(p1, start)

  // Iterate until the orbit is on the attractor.
  for i = 0; i < ignore; i++ { p1 = fn(p1) }

  // Select (almost any) nearby point (separated by d0).
  // it's about 1000 MIN_FLOATs away from d
  // sqrt 1000 = 31
  copy(p2 ,p1)
  for i = 0; i < 31; i++ {
		for dim = 0; dim < len(p2); dim++ {
			p2[dim] = math.Nextafter(p2[dim], 1000)
		}
  }

  // Get a d0 of the right scale
  d0 = p1[0] - p2[0]

  //fmt.Println(x1, y1, " - ", x2, y2, x1 - x2, y1 - y2)

  for its = 0; its < 10000; its++ {
		for dim = 0; dim < len(p2); dim++ {
			if math.Fabs(p1[dim]) > 10000 || math.Fabs(p2[dim]) > 10000 ||
				math.IsNaN(p1[dim]) || math.IsNaN(p2[dim]) {
				return 0 
			}
		}
  
		// Advance both orbits one iteration and calculate the new separation d1.
    //fmt.Println("L1", its, p1[dim], y1)
    //fmt.Println("L2", its, x2, y2)
    p1 = fn(p1)
    p3 = fn(p2)

		for dim = 0; dim < len(p2); dim++ {
			if feq(p1[dim], p3[dim], 2) { continue }
		}

		dim_sum = 0
		for dim = 0; dim < len(p2); dim++ {
			dim_sum += math.Pow(p3[dim] - p1[dim], 2)
		}
    d1 = math.Sqrt(dim_sum)

    // Evaluate log |d1/d0| in any convenient base.
    ltot += math.Log2(math.Fabs(d1 / d0))

    // Readjust one orbit so its separation is d0 in the same direction as d1.
		for dim = 0; dim < len(p2); dim++ {
			p2[dim] = p1[dim] + d0 * (p3[dim] - p1[dim]) / d1
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

func make_plot_fn(fn F, start []float64, disguard, plot int) func() {
	return func() {
    var (
      p []float64 = make([]float64, len(start))
      i int
    )
		copy(p, start)

    for i = 0; i < disguard; i++ { p = fn(p) }
    for i = 0; i < plot; i++ {
      p = fn(p)
      gl.Vertex3d(p[0], p[1], p[2])
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

func ncoeffs(order, dimension int) int { return (order + 1)*(order + 2) }

func multiply_out_terms(order int, x, y float64, mults []float64) {
		
	var coeff int = -1

  for i := 0; i < ncoeffs(order, 2) / 2; i++ { mults[i] = 1 }
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

func apply_quadric_coeffs(p []float64, coeffs []float64) ([]float64) {
	var (	x, y float64 = p[0], p[1] )
	p[0] = coeffs[0] + 
		coeffs[1] * y + 
		coeffs[2] * y * y +
		coeffs[3] * x +
		coeffs[4] * x * y + 
		coeffs[5] * x * x
	p[1] = coeffs[6] + 
		coeffs[7] * y + 
		coeffs[8] * y * y +
		coeffs[9] * x +
		coeffs[10] * x * y + 
		coeffs[11] * x * x 
	return p
}


func apply_cubic_coeffs(p []float64, coeffs []float64) ([]float64) {
	var (	x, y float64 = p[0], p[1] )
	p[0] = coeffs[0] + 
		coeffs[1] * y + 
		coeffs[2] * y * y +
		coeffs[3] * y * y * y +
		coeffs[4] * x +
    coeffs[5] * x * y + 
    coeffs[6] * x * y * y + 
    coeffs[7] * x * x +
    coeffs[8] * x * x * y +
    coeffs[9] * x * x * x
	p[1] = coeffs[10] + 
		coeffs[11] * y + 
		coeffs[12] * y * y +
		coeffs[13] * y * y * y +
		coeffs[14] * x +
    coeffs[15] * x * y + 
    coeffs[16] * x * y * y + 
    coeffs[17] * x * x +
    coeffs[18] * x * x * y +
    coeffs[19] * x * x * x
	return p
}

func apply_quartic_coeffs(p []float64, coeffs []float64) ([]float64) {
	var (	x, y float64 = p[0], p[1] )
	p[0] = coeffs[0] + 
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
    coeffs[14] * x * x * x * x
	p[1] = coeffs[15] + 
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
	return p
}

func apply_quintic_coeffs(p []float64, coeffs []float64) ([]float64) {
	var (	x, y float64 = p[0], p[1] )
	p[0] = coeffs[0] + 
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
    coeffs[20] * x * x * x * x * x
	p[1] = coeffs[21] + 
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
	return p
}

func make_map_fn(order, dimension int, coeffs []float64) F {

	switch dimension {
	case 2:
		switch order {
		case 2:	return func(p []float64) []float64 {
				return apply_quadric_coeffs(p, coeffs) 
			}
		case 3:	return func(p []float64) []float64 {
				return apply_cubic_coeffs(p, coeffs) 
			}
		case 4:	return func(p []float64) []float64 {
				return apply_quartic_coeffs(p, coeffs) 
			}
		case 5:	return func(p []float64) []float64 {
				return apply_quintic_coeffs(p, coeffs) 
			}
		}
		//case 3
	}

	var	(
		nterms = ncoeffs(order, dimension) / dimension
		term_mults []float64 = make([]float64, nterms)
	)
	
	return func(p []float64) []float64 {
		var (
			xnext float64 = 0
			ynext float64 = 0
		)
		multiply_out_terms(order, p[0], p[1], term_mults)
		for i := 0; i < nterms; i++ {
			xnext += term_mults[i] * coeffs[i]
			ynext += term_mults[i] * coeffs[i + nterms]
		}
		p[0], p[1] = xnext, ynext
		return p
	}
}

func find_map_with_L(order, dimension int, min, max float64) (coeffs []float64, start []float64) {
  var (
    L float64 = 0
    rejected int = -1
	)
	coeffs = make([]float64, ncoeffs(order, dimension))
	start = make([]float64, int(math.Fmax(3, float64(dimension))))
	for i := 0; i < dimension; i++ { start[i] = rand.Float64() }

  for L < min || L > max {
    rejected++
    for i := range coeffs {
      coeffs[i] = 1.2 - 2.4 * rand.Float64()
    }
    L = lyapunov_exponent(make_map_fn(order, dimension, coeffs), start)
  }
  fmt.Println("Rejected:", rejected)
  fmt.Println("L:", L)
  fmt.Println(coeffs)

  return
}

func make_henon_fn() F {
  //  a1 = 1 a3 = -1.4 a5 = 0.3 a8 = 1 
  coeffs := []float64{1, 0, -1.4, 0, 0.3, 0, 0, 1, 0,0,0,0}
  return make_map_fn(2, 2, coeffs)
}

func testPlot(order, dimension int) {
  var (
    iterations int64 = 0
    start_t = time.Nanoseconds()
    split int64 = start_t
    total int64

    new_attractor, redraw bool = true, true
		npoints int = 1e6

    coeffs []float64
    //offsets, offset_coeffs []float64
		start []float64 = make([]float64, int(math.Fmax(3, float64(dimension))))

		attractor = gl.GenLists(1);
  )

	//offsets = make([]float64, ncoeffs(order))
	//offset_coeffs = make([]float64, ncoeffs(order))
	coeffs = make([]float64, ncoeffs(order, 2))

  //for i := range offsets { offsets[i] += rand.Float64() }

  for handleEvents(&new_attractor, &redraw, &npoints) {

    if new_attractor { 
      coeffs, start = find_map_with_L(order, dimension, 0.1, 0.4)
			redraw = true
      new_attractor = false
		}
		if redraw {
			gl.NewList(attractor, gl.COMPILE);
			generate_list(
				make_plot_fn(make_map_fn(order, dimension, coeffs), start, 500, npoints))
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
    gl.Clear(gl.COLOR_BUFFER_BIT)

    
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

  gl.Viewport(0, 0, int(screen.W), int(screen.H))
  gl.LoadIdentity()
  gl.Ortho(0, float64(screen.W), float64(screen.H), 0, -1.0, 1.0)

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
						xoff, yoff, zoff = 0,0,0
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

type sliceSliceInt [][]int
func (ssi sliceSliceInt) Len() int { return len(ssi) }
// Less returns whether the element with index i should sort
// before the element with index j.
func (ssi sliceSliceInt) Less(i, j int) bool {
	for x := range ssi[i] {
		if ssi[i][x] != ssi[j][x] { return ssi[i][x] < ssi[j][x] }
	}
	return true
}
// Swap swaps the elements with indexes i and j.
func (ssi sliceSliceInt) Swap(i, j int) {
	tmp := ssi[i]
	ssi[i] = ssi[j]
	ssi[j] = tmp
}

func allCombinations(dimension, order int) [][]int {
	var (
		result [][]int = make([][]int, 0)
		fn func(int, int, []int) 
		pref = make([]int, 0)
	)

	fn = func(dimension, order int, prefix []int) {
		if order == 0 { 
			result = append(result, make([]int, len(prefix)))
			copy(result[len(result) - 1], prefix)
			sort.IntSlice(result[len(result) - 1]).Sort()
			return 
		}
		for d := 0; d <= dimension; d++ { 
			fn(dimension, order - 1, append(prefix, d))
		} 
	}
	
	fn(dimension, order, pref)
	return result
}

func uniqueCombinations(dimension, order int) [][]int {
	var (
		result [][]int = allCombinations(dimension, order)
		last, end int = 0, len(result)
	)
	sort.Sort(sliceSliceInt(result))

	for i := 1; i < end; i++ {
		for elem := range result[i] {
			if result[i][elem] != result[last][elem] {
				copy(result[last+1:], result[i:])
				end -= i - last - 1
				result = result[:end]
				last += 1
				i = last
				break
			}
		}
	}
		
	return result[:end]
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

	testPlot(order, dimension)
}