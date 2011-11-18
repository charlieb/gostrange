package main

import (
	//"fmt"
	"math"
)

type Matrix interface {
	MakeFrom([][]float64)

	Cols() [][]float64
	FlatCols() []float64

	Width() int
	Height() int

	Apply(n, mn Matrix)
	Transpose()
}

type matrix_s struct {
	data []float64
	cols [][]float64
}

func (m *matrix_s) Transpose() {
	mt := makeMatrixS(m.Height(), m.Width())
	for i := 0; i < m.Width(); i++ {
		for j := 0; j < m.Height(); j++ {
			mt.cols[j][i] = m.cols[i][j]
		}
	}
	m.data = mt.data
	m.cols = mt.cols
}

func (m *matrix_s) Width() int { return len(m.cols) }
func (m *matrix_s) Height() int { return len(m.cols[0]) }

func (m *matrix_s) Cols() [][]float64 { return m.cols }
func (m *matrix_s) FlatCols() []float64 { return m.data }

func (m *matrix_s) MakeFrom(c [][]float64) {
	new_m := makeMatrixS(len(c), len(c[0]))
	for i := 0; i < len(c); i++ {
		for j := 0; j < len(c[0]); j++ {
			new_m.cols[i][j] = c[i][j]
		}
	}
	m.data = new_m.data
	m.cols = new_m.cols
}


func makeMatrixS(x, y int) *matrix_s {
	var m matrix_s = matrix_s{ make([]float64, x * y), make([][]float64, x) }
	for i := 0; i < x; i++ { 
		m.cols[i] = m.data[i*y : (i+1)*y]
	}
	return &m
}
func MakeMatrix(x, y int) Matrix { return makeMatrixS(x,y) }

func (m *matrix_s) Apply(n, mn Matrix) {
	var (
		tot float64 = 0
		nr, mr, mnr [][]float64 = n.Cols(), m.Cols(), mn.Cols()
		nx int = n.Width()
		my int = m.Height()
	)

	/*

	 j1i1 j2i1 j3i1   i1k1 i1k2 i1k3   j1k1 j1k2 j1k3
	 j1i2 j2i2 j3i3 * i2k1 i2k2 i2k3 = j2k1 j2k2 j2k3
	 j1i3 j2i3 j3i4   i3k1 i3k2 i3k3   j3k1 j3k2 j3k3
	 
	 */


	for i := 0; i < nx; i++ {     	 // for each column in N
		for k := 0; k < my; k++ {    	 //  for each row in M 
			tot = 0                      //   reset value in MN
			for j := 0; j < my; j++ {    //   for each row in N (= col in M)
				tot += nr[i][j] * mr[j][k] //     multiply M column by N row
			}
			mnr[i][k] = tot              //   final location is M column by N row
		}
	}
}

func Identity(size int) (m Matrix) {
	m = MakeMatrix(size, size)
	for i,c := range m.Cols() { c[i] = 1 }
	return
}

func ApplyW(m, mw Matrix) {
	var (
		mc, mwc [][]float64 = m.Cols(), mw.Cols()
		mx = m.Width()
	)

	for i := 0; i < mx; i++ {
		//if mc[i][3] > 0.3 {mc[i][3] = 0.3}
		mwc[i][0] = mc[i][0] - mc[i][0] * mc[i][3]
		mwc[i][1] = mc[i][1] - mc[i][1] * mc[i][3]
		mwc[i][2] = mc[i][2] - mc[i][2] * mc[i][3]
		mwc[i][3] = mc[i][3] - mc[i][3] * mc[i][3]
	}
}

func RotationXY(angle float64) Matrix {
	var (
		m Matrix = Identity(4)
		c [][]float64 = m.Cols()
	)
	
	c[0][0] = math.Cos(angle)
	c[0][1] = math.Sin(angle)
	c[1][0] = -math.Sin(angle)
	c[1][1] = math.Cos(angle)
	return m
}

func RotationXZ(angle float64) Matrix {
	var (
		m Matrix = Identity(4)
		c [][]float64 = m.Cols()
	)
	
	c[0][0] = math.Cos(angle)
	c[0][2] = -math.Sin(angle)
	c[2][0] = math.Sin(angle)
	c[2][2] = math.Cos(angle)
	return m
}


func RotationYZ(angle float64) Matrix {
	var (
		m Matrix = Identity(4)
		c [][]float64 = m.Cols()
	)
	
	c[1][1] = math.Cos(angle)
	c[1][2] = math.Sin(angle)
	c[2][1] = -math.Sin(angle)
	c[2][2] = math.Cos(angle)
	return m
}

func RotationXW(angle float64) Matrix {
	var (
		m Matrix = Identity(4)
		c [][]float64 = m.Cols()
	)
	
	c[0][0] = math.Cos(angle)
	c[0][3] = math.Sin(angle)
	c[3][0] = -math.Sin(angle)
	c[3][3] = math.Cos(angle)
	return m
}

func RotationYW(angle float64) Matrix {
	var (
		m Matrix = Identity(4)
		c [][]float64 = m.Cols()
	)
	
	c[1][1] = math.Cos(angle)
	c[1][3] = -math.Sin(angle)
	c[3][1] = math.Sin(angle)
	c[3][3] = math.Cos(angle)
	return m
}

func RotationZW(angle float64) Matrix {
	var (
		m Matrix = Identity(4)
		c [][]float64 = m.Cols()
	)
	
	c[2][2] = math.Cos(angle)
	c[2][3] = -math.Sin(angle)
	c[3][2] = math.Sin(angle)
	c[3][3] = math.Cos(angle)
	return m
}