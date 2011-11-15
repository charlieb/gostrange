package main

import (
	//"fmt"
	"math"
)

type Matrix [][]float64

/*

*/
func MakeMatrix(x, y int) Matrix { 
	var ( m Matrix = make([][]float64, x) )
	for i := 0; i < x; i++ { m[i] = make([]float64, y) }
	return m
}

func (m Matrix) Mult(n, mn Matrix) {
	var (
		tot float64 = 0
	)

	/*

	 [ j1i1 j2i1 j3i1 ]   [ j1k1 j1k2 j1k3 ]   [ i1k1 i1k2 i1k3 ]
	 [ j1i2 j2i2 j3i3 ] * [ j2k1 j2k2 j2k3 ] = [ i2k1 i2k2 i2k3 ]
	 [ j1i3 j2i3 j3i4 ]   [ j3k1 j3k2 j3k3 ]   [ i3k1 i3k2 i3k3 ]
	 
	 */

	for i := 0; i < len(m); i++ {
		for k := 0; k < len(n[0]); k++ {
			tot = 0
			for j := 0; j < len(m); j++ {
				tot += m[i][j] * n[j][k]
			}
			mn[i][k] = tot
		}
	}
}

func Identity(size int) (m Matrix) {
	m = make([][]float64, size)
	for i := 0; i < size; i++ {	
		m[i] = make([]float64, size) 
		m[i][i] = 1
	}
	return
}

func Vector(size int) (m Matrix) {
	m = make([][]float64, 1)
	m[0] = make([]float64, size) 
	return
}

func ApplyW(m, mw Matrix) {
	for i := 0; i < len(m[0]); i++ {
		if m[3][i] > 0.3 {m[3][i] = 0.3}
		mw[i][0] = m[i][0] - m[i][0] * m[i][3]
		mw[i][1] = m[i][1] - m[i][1] * m[i][3]
		mw[i][2] = m[i][2] - m[i][2] * m[i][3]
		mw[i][3] = m[i][3] - m[i][3] * m[i][3]
	}
}

func RotationXY(angle float64) (m Matrix) {
	m = Identity(4)
	m[0][0] = math.Cos(angle)
	m[0][1] = math.Sin(angle)
	m[1][0] = -math.Sin(angle)
	m[1][1] = math.Cos(angle)
	return
}

func RotationXZ(angle float64) (m Matrix) {
	m = Identity(4)
	m[0][0] = math.Cos(angle)
	m[0][2] = -math.Sin(angle)
	m[2][0] = math.Sin(angle)
	m[2][2] = math.Cos(angle)
	return
}


func RotationYZ(angle float64) (m Matrix) {
	m = Identity(4)
	m[1][1] = math.Cos(angle)
	m[1][2] = math.Sin(angle)
	m[2][1] = -math.Sin(angle)
	m[2][2] = math.Cos(angle)
	return
}

func RotationXW(angle float64) (m Matrix) {
	m = Identity(4)
	m[0][0] = math.Cos(angle)
	m[0][3] = math.Sin(angle)
	m[3][0] = -math.Sin(angle)
	m[3][3] = math.Cos(angle)
	return
}

func RotationYW(angle float64) (m Matrix) {
	m = Identity(4)
	m[1][1] = math.Cos(angle)
	m[1][3] = -math.Sin(angle)
	m[3][1] = math.Sin(angle)
	m[3][3] = math.Cos(angle)
	return
}

func RotationZW(angle float64) (m Matrix) {
	m = Identity(4)
	m[2][2] = math.Cos(angle)
	m[2][3] = -math.Sin(angle)
	m[3][2] = math.Sin(angle)
	m[3][3] = math.Cos(angle)
	return
}