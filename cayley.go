package cayley
import (
	"fmt"
	"math"
)
func Conj(x complex128) complex128 {
	return complex(real(x),-1*imag(x))
}
func SumSqr(x complex128) float64 {
	return (real(x)*real(x)) + (imag(x)*imag(x))
}

type CDC struct {
	left,right complex128
}
type CDN struct {
	left,right CayleyDickson
}
type CayleyDickson interface {
	Complex128() complex128
	Real() float64
	Imag() CayleyDickson
	IsZero(...float64) bool
	IsDirection(...float64) bool
	Direction(...float64) int
	Equal(CayleyDickson,...float64) bool
	Resize(int) CayleyDickson
	Translate(complex128) CayleyDickson
	Length() int
	Left() CayleyDickson
	Right() CayleyDickson
	Conj() CayleyDickson
	Neg() CayleyDickson
	SumSqr() float64
	Scale(float64) CayleyDickson
	String(...int) string
	Values() []float64
	Matrix() [][]float64
}

func getDir(dir,sz int) string {
	if dir == 0 {
		return ""
	} else if sz <= 2 {
		return string(rune(dir-1+'i'))
	} else if dir <= 26 {
			return string(rune(dir-1+'a'))
	} else if dir > 52 {
		k := dir-1
		return getDir(k/52,sz)+getDir((k%52)+1,sz)
	} else {
		return string(rune(dir-27+'A'))
	}}
func NewZero(sz int) CayleyDickson {
	if sz == 0 {
		return CDC{0,0}
    } else if sz == 1 {
		return CDN{CDC{0,0},CDC{0,0}}
	} else {
		return CDN{NewZero(sz-1),NewZero(sz-1)}
	}}
func NewScalar(dir int, x complex128) CayleyDickson {
	if dir == 0 {
		return CDC{x,0}
	} else if dir == 1 {
		return CDC{0,x}
	} else {
		_,e := math.Frexp(float64(dir))
		sz := (e-2)/2
		return CDN{NewZero(sz),NewScalar(dir - (1 << (e-1)),x).Resize(sz)}
	}}
func NewPosition(vals []float64) CayleyDickson {
	if k := len(vals); k == 0 {
		return CDC{0,0}
	} else if k == 1 {
		return CDC{complex(vals[0],0),0}
	} else if k == 2 {
		return CDC{complex(vals[0],vals[1]),0}
	} else if k == 3 {
		return CDC{complex(vals[0],vals[1]),complex(vals[2],0)}
	} else if k == 4 {
		return CDC{complex(vals[0],vals[1]),complex(vals[2],vals[3])}
	} else {
		_,e := math.Frexp(float64(k-1));
		q := 2 << (e-1)
		return CDN{NewPosition(vals[:q]),NewPosition(vals[q:]).Resize((e-2)/2)}
	}}
func Add(x,y CayleyDickson, e ...float64) CayleyDickson {
	a,b := x.Length(),y.Length()
	if x.IsZero(e...) {
		return y
	} else if y.IsZero(e...) {
		return x
	} else if a == b {
		if a == 0 {
			return CDC{x.(CDC).left+y.(CDC).left,x.(CDC).right+y.(CDC).right}
		} else {
			return CDN{Add(x.Left(),y.Left(),e...),Add(x.Right(),y.Right(),e...)}
		}
	} else if a < b {
		return CDN{Add(x,y.Left(),e...),y.Right()}
	} else {
		return CDN{Add(x.Left(),y,e...),x.Right()}
	}}
func FMA(x CayleyDickson, y,z float64) CayleyDickson {
	if x.Length() == 0 {
		v := x.(CDC)
		k := complex(y,0)
		return CDC{(v.left*k)+complex(z,0),v.right*k}
	} else {
		return CDN{FMA(x.Left(),y,z),x.Right().Scale(y)}
	}}
func Mul(x,y CayleyDickson,e ...float64) CayleyDickson {
	p,q := x.Length(),y.Length()
	A := func() CayleyDickson {
		return Mul(x.Left(),y.Left())
	}
	B := func() CayleyDickson {
		return Mul(y.Right().Conj(),x.Right()).Neg()
	}
	C := func() CayleyDickson {
		return Mul(y.Right(),x.Left())
	}
	D := func() CayleyDickson {
		return Mul(x.Right(),y.Left().Conj())
	}
	var ab,cd CayleyDickson
	lx,rx := x.Left().IsZero(e...),x.Right().IsZero(e...)
	ly,ry := y.Left().IsZero(e...),y.Right().IsZero(e...)
	switch {
	case (lx && rx) || (ly && ry):
		if p >= q {
			return NewZero(p)
		} else {
			return NewZero(q)
		}
	case q == 0 && y.(CDC).right == 0 && imag(y.(CDC).left) == 0:
		return x.Scale(real(y.(CDC).left))
	case p == 0 && x.(CDC).right == 0 && imag(x.(CDC).left) == 0:
		return y.Scale(real(x.(CDC).left))
	case p == q && p == 0:
		v,w := x.(CDC),y.(CDC)
		a := v.left * w.left
		b := Conj(w.right) * v.right
		c := w.right * v.left
		d := v.right * Conj(w.left)
		return CDC{a-b,c+d}
	case p > q:
		return Mul(x,y.Resize(p))
	case p < q:
		if ry {
			return Mul(x,y.Left())
		} else {
			return Mul(x.Resize(q),y)
		}
	case rx && ry:
		return A().Resize(p)
	case rx || ry:
		if lx || ly {
			ab = NewZero(p-1)
		} else {
			ab = A()
		}
		if rx {
			cd = C()
		} else {
			cd = D()
		}
	case lx || ly:
		ab = B()
		if lx && ly {
			cd = NewZero(p-1)
		} else if ly {
			cd = C()
		} else {
			cd = D()
		}
	default:
		ab = Add(A(), B())
		cd = Add(C(), D())
	}
	return CDN{ab,cd}
}
func Norm(x CayleyDickson) float64 {
	return math.Sqrt(x.SumSqr())
}
func InvMul(x CayleyDickson) CayleyDickson {
	n := Norm(x)
	return x.Conj().Scale(1/(n*n))
}
func Versor(x CayleyDickson) CayleyDickson {
	return x.Scale(1/Norm(x))
}
func Sqrt(x CayleyDickson) CayleyDickson {
	n := Norm(x)
	a := x.Real()
	r := math.Sqrt((n+a)/2)
	s := math.Sqrt((n-a)/2)
	v := x.Imag()
	return FMA(v,s/Norm(v),r)
}
func Exp(x CayleyDickson) CayleyDickson {
	v := x.Imag()
	n := Norm(v)
	q := math.Exp(x.Real())
	return FMA(v,math.Sin(n)*q/n,math.Cos(n)*q)
}
func Log(x CayleyDickson) CayleyDickson {
	n := Norm(x)
	a := x.Real()
	v := x.Imag()
	return FMA(v,math.Acos(a/n)/Norm(v),math.Log(n))
}
func Pow(x,y CayleyDickson) CayleyDickson {
	return Exp(Mul(Log(x),y))
}
func (c CDC) Real() float64 {
	return real(c.left)
}
func (c CDC) Imag() CayleyDickson {
	return CDC{complex(0,imag(c.left)),c.right}
}
func (c CDC) Complex128() complex128 {
	return c.left
}
func (c CDC) IsZero(e ...float64) bool {
	if len(e) == 0 {
		return (c.left == 0) && (c.right == 0)
	}
	return (math.Abs(real(c.left)) <= e[0]) && (math.Abs(real(c.right)) <= e[0]) &&
		   (math.Abs(imag(c.left)) <= e[0]) && (math.Abs(imag(c.right)) <= e[0])
}
func (c CDC) IsDirection(e ...float64) bool {
	t := (len(e) == 0)
	if t && c.left == 0 {
		return (c.right != 0)
	} else if !t && math.Abs(real(c.left)) <= e[0] && math.Abs(imag(c.left)) <= e[0] {
		return (math.Abs(real(c.right)) > e[0]) || (math.Abs(imag(c.right)) > e[0])
	} else if t {
		return (c.right == 0)
	}
	return (math.Abs(real(c.right)) <= e[0]) && (math.Abs(imag(c.right)) <= e[0])
}
func (c CDC) Direction(e ...float64) int {
	t := (len(e) == 0)
	if t && c.left == 0 {
		if c.right == 0 {
			return -1
		} else {
			return 1
		}
	} else if !t && (math.Abs(real(c.left)) <= e[0]) && (math.Abs(imag(c.left)) <= e[0]) {
		if math.Abs(real(c.right)) <= e[0] && math.Abs(imag(c.right)) <= e[0] {
			return -1
		} else {
			return 1
		}
	}
	return 0
}
func (c CDC) Left() CayleyDickson {
	return CDC{c.left,0}
}
func (c CDC) Right() CayleyDickson {
	return CDC{c.right,0}
}
func (c CDC) Conj() CayleyDickson {
	return CDC{Conj(c.left),-1*c.right}
}
func (c CDC) Neg() CayleyDickson {
	return CDC{-1*c.left,-1*c.right}
}
func (c CDC) SumSqr() float64 {
	r := (real(c.left)*real(c.left)) + (real(c.right)*real(c.right))
	i := (imag(c.left)*imag(c.left)) + (imag(c.right)*imag(c.right))
	return r+i
}
func (c CDC) Length() int {
	return 0
}
func (c CDC) Resize(sz int) CayleyDickson {
	if sz < 0 {
		panic("Resize: invalid size")
	} else if sz == 0 {
		return c
	} else if sz == 1 {
		return CDN{c,CDC{0,0}}
	}
	return CDN{c.Resize(sz-1),NewZero(sz-1)}
}
func (c CDC) Translate(x complex128) CayleyDickson {
	return CDC{c.left+x,c.right+x}
}
func (c CDC) Equal(x CayleyDickson, e ...float64) bool {
	var k float64 = 0
	if len(e) > 0 {
		k = e[0]
	}
	return (x.Length() == 0) &&
		(math.Abs(real(c.left) - real(x.(CDC).left)) <= k) &&
		(math.Abs(imag(c.left) - imag(x.(CDC).left)) <= k) &&
		(math.Abs(real(c.right) - real(x.(CDC).right)) <= k) &&
		(math.Abs(imag(c.right) - imag(x.(CDC).right)) <= k)
}
func (c CDC) String(k ...int) string {
	s := "+"
	sr := "+"
	si := "+"
	if imag(c.left) < 0 {
		sr = ""
	}
	if real(c.right) < 0 {
		s = ""
	}
	if imag(c.right) < 0 {
		si = ""
	}
	q := "k"
	if len(k) > 1 {
		sr = getDir(k[0],k[1])+sr
		s = getDir(k[0]+1,k[1])+s
		si = getDir(k[0]+2,k[1])+si
		q = getDir(k[0]+3,k[1])
	} else {
		s = "i"+s
		si = "j"+si
	}
	return fmt.Sprintf("%g%s%g%s%g%s%g%s", real(c.left),sr,imag(c.left),s,real(c.right),si,imag(c.right),q)
}
func (c CDC) Scale(x float64) CayleyDickson {
	v := complex(x,0)
	return CDC{c.left*v,c.right*v}
}
func (c CDC) Values() []float64 {
	return []float64{real(c.left),imag(c.left),real(c.right),imag(c.right)}
}
func (c CDC) Matrix() [][]float64 {
	if c.Right().IsZero() {
		return [][]float64{[]float64{real(c.left),-1*imag(c.left)},
						   []float64{imag(c.left),real(c.left)}}
	} else {
		return [][]float64{[]float64{real(c.left),-1*imag(c.left),real(c.right),-1*imag(c.right)},
						   []float64{imag(c.left),real(c.left),imag(c.right),real(c.right)},
						   []float64{real(c.right),-1*imag(c.right),real(c.left),-1*imag(c.left)},
						   []float64{imag(c.right),real(c.right),imag(c.left),real(c.left)}}
	}}

func (c CDN) Values() []float64 {
	return append(c.left.Values(),c.right.Values()...)
}
func (c CDN) Matrix() [][]float64 {
	l := c.left.Matrix()
	r := c.right.Matrix()
	k := len(l)
	ret := make([][]float64,k+k)
	for i:=0;i<k;i++ {
		ret[i] = append(l[i],r[i]...)
		ret[i+k] = append(r[i],l[i]...)
	}
	return ret
}
func (c CDN) Scale(x float64) CayleyDickson {
	return CDN{c.left.Scale(x),c.right.Scale(x)}
}
func (c CDN) String(k ...int) string {
	sz := c.Length()
	p := 2 << sz
	s := "+"
	if c.right.Real() < 0 {
		s = ""
	}
	x := 0
	if len(k) > 1 {
		x = k[0]
		sz = k[1]
	}
	y := x+p
	return fmt.Sprintf("%v%s%v",c.left.String(x,sz),s,c.right.String(y,sz))
}
func (c CDN) Equal(x CayleyDickson, e ...float64) bool {
	return (c.Length() == x.Length()) && c.left.Equal(x.Left(),e...) && c.right.Equal(x.Right(),e...)
}
func (c CDN) Length() int {
	return c.left.Length()+1
}
func (c CDN) Resize(sz int) CayleyDickson {
	if sz < 0 {
		panic("Resize: invalid size")
	} else if k := c.Length(); sz == k {
		return c
	} else if sz == k-1 {
		return c.left
	} else if sz < k-1 {
		return c.left.Resize(sz)
	}
	return CDN{c.Resize(sz-1),NewZero(sz-1)}
}
func (c CDN) Translate(x complex128) CayleyDickson {
	return CDN{c.left.Translate(x),c.right.Translate(x)}
}
func (c CDN) Real() float64 {
	return c.left.Real()
}
func (c CDN) Imag() CayleyDickson {
	return CDN{c.left.Imag(),c.right}
}
func (c CDN) Complex128() complex128 {
	return c.left.Complex128()
}
func (c CDN) IsZero(e ...float64) bool {
	return c.left.IsZero(e...) && c.right.IsZero(e...)
}
func (c CDN) IsDirection(e ...float64) bool {
	if c.left.IsDirection(e...) {
		return c.right.IsZero(e...)
	} else if c.left.IsZero(e...) {
		return c.right.IsDirection(e...)
	}
	return false
}
func (c CDN) Direction(e ...float64) int {
	if c.left.IsZero(e...) {
		return c.right.Direction(e...)+(2 << c.left.Length()+1)-1
	}
	return c.left.Direction(e...)
}
func (c CDN) Left() CayleyDickson {
	return c.left
}
func (c CDN) Right() CayleyDickson {
	return c.right
}
func (c CDN) Conj() CayleyDickson {
	return CDN{c.left.Conj(),c.right.Neg()}
}
func (c CDN) Neg() CayleyDickson {
	return CDN{c.left.Neg(),c.right.Neg()}
}
func (c CDN) SumSqr() float64 {
	return c.left.SumSqr()+c.right.SumSqr()
}
