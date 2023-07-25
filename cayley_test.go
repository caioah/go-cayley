package cayley
import (
	"testing"
	"caioah/mtx"
	"fmt"
)
var x CayleyDickson = NewPosition([]float64{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15})
var y CayleyDickson = NewPosition([]float64{16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31})
var z CayleyDickson = CDN{x,y}
var sz int = 3
type Stringer interface {
	String(...int) string
}
type Value struct {
	val interface{}
}
func (v Value) String(_ ...int) string {
	return fmt.Sprint(v.val)
}
func testCnd(t *testing.T,name,cmp string,op func(interface{},interface{}) bool,x,y Stringer) {
	if !op(x,y) {
		t.Errorf("Error in %v: %v %v %v", name, x.String(), cmp, y.String())
	} else {
		t.Logf("%v passed",name)
	}
}
func testEq(t *testing.T,name string,x,y CayleyDickson, e ...float64) {
	testCnd(t,name,"!=",func(x,y interface{}) bool { return x.(CayleyDickson).Equal(y.(CayleyDickson),e...) },x,y)
}
func NotEqual(x,y CayleyDickson, e ...float64) bool {
	return !x.Equal(y,e...)
}
func TestCDC (t *testing.T) {
	v := z.Resize(0)
	testEq(t,"Resize",v,CDC{0+1i,2+3i})
	testCnd(t,"Real","!=",func(x,y interface{}) bool { return x == y },Value{v.Real()},Value{x.Real()})
	testCnd(t,"Complex128","!=",func(x,y interface{}) bool { return x == y },Value{v.Complex128()}, Value{x.Complex128()})
	testEq(t,"Left",z.Left(),x)
	testEq(t,"Right",z.Right(),y)
	testEq(t,"Conj",Mul(z.Conj(),z),CDC{complex(z.SumSqr(),0),0}.Resize(sz),1e-16)
	testEq(t,"InvMul",Mul(z,InvMul(z)),CDC{1,0}.Resize(sz),1e-16)
	testEq(t,"InvMul2",Mul(InvMul(z),z),CDC{1,0}.Resize(sz),1e-16)
	w := Add(x,y)
	testEq(t,"Add",x,Add(w,y.Neg()))
	testEq(t,"Add2",y,Add(w,x.Neg()))
	w = Mul(x,y)
	testEq(t,"Mul",x,Mul(w,InvMul(y)),1e-14)
	testCnd(t,"Mul","==",func(x,y interface{}) bool { return !x.(CayleyDickson).Equal(y.(CayleyDickson),5e-1) },y,Mul(w,InvMul(x)))
	w = Mul(y,x)
	testEq(t,"Mul2",y,Mul(w,InvMul(x)),1e-14)
	testCnd(t,"Mul2","==",func(x,y interface{}) bool { return !x.(CayleyDickson).Equal(y.(CayleyDickson),5e-1) },x,Mul(w,InvMul(y)))
	w = Sqrt(z)
	testEq(t,"Sqrt",Mul(w,w),z,1e-14)
	v = Versor(z)
	w = Exp(v)
	testEq(t,"Log",Log(w),v,1e-14)
	w = Log(z)
	testEq(t,"Exp",Exp(w),z,1e-13)
	w = Pow(z,CDC{2,0})
	testEq(t,"Pow",w,Mul(z,z),1e-11)
	t.Logf("%v",z.Resize(0).String())
	t.Logf("%v",z.Resize(1).String())
	t.Logf("%v",z.Resize(2).String())
	t.Logf("%v",z.Resize(3).String())
	t.Logf("%v",z.String())
	w = Add(z,CDN{NewZero(sz),z.Translate(complex(float64(int(2)<<(sz+1)),float64(int(2)<<(sz+1))))})
	t.Logf("%v",w.String())
	t.Logf("%v",z.Values())
	t.Logf("%v",z.Matrix())
	w = NewPosition(mtx.VMul(z.Matrix(),CDC{2,0}.Resize(sz).Values()))
	t.Logf("%v",w.String())
	t.Logf("%v",NewScalar(32,1).String())
	t.Logf("%v",NewScalar(32,1).Direction())
}
func BenchmarkMul(b *testing.B) {
	//b.ResetTimer()
	for i:=0;i<b.N;i++ {
		Mul(CDC{1+2i,2+3i},CDC{3+4i,4+5i})
	}
}
