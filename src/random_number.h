#ifndef H_RANDOM_NUMBER
#include <random>
#include <climits>
#include <cfloat>
class XORShift{
 public:
  static constexpr unsigned long min(){return 0uL;}
  static constexpr unsigned long max(){return ULONG_MAX;}

  XORShift(){
    std::random_device rd;
    w = rd();
    printf("%lu\n",w);
  }
  XORShift(const unsigned long s){ w = s; }

  unsigned long rand(){
    unsigned long t = (x^(x<<11));
    x=y; y=z; z=w;
    return w = (w^(w>>19))^(t^(t>>8)) ;
  }
  double drand(){
    union{
      double d;
      unsigned long l;
    }rn;
    rn.l = random();
    return rn.d / DBL_MAX;
  }
 private:
  unsigned long x = 123456789u, y = 362436069u, z = 521288629u, w;
};

#define _LN_2 6.9314718055994528623E-1f
#define _2_TO_MINUS_31 4.6566128730773925781E-10f
#define _2_TO_MINUS_32 2.3283064365386962891E-10f

#define _TEA_K0 0xA341316C
#define _TEA_K1 0xC8013EA4
#define _TEA_K2 0xAD90777D
#define _TEA_K3 0x7E95761E
#define _TEA_DT 0x9E3779B9

class TEA{
 private:
  float bound(const float& x,const float& min,const float& max) const {
    assert(min<max);
    //if(x > max) printf("generated x(=%f) > max\n",x);
    //if(x < min) printf("generated x(=%f) < min\n",x);
    return (x > max) ? max : ((min > x) ? min : x);
  }
  template<int N>
  void __TEA_core(unsigned int& v0, unsigned int& v1, unsigned int sum=0){
    sum += _TEA_DT;
    v0 += ( ( v1 << 4 ) + _TEA_K0 ) ^ ( v1 + sum ) ^ ( ( v1 >> 5 ) + _TEA_K1 );
    v1 += ( ( v0 << 4 ) + _TEA_K2 ) ^ ( v0 + sum ) ^ ( ( v0 >> 5 ) + _TEA_K3 );
    __TEA_core < N - 1 > ( v0, v1, sum );
  }


  template<int N>
  float gaussian_TEA( bool pred, int u, int v ){
    unsigned int v0 =  pred ? u : v;
    unsigned int v1 = !pred ? u : v;
    __TEA_core<N>(v0,v1);
    float f = cos( 2.f * M_PI * ( v0 & 0X7FFFFFFF ) * _2_TO_MINUS_31) * ( ( v0 & 0X80000000 ) ? 1.0f : -1.0f );
    float r = sqrtf( -2.0f * _LN_2 * log2f(v1 * _2_TO_MINUS_32) );
    return bound( r * f, -4.0f, 4.0f );
  }
  unsigned int __float_as_uint( float r ){
    union{
      float f;
      unsigned int u;
    } ret;
    ret.f = r;
    return ret.u;
  }
 public:
  float drand(const float& v0,const float& v1){
    unsigned int u0 = __float_as_uint(v0);
    unsigned int u1 = __float_as_uint(v1);
    return gaussian_TEA<4>(u0>u1,u0,u1);
  }
  float drand(const int& v0,const int& v1){
    return gaussian_TEA<4>(v0>v1,v0,v1);
  }
};

template<>
void TEA::__TEA_core<0>( unsigned int &v0, unsigned int &v1, unsigned int sum) {}

#endif
