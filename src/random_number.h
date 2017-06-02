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
#endif
