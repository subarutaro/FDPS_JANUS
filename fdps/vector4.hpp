#pragma once

#include<iostream>
#include<iomanip>

namespace ParticleSimulator{
    template<class T>
    class Vector4{
    public:
      T x, y, z, w;
      Vector4() : x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}
      Vector4(const T _x, const T _y, const T _z, const T _w) : x(_x), y(_y), z(_z), w(_w){}
      Vector4(const T s) : x(s), y(s), z(s), w(s){}
      Vector4(const Vector4 & src) : x(src.x), y(src.y), z(src.z), w(src.w){}

      typedef T DataType;
      static const int DIM = 4;

      const Vector4 & operator = (const Vector4 & rhs){
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
	w = rhs.w;
	return (*this);
      }

      const Vector4 & operator = (const T s){
	x = y = z = w = s;
	return (*this);
      }

      Vector4 operator + (const Vector4 & rhs) const{
	return Vector4(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w);
      }
      const Vector4 & operator += (const Vector4 & rhs){
	(*this) = (*this) + rhs;
	return (*this);
      }
      Vector4 operator - (const Vector4 & rhs) const{
	return Vector4(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
      }

      const Vector4 & operator -= (const Vector4 & rhs){
	(*this) = (*this) - rhs;
	return (*this);
      }

      Vector4 operator * (const T s) const{
	return Vector4(x * s, y * s, z * s, w * s);
      }
      const Vector4 & operator *= (const T s){
	(*this) = (*this) * s;
	return (*this);
      }
      friend Vector4 operator * (const T s, const Vector4 & v){
	return (v * s);
      }
      Vector4 operator / (const T s) const{
	return Vector4(x / s, y / s, z / s, w / s);
      }
      const Vector4 & operator /= (const T s){
	(*this) = (*this) / s;
	return (*this);
      }

      const Vector4 & operator + () const {
	return (* this);
      }

      const Vector4 operator - () const {
	return Vector4(-x, -y, -z, -w);
      }

      T operator * (const Vector4 & rhs) const{
	return (x * rhs.x) + (y * rhs.y) + (z * rhs.z) + (w * rhs.w);
      }
      /*
      Vector4 operator ^ (const Vector4 & rhs) const{
	return Vector4( (y * rhs.z - z * rhs.y),
			(z * rhs.x - x * rhs.z),
			(x * rhs.y - y * rhs.x) );
      }
      //*/
      template <typename U>
      operator Vector4<U> () const {
	return Vector4<U> (static_cast<U>(x),
			   static_cast<U>(y),
			   static_cast<U>(z),
			   static_cast<U>(w));
      }

      T getMax() const {
	T max_val0 = (x > y) ? x : y;
	T max_val1 = (z > w) ? z : w;
	max_val0 = (max_val0 > max_val1 ) ? max_val0 : max_val1;
	return max_val0;
      }

      T getMin() const {
	T min_val0 = (x < y) ? x : y;
	T min_val1 = (z < w) ? z : w;
	min_val0 = (min_val0 < min_val1 ) ? min_val0 : min_val1;
	return min_val0;
      }

      template <class F>
      Vector4 applyEach(F f) const {
	return Vector4(f(x), f(y), f(z), f(w));
      }

      template <class F>
      friend Vector4 ApplyEach(F f, const Vector4 & arg1, const Vector4 & arg2){
	return Vector4( f(arg1.x, arg2.x), f(arg1.y, arg2.y), f(arg1.z, arg2.z), f(arg1.w,arg2.w));
      }

      friend std::ostream & operator <<(std::ostream & c, const Vector4 & u){
	c<<u.x<<"   "<<u.y<<"    "<<u.z<<"   "<<u.w;
	return c;
      }

      friend std::istream & operator >>(std::istream & c, Vector4 & u){
	c>>u.x; c>>u.y; c>>u.z; c>>u.w;
	return c;
      }

#if 0
      const T & operator[](const int i) const {
#ifdef PARTICLE_SIMULATOR_VECTOR_RANGE_CHECK
	if(i >= DIM || i < 0){
	  std::cout<<"PS_ERROR: Vector invalid access. \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;		
	  std::cerr<<"Vector element="<<i<<" is not valid."<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	  MPI::COMM_WORLD.Abort(-1);
#else
	  exit(-1);
#endif		
	}
#endif
	return (&x)[i];
      }

      T & operator[](const int i){
#ifdef PARTICLE_SIMULATOR_VECTOR_RANGE_CHECK
	if(i >= DIM || i < 0){
	  std::cout<<"PS_ERROR: Vector invalid access. \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;		
	  std::cerr<<"Vector element="<<i<<" is not valid."<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	  MPI::COMM_WORLD.Abort(-1);
#else
	  exit(-1);
#endif		
	}
#endif
	return (&x)[i];
      }
#else
      const T & operator[](const int i) const {
	//std::cerr<<"operator []"<<std::endl;
	if(0==i) return x;
	if(1==i) return y;
	if(2==i) return z;
	if(3==i) return w;
	std::cout<<"PS_ERROR: Vector invalid access. \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;		
	std::cerr<<"Vector element="<<i<<" is not valid."<<std::endl;
#  ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	MPI::COMM_WORLD.Abort(-1);
#  endif		
	exit(-1);
      }
      T & operator[](const int i){
	//std::cerr<<"operator []"<<std::endl;
	if(0==i) return x;
	if(1==i) return y;
	if(2==i) return z;
	if(3==i) return w;
	std::cout<<"PS_ERROR: Vector invalid access. \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;		
	std::cerr<<"Vector element="<<i<<" is not valid."<<std::endl;
#  ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	MPI::COMM_WORLD.Abort(-1);
#  endif		
	exit(-1);
      }
#endif
      /*
      T getDistanceSQ(const Vector4 & u) const {
	T dx = x - u.x;
	T dy = y - u.y;
	T dz = z - u.z;
	T dw = w - u.w;
	return dx*dx + dy*dy + dz*dz + dw*dw;
      }
      //*/
      bool operator == (const Vector4 & u) const {
	return ( (x==u.x) && (y==u.y) && (z==u.z) && (w==u.w));
      }
      bool operator != (const Vector4 & u) const {
	return ( (x!=u.x) || (y!=u.y) || (z!=u.z) || (w!=u.w));
      }
	
      /*
	Vector4 getDiagonal (const Vector4 & u) const {
	return Vector4(x*u.x, y*u.y, z*u.z);
	}
	Vector4 getDiagonal (const Vector4<int> & u) const {
	return Vector4(x*(T)(u.x), y*(T)(u.y), z*(T)(u.z));
	}
      */
    };

  template <>
  inline Vector4<float> Vector4<float>::operator / (const float s) const {
    const float inv_s = 1.0f/s;
    return Vector4(x * inv_s, y * inv_s, z * inv_s, w * inv_s);
  }
  template <>
  inline Vector4<double> Vector4<double>::operator / (const double s) const {
    const double inv_s = 1.0/s;
    return Vector4(x * inv_s, y * inv_s, z * inv_s, w * inv_s);
  }
}
