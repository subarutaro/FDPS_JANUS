
#pragma once

#include<iostream>
#include<iomanip>
#include<cmath>
#include"vector4.hpp"

namespace ParticleSimulator{
    template<class T>
    class Matrix4{
    public:
      //constructor
      T xx, xy, xz, xw;
      T yx, yy, yz, yw;
      T zx, zy, zz, zw;
      T wx, wy, wz, ww;
      Matrix4()
	: xx(T(0)), xy(T(0)), xz(T(0)), xw(T(0)),
	  yx(T(0)), yy(T(0)), yz(T(0)), yw(T(0)),
	  zx(T(0)), zy(T(0)), zz(T(0)), zw(T(0)),
	  wx(T(0)), wy(T(0)), wz(T(0)), ww(T(0)) {}
      Matrix4(const T _xx, const T _xy, const T _xz, const T _xw,
	      const T _yx, const T _yy, const T _yz, const T _yw,
	      const T _zx, const T _zy, const T _zz, const T _zw,
	      const T _wx, const T _wy, const T _wz, const T _ww)
	: xx(_xx), xy(_xy), xz(_xz), xw(_xw),
	  yx(_yx), yy(_yy), yz(_yz), yw(_yw),
	  zx(_zx), zy(_zy), zz(_zz), zw(_zw),
	  wx(_wx), wy(_wy), wz(_wz), ww(_ww) {}
      Matrix4(const T s) : xx(s), xy(s), xz(s), xw(s),
			   yx(s), yy(s), yz(s), yw(s),
			   zx(s), zy(s), zz(s), zw(s),
			   wx(s), wy(s), wz(s), ww(s) {}
      Matrix4(const Matrix4 & src)
	: xx(src.xx), xy(src.xy), xz(src.xz), xw(src.xw),
	  yx(src.yx), yy(src.yy), yz(src.yz), yw(src.yw),
	  zx(src.zx), zy(src.zy), zz(src.zz), zw(src.zw),
	  wx(src.wx), wy(src.wy), wz(src.wz), ww(src.ww) {}

      const Matrix4 & operator = (const Matrix4 & rhs) {
	xx = rhs.xx; xy = rhs.xy; xz = rhs.xz; xw = rhs.xw;
	yx = rhs.yx; yy = rhs.yy; yz = rhs.yz; yw = rhs.yw;
	zx = rhs.zx; zy = rhs.zy; zz = rhs.zz; zw = rhs.zw;
	wx = rhs.wx; wy = rhs.wy; wz = rhs.wz; ww = rhs.ww;
	return (*this);
      }
      const Matrix4 & operator = (const T s) {
	xx = yx = zx = wx = s;
	xy = yy = zy = wy = s;
	xz = yz = zz = wz = s;
	xw = yw = zw = ww = s;
	return (*this);
      }

      Matrix4 operator + (const Matrix4 & rhs) const {
	return Matrix4(xx + rhs.xx, xy + rhs.xy, xz + rhs.xz, xw + rhs.xw,
		       yx + rhs.yx, yy + rhs.yy, yz + rhs.yz, yw + rhs.yw,
		       zx + rhs.zx, zy + rhs.zy, zz + rhs.zz, zw + rhs.zw,
		       wx + rhs.wx, wy + rhs.wy, wz + rhs.wz, ww + rhs.ww);
      }
      /*
	Matrix4 operator + (const T & rhs) const {
	  return Matrix4(xx + rhs, xy + rhs, xz + rhs, xw + rhs,
			 yx + rhs, yy + rhs, yz + rhs, yw + rhs,
			 zx + rhs, zy + rhs, zz + rhs, zw + rhs,
			 wx + rhs, wy + rhs, wz + rhs, ww + rhs);
        }
      Matrix4 operator - (const T & rhs) const {
	return Matrix4(xx - rhs, yy - rhs, xy, yx);
      }
      const Matrix4 & operator += (const Matrix4 & rhs) {
	(*this) = (*this) + rhs;
	return (*this);
      }
        Matrix4 operator - (const Matrix4 & rhs) const {
            return Matrix4(xx - rhs.xx, yy - rhs.yy, xy - rhs.xy, yx - rhs.yx);
        }
        const Matrix4 & operator -= (const Matrix4 & rhs) {
            (*this) = (*this) - rhs;
            return (*this);
        }
      //*/
      Matrix4 operator * (const T & rhs) const {
	return Matrix4(xx * rhs, yy * rhs, xy * rhs, yx * rhs);
      }
      Vector4<T> operator * (const Vector4<T> & rhs) const {
	return Vector4<T>(Vector4<T>(xx,xy,xz,xw)*rhs,
			  Vector4<T>(yx,yy,yz,yw)*rhs,
			  Vector4<T>(zx,zy,zz,zw)*rhs,
			  Vector4<T>(wx,wy,wz,ww)*rhs);
      }
      Matrix4 operator * (const Matrix4 & rhs) const {
	const Vector4<T> cx(xx,xy,xz,xw);
	const Vector4<T> cy(yx,yy,yz,yw);
	const Vector4<T> cz(zx,zy,zz,zw);
	const Vector4<T> cw(wx,wy,wz,ww);
	const Vector4<T> rx(rhs.xx,rhs.yx,rhs.zx,rhs.wx);
	const Vector4<T> ry(rhs.xy,rhs.yy,rhs.zy,rhs.wy);
	const Vector4<T> rz(rhs.xz,rhs.yz,rhs.zz,rhs.wz);
	const Vector4<T> rw(rhs.xw,rhs.yw,rhs.zw,rhs.ww);
	return Matrix4(cx*rx,cx*ry,cx*rz,cx*rw,
		       cy*rx,cy*ry,cy*rz,cy*rw,
		       cz*rx,cz*ry,cz*rz,cz*rw,
		       cw*rx,cw*ry,cw*rz,cw*rw);
      }
      /*
      Matrix4 operator / (const T & rhs) const {
	return Matrix4(xx / rhs, yy / rhs, xy / rhs, yx / rhs);
      }
      const Matrix4 & operator *= (const T & rhs) {
      (*this) = (*this) * rhs;
      return (*this);
      }
      friend Matrix4 operator * (const T s, const Matrix4 & m) {
      return (m * s);
      }
      const Matrix4 & operator /= (const T & rhs) {
      (*this).xx /= rhs;
      (*this).yy /= rhs;
      (*this).xy /= rhs;
      (*this).yx /= rhs;
      return (*this);
      }

      const Matrix4 & operator + () const {
      return (* this);
      }

      const Matrix4 operator - () const {
      return Matrix4(-xx, -yy, -xy, -yx);
      }


      T getTrace() const {
      return (xx + yy);
      }
      T getDeterminant() const {
      return (xx * yy - xy * yx);
      }
      T getSecondInvariant() const {
      return 0.5 * (- this->getTrace() * this->getTrace() + (*this * *this).getTrace());
      }
      //*/
      Matrix4 getTrans() const {
	return Matrix4(xx, yx, zx, wx,
		       xy, yy, zy, wy,
		       xz, yz, zz, wz,
		       xw, yw, zw, ww);
      }
      /*
      template <typename U>
      operator Matrix4<U> () const {
      return Matrix4<U>( static_cast<U>(xx), static_cast<U>(yy), static_cast<U>(xy), static_cast<U>(yx) );
      }

      friend std::ostream& operator << (std::ostream& c, const Matrix4<T> & mat){
      c<<mat.xx<<"   "<<mat.xy<<std::endl;
      c<<mat.yx<<"   "<<mat.yy<<std::endl;
      return c;
      }
      //*/
    };
}
