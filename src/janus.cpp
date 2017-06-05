
//
// vdwtest.cpp
//
// simple short-range MD test program
//
// Jun Makino March 9 2015
//
// Known problem as of Mar 13 2015
//   -- acc and phi are consistent only for the case of m=1
// This has been fixed as of Mar 15 2015


#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>

#ifdef DISSIPATIVE_RANDOM
#include "random_number.h"
#endif

#include <cmath>

#define LOWEST_LIMIT 1e-4

// patch informations
#if 1
const int Npatch = 3;
const PS::F64vec3 patch[Npatch] = {PS::F64vec3( 0.0,           1.0, 0.0),
				   PS::F64vec3( 0.86602540378,-0.5, 0.0),
				   PS::F64vec3(-0.86602540378,-0.5, 0.0)};
const PS::F64 coef_r = 396.; // repulsive coefficient
const PS::F64 coef_a[Npatch][Npatch] ={{ 220., 220., 220.},
				       { 220., 220., 220.},
				       { 220., 220., 220.}}; // attractive coefficient of patch
const PS::F64 coef_v = 0.5;  // exponent of f(n,n,r), nu

const PS::F64 tm[Npatch] = {45.0 / 180.0 * M_PI,
			    45.0 / 180.0 * M_PI,
			    45.0 / 180.0 * M_PI}; // theta_m
//const PS::F64 solvent_ratio = 0.5;
const PS::F64 solvent_ratio = 0.0;
#endif

#if 0
const int Npatch = 2;
const PS::F64vec3 patch[Npatch] = {PS::F64vec3( 0.0, 0.0, 1.0),
				   PS::F64vec3( 0.0, 0.0,-1.0)};
const PS::F64 coef_r = 396.;
const PS::F64 coef_a[Npatch][Npatch] = {{88.,88.},
					{88.,88.}};
const PS::F64 coef_v = 0.5;

const PS::F64 tm[Npatch] = {60.0 / 180.0 * M_PI,
			    60.0 / 180.0 * M_PI};

const PS::F64 solvent_ratio = 0.0;
#endif

#if 0
const int Npatch = 1;
const PS::F64vec3 patch[Npatch] = {PS::F64vec3( 0.0, 0.0, 1.0)};

const PS::F64 coef_r = 396.;
const PS::F64 coef_a[Npatch] = {220.};
const PS::F64 coef_v = 0.5;

const PS::F64 tm[Npatch] = {120.0 / 180.0 * M_PI};
const PS::F64 solvent_ratio = 0.0;
#endif

#if (NANOSLIT || NANOTUBE)
const PS::F64 Rwall = 1.0; // width of slit or diameter of tube
const PS::F64 density_wall = 10.0; // density of wall
const PS::F64 coef_a_wall[2] = {100.,100.}; // attractive coefficient of {patchy, solvent}
#endif

// Heat bath parameter
#ifdef NOSE_HOOVER
const PS::F64 Q = 50.0;
#endif

#ifdef DISSIPATIVE_RANDOM
const PS::F64 gamma_dpd = 4.5;
const PS::F64 rn_max = sqrt(3.0);
#endif
//------------------

class Quaternion{
public:
  PS::F64 x,y,z,w;

  Quaternion() : x(1.0),y(0.0),z(0.0),w(0.0){};
  Quaternion(const PS::F64& _s,const PS::F64vec3& _v) : x(_s), y(_v.x),z(_v.y),w(_v.z) {}
  Quaternion(const PS::F64& _x,const PS::F64& _y,const PS::F64& _z, const PS::F64& _w)
    : x(_x), y(_y), z(_z), w(_w) {}
  Quaternion(const Quaternion& src) : x(src.x), y(src.y), z(src.z), w(src.w){}
  Quaternion(const PS::F64& src) : x(src), y(src), z(src), w(src){}

  Quaternion operator+(const Quaternion& rhs) const {
    return Quaternion(x+rhs.x,
		      y+rhs.y,
		      z+rhs.z,
		      w+rhs.w);
  }
  const Quaternion& operator+=(const Quaternion& rhs){
    (*this) = (*this) + rhs;
    return (*this);
  }
  Quaternion operator-(const Quaternion& rhs) const {
    return Quaternion(x-rhs.x,
		      y-rhs.y,
		      z-rhs.z,
		      w-rhs.w);
  }
  const Quaternion& operator-=(const Quaternion& rhs){
    (*this) = (*this) - rhs;
    return (*this);
  }
  Quaternion operator*(const Quaternion& rhs) const {
    const PS::F64 s0 = x;
    const PS::F64 s1 = rhs.x;
    const PS::F64vec3 v0 = PS::F64vec3(y,z,w);
    const PS::F64vec3 v1 = PS::F64vec3(rhs.y,rhs.z,rhs.w);
    return Quaternion(s0 * s1 - v0 * v1, s0 * v1 + v0 * s1 + (v0^v1));
  }
  Quaternion operator*(const PS::F64vec3& _rhs) const { 
    Quaternion rhs(0.0,_rhs);
    return (*this)*rhs;
  }
  friend Quaternion operator*(const PS::F64vec3& _lhs,const Quaternion& rhs){
    Quaternion lhs(0.0,_lhs);
    return lhs*rhs;
  }
  Quaternion operator*(const PS::F64& rhs) const { 
    return Quaternion(x*rhs,
		      y*rhs,
		      z*rhs,
		      w*rhs);
  }
  const Quaternion& operator*=(const PS::F64& rhs){
    (*this) = (*this) * rhs;
    return (*this);
  }
  friend Quaternion operator*(const PS::F64mat4& lhs,const Quaternion& rhs){ 
    return Quaternion(Quaternion(lhs.xx,lhs.xy,lhs.xz,lhs.xw) % rhs,
		      Quaternion(lhs.yx,lhs.yy,lhs.yz,lhs.yw) % rhs,
		      Quaternion(lhs.zx,lhs.zy,lhs.zz,lhs.zw) % rhs,
		      Quaternion(lhs.wx,lhs.wy,lhs.wz,lhs.ww) % rhs);
  }
  PS::F64 operator%(const Quaternion& rhs) const {
    return x*rhs.x + y*rhs.y + z*rhs.z + w*rhs.w;
  }

  const Quaternion& operator=(const Quaternion& rhs){ 
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    w = rhs.w;

    return (*this);
  }
  const Quaternion& operator=(const PS::F64& rhs){
    x = rhs;
    y = rhs;
    z = rhs;
    w = rhs;

    return (*this);
  }

  friend Quaternion normalize(const Quaternion& q){
    PS::F64 ni = 1.0/norm(q);
    return Quaternion(q.x*ni,
		      q.y*ni,
		      q.z*ni,
		      q.w*ni);
  }

  friend Quaternion conj(const Quaternion& q){
    return Quaternion( q.x,
		      -q.y,
		      -q.z,
		      -q.w);
  }
  friend PS::F64 norm(const Quaternion& q){return sqrt(q%q);}

  friend PS::F64mat3asym Rotate(const Quaternion& q){
    return PS::F64mat3asym(q.x*q.x + q.y*q.y - q.z*q.z - q.w*q.w,
			   2.0*(q.y*q.z - q.x*q.w),
			   2.0*(q.y*q.w + q.x*q.z),

			   2.0*(q.y*q.z + q.x*q.w),
			   q.x*q.x-q.y*q.y+q.z*q.z-q.w*q.w,
			   2.0*(q.z*q.w - q.x*q.y),

			   2.0*(q.y*q.w - q.x*q.z),
			   2.0*(q.z*q.w + q.x*q.y),
			   q.x*q.x-q.y*q.y-q.z*q.z+q.w*q.w);
  }
};

template <class T>
T P0(const T a){
  return a;
}
template <class T>
T P1(const T a){
  return T(-a.y, a.x, a.w,-a.z);
}
template <class T>
T P2(const T a){
  return T(-a.z,-a.w, a.x, a.y);
}
template <class T>
T P3(const T a){
  return T(-a.w, a.z,-a.y, a.x);
}

template <class T>
PS::F64mat4 S(const T a){
  return PS::F64mat4(a.x,-a.y,-a.z,-a.w,
		     a.y, a.x, a.w,-a.z,
		     a.z,-a.w, a.x, a.y,
		     a.w, a.z,-a.y, a.x);
}

class BinaryHeader{
public:
  PS::S32 n_ptcl;
  BinaryHeader(const PS::S32 n_ptcl_) : n_ptcl(n_ptcl_){}

  void writeBinary(FILE* fp) const {
    fwrite(&n_ptcl,sizeof(PS::S32),1,fp);
  }
  PS::S32 readBinary(FILE* fp){
    fread(&n_ptcl,sizeof(PS::S32),1,fp);
    return n_ptcl;
  }
};

class CDVHeader{
public:
  PS::F64vec s,e;
  CDVHeader(){}
  CDVHeader(PS::F64vec _s,PS::F64vec _e) : s(_s),e(_e) {}
  PS::S32 readAscii(FILE * fp){
    return 0;
  }
  void writeAscii(FILE* fp) const{
#ifdef NANOSLIT
    fprintf(fp, "'box_sx=%lf,box_sy=%lf,box_sz=%lf,box_ex=%lf,box_ey=%lf,box_ez=%lf\n",
	    s.x,s.y,-0.5*Rwall, e.x,e.y,0.5*Rwall);
#elif defined NANOTUBE
    fprintf(fp, "'box_sx=%lf,box_sy=%lf,box_sz=%lf,box_ex=%lf,box_ey=%lf,box_ez=%lf\n",
	    -0.5*Rwall,-0.5*Rwall,s.z, 0.5*Rwall,0.5*Rwall,e.z);
#else
    fprintf(fp, "'box_sx=%lf,box_sy=%lf,box_sz=%lf,box_ex=%lf,box_ey=%lf,box_ez=%lf\n",
	    s.x,s.y,s.z, e.x,e.y,e.z);
#endif
    fprintf(fp, "'r0=0.35\n");
    fprintf(fp, "'r1=0.2\n");
    fprintf(fp, "'r2=0.1\n");
  }
};

class Force{
public:
  PS::F64vec3 force;
  PS::F64vec3 torque;
  PS::F64 pot;
  void clear(){
    force  = 0.0;
    torque = 0.0;
    pot    = 0.0;
  }
};

class FP{
public:
  PS::S32 id;
  PS::S32 type;
  PS::F64 mass;
  PS::F64vec3 pos;
  PS::F64vec3 vel;
  PS::F64vec3 force;

  Quaternion  angle;
  PS::F64vec3 angvel;
  PS::F64vec3 torque;

  PS::F64 pot;
  PS::F64 search_radius;
  PS::F64 getRsearch() const {
    return this->search_radius;
  }
  PS::F64vec getPos() const { return pos; }
  void copyFromForce(const Force & f){
    force  = f.force;
    torque = f.torque;
    pot    = f.pot;
  }

  // writeXXX must be a const member function
  void writeAscii(FILE* fp) const {
    fprintf(fp, "%d %d %lf %lf %lf %lf %lf %lf %lf\n",
	    (Npatch+1)*id, type, pos.x, pos.y, pos.z, angle.x, angle.y, angle.z, angle.w);
    if(type != 1)
      for(int i=0;i<Npatch;i++){
	const PS::F64vec n = 0.2 * (Rotate(angle)*patch[i]);
	fprintf(fp, "%d %d %lf %lf %lf\n",
		(Npatch+1)*id+i+1, 2, pos.x+n.x, pos.y+n.y, pos.z+n.z);
      }
  }
  void readAscii(FILE* fp){
    fscanf(fp, "%d %d %lf %lf %lf\n",
	   &id, &type, &pos.x, &pos.y, &pos.z);
  }

  void writeBinary(FILE* fp) const {
    fwrite(this,sizeof(FP),1,fp);
  }
  void readBinary(FILE* fp){
    fread(this,sizeof(FP),1,fp);
  }

  Quaternion RichardsonMethod(const Quaternion& angle,const PS::F64vec3& angvel,const PS::F64& dth){
    //const Quaternion angvel1 = S(angle) * Quaternion(0.0,Rotate(angle) * Trans(Rotate(angle)) * angvel) * 0.5;
    const Quaternion angvel1 = S(angle) * Quaternion(0.0, angvel) * 0.5;
    const Quaternion angle1  = normalize(angle + angvel1 * dth * 2.0);
    const Quaternion angle1h = normalize(angle + angvel1 * dth);
    //const Quaternion angvel2 = S(angle1h) * Quaternion(0.0,Rotate(angle1h) * Trans(Rotate(angle1h)) * angvel) * 0.5;
    const Quaternion angvel2 = S(angle1h) * Quaternion(0.0, angvel) * 0.5;
    const Quaternion angle2  = normalize(angle1h + angvel2*dth);

    return normalize(angle2 * 2.0 - angle1);
  }
#ifdef NOSE_HOOVER
  void ApplyHeatBath(const PS::F64& zeta,const PS::F64& dt){
    const PS::F64 factor = exp(-zeta * dt*0.5);
    vel    *= factor;
    angvel *= factor;
  }
#endif
  void IntegrateBeforeForceCalc(const PS::F64 dt,const PS::F64 lh){
    const PS::F64 dth = dt * 0.5;
    vel += force * dth;
    if(type != 1) angvel += torque * dth;

    pos += vel * dt;
#ifndef NANOTUBE
    if(pos.x >   lh) pos.x -= 2.0*lh;
    if(pos.x <= -lh) pos.x += 2.0*lh;
    if(pos.y >   lh) pos.y -= 2.0*lh;
    if(pos.y <= -lh) pos.y += 2.0*lh;
#endif
#ifndef NANOSLIT
    if(pos.z >   lh) pos.z -= 2.0*lh;
    if(pos.z <= -lh) pos.z += 2.0*lh;
#endif
    if(type != 1) angle = RichardsonMethod(angle,angvel,dth);
  }

  void IntegrateAfterForceCalc(const PS::F64 dt){
    const PS::F64 dth = dt * 0.5;
    vel += force * dth;
    if(type != 1) angvel += torque * dth;
  }

  void CalcWallForce(){
#ifdef NANOSLIT
    const PS::F64 Rup = 0.5*Rwall - pos.z;
    const PS::F64 Rdw = 0.5*Rwall + pos.z;

    if(Rup <= 1.0)
      force.z -= 0.16666666666 * M_PI * density_wall * coef_a_wall[type]
	* (1 -2.0 * Rup + 2.0 * Rup*Rup*Rup - Rup*Rup*Rup*Rup) * pos.z / fabs(pos.z);
    if(Rdw <= 1.0)
      force.z -= 0.16666666666 * M_PI * density_wall * coef_a_wall[type]
	* (1 -2.0 * Rdw + 2.0 * Rdw*Rdw*Rdw - Rdw*Rdw*Rdw*Rdw) * pos.z / fabs(pos.z);
#endif

#ifdef NANOTUBE
  const PS::F64 R = 0.5*Rwall - sqrt(pos.x*pos.x + pos.y*pos.y);


  if(R <= 1.0){
    const PS::F64 ftmp = 0.16666666666 * M_PI * density_wall * coef_a_wall[type] * (1 -2.0 * R + 2.0 * R*R*R - R*R*R*R)/R;
    force.x -=  ftmp * pos.x;
    force.y -=  ftmp * pos.y;
  }
#endif
  }
};

class EPI{
public:
  PS::S32 type;
  PS::F64vec pos;
  Quaternion angle;
#ifdef DISSIPATIVE_RANDOM
  PS::F64vec vel;
#endif

  PS::F64vec getPos() const { return pos;}
  void copyFromFP(const FP & fp){
    pos = fp.pos;
#ifdef DISSIPATIVE_RANDOM
    vel = fp.vel;
#endif
    angle = fp.angle;
    type  = fp.type;
  }
};

class EPJ{
public:
  PS::S32 type;
  PS::F64vec3 pos;
  Quaternion angle;
#ifdef DISSIPATIVE_RANDOM
  PS::F64vec vel;
#endif
  PS::F64 search_radius;
  void copyFromFP(const FP & fp){ 
    pos = fp.pos;
#ifdef DISSIPATIVE_RANDOM
    vel = fp.vel;
#endif
    angle = fp.angle;
    type = fp.type;
    search_radius = fp.search_radius;
  }
  PS::F64 getRSearch() const{
    return this->search_radius;
  }
  PS::F64vec getPos() const { return pos; }
  void setPos(const PS::F64vec & pos_new){
    pos = pos_new;
  }
  //PS::F64 getCharge() const { return charge; }
};

//#define APPROXIMATE_SINSIN
struct CalcForceEpEp{
#ifdef DISSIPATIVE_RANDOM
  const PS::F64 dt;
  const PS::F64 temperature;
#endif

  CalcForceEpEp(
#ifdef DISSIPATIVE_RANDOM
		const PS::F64 _dt,
		const PS::F64 _temperature
#endif
		)
#ifdef DISSIPATIVE_RANDOM
    : dt(_dt), temperature(_temperature)
#endif
  {}
  void operator () (const EPI * ep_i,
		    const PS::S32 n_ip,
		    const EPJ * ep_j,
		    const PS::S32 n_jp,
		    Force * force){
    const PS::F64 ph = M_PI*0.5;

#ifdef APPROXIMATE_SINSIN
    const PS::F64 c0 = ph*tmi;
    const PS::F64 c1 = (c0 - c0*c0*c0)/6.0;
    const PS::F64 c2 = c0*(3.0*c0*c0*c0 - 10.0*c0*c0 + 7.0)/360.0;
    const PS::F64 c3 = (-3.0*powf(c0,7)+21*powf(c0,5)-49.0*c0*c0*c0 + 31.0*c0)/15120.0;
    const PS::F64 c4 = c0*(5.0*powf(c0,8)-60.0*powf(c0,6)+294.0*c0*c0*c0*c0 - 620.0*c0*c0 + 381.0)/1814400.0;
#endif

#ifdef DISSIPATIVE_RANDOM
    const PS::F64 sqrtdti = 1.0 / sqrt(dt);
    const PS::F64 sigma_dpd = sqrt(2.0 * temperature * gamma_dpd);
    static XORShift rn;
#endif

    for(int i=0;i<n_ip;i++){
      const PS::F64vec3 ri = ep_i[i].pos;
      const Quaternion  ai = ep_i[i].angle;
      const PS::S64 type_i = ep_i[i].type;
      PS::F64vec3 force_i  = 0.0;
      PS::F64vec3 torque_i = 0.0;
      PS::F64 pot_i = 0.0;

      for(int j=0;j<n_jp;j++){
	const PS::F64vec3 rj = ep_j[j].pos;
	const PS::F64vec3 dr = ri - rj;
	const PS::F64 r2 = dr*dr;
	if(r2 > 1.0 || r2 == 0.0) continue;
	const PS::F64 rinv = 1.0 / sqrt(r2);
	const PS::F64 r2i = rinv*rinv;
	const PS::F64 r   = r2*rinv;
	// repulsive force
	force_i += dr * (coef_r*(1.0 - r) * rinv);
	pot_i += 0.5 * coef_r * (1.0 - r) * (1.0 - r);
#ifdef DISSIPATIVE_RANDOM
	// dissipative force
	const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
	const PS::F64 wij = 1.0 - r2;
	const PS::F64 fd = gamma_dpd * wij*wij * (dv*dr) * rinv;
	force_i -= fd * dr;
	// random force
	const PS::F64 fr = sigma_dpd * wij * rinv * rn.drand() * rn_max * sqrtdti;
	force_i += fr * dr;
#endif
	if(type_i == 1 || ep_j[j].type == 1) continue;
	// attractive force and torque
	const Quaternion  aj = ep_j[j].angle;
	for(int k=0;k<Npatch;k++){
	  const PS::F64vec3 ni = Rotate(ai)*patch[k];
	  const PS::F64 cos_tm_i = cos(tm[k]);
	  for(int l=0;l<Npatch;l++){
	    const PS::F64vec3 nj = Rotate(aj)*patch[l];
	    const PS::F64 cos_ti = -(ni * dr) * rinv;
	    const PS::F64 cos_tj =  (nj * dr) * rinv;
	    const PS::F64 cos_tm_j = cos(tm[l]);

	    if(cos_ti < cos_tm_i || cos_tj < cos_tm_j) continue;
	    const PS::F64 ti = acos(cos_ti);
	    const PS::F64 tj = acos(cos_tj);

	    const PS::F64 cos_phtitm = cos(ph*ti/tm[k]);
	    const PS::F64 cos_phtjtm = cos(ph*tj/tm[l]);

	    const PS::F64 ff = cos_phtitm * cos_phtjtm;
	    const PS::F64 fv  = (ff==0.0) ? 0.0 : powf(ff, coef_v);
	    const PS::F64 fvi = (ff==0.0) ? 0.0 : powf(ff, coef_v - 1.0);

	    const PS::F64vec dcostidr = -(ni - dr*((ni*dr)*r2i))*rinv;
	    const PS::F64vec dcostjdr =  (nj - dr*((nj*dr)*r2i))*rinv;

#ifdef APPROXIMATE_SINSIN
	    PS::F64 dtidcossin;
	    if(ti > LOWEST_LIMIT)
	      dtidcossin = - sin(ph*ti/tm[k]) / sin(ti);
	    else{
	      const PS::F64 ti2 = ti*ti;
	      dtidcossin = -(c0 + ti2 * (c1 + ti2 * (c2 + ti2 * (c3 + ti2*c4))));
	    }

	    PS::F64 dtjdcossin;
	    if(tj > LOWEST_LIMIT)
	      dtjdcossin = - sin(ph*tj/tm[l]) / sin(tj);
	    else{
	      const PS::F64 tj2 = tj*tj;
	      dtjdcossin = -(c0 + tj2 * (c1 + tj2 * (c2 + tj2 * ( c3 + tj2*c4))));
	    }
#else
	    PS::F64 dtidcossin;
	    if(cos_ti*cos_ti != 1.0)
	      dtidcossin = - sin(ph*ti/tm[k]) / sqrt(1.0 - cos_ti*cos_ti);
	    else
	      dtidcossin = - ph/tm[k];

	    PS::F64 dtjdcossin;
	    if(cos_tj*cos_tj != 1.0)
	      dtjdcossin = - sin(ph*tj/tm[l]) / sqrt(1.0 - cos_tj*cos_tj);
	    else
	      dtjdcossin = - ph/tm[l];
#endif
	    const PS::F64 tmpi = dtidcossin * cos_phtjtm;
	    const PS::F64 tmpj = dtjdcossin * cos_phtitm;

	    pot_i   -= 0.5*fv*coef_a[k][l]*(r - r2);
	    force_i += (coef_a[k][l] * fv * (0.5 - r)*rinv)*dr;
	    force_i -= (0.5*coef_a[k][l] * (r - r2) * coef_v * fvi * ph) * (tmpi*dcostidr/tm[k] + tmpj*dcostjdr/tm[l]);
	    const PS::F64 ttmp = 0.5*ph*coef_a[k][l]/tm[k]*(1.0 - r)*coef_v*fvi*tmpi;
	    const PS::F64vec g = ttmp * dr;
	    torque_i += ni ^ g;
	  }
	}
      }
      force[i].pot    = pot_i;
      force[i].force  = force_i;
      force[i].torque = torque_i;
    }
  }
};

void MakeFaceCubicCenter(const PS::S32 n_tot,
			 PS::F64 *&mass,
			 PS::F64vec *&pos,
			 PS::F64vec *&vel,
			 Quaternion *&angle,
			 PS::F64vec *&angvel,
			 const double density,
			 const int seed = 0){
  //static const double PI = atan(1.0) * 4.0;
  PS::MTTS mt;
  double cell_size = pow((double)n_tot/density,1./3.);
  int nunit = 1;
  while(4*nunit*nunit*nunit < n_tot) nunit++;
  if (n_tot != 4*nunit*nunit*nunit){
    std::cerr << "MakeFaceCubicCenter: n_tot and 4*nunit^3 must be the same. "
	      << n_tot << "!= " << 4*nunit*nunit*nunit <<std::endl;
    PS::Abort();
  }
  mt.init_genrand(PS::Comm::getRank()*PS::Comm::getNumberOfThread()+PS::Comm::getThreadNum());

  double unit_size = cell_size/(double)nunit;
  double ush = unit_size * 0.5;
  PS::F64vec unit[4];
  unit[0].x = 0.0; unit[1].x = ush; unit[2].x = 0.0; unit[3].x = ush;
  unit[0].y = 0.0; unit[1].y = ush; unit[2].y = ush; unit[3].y = 0.0;
  unit[0].z = 0.0; unit[1].z = 0.0; unit[2].z = ush; unit[3].z = ush;

  int ip=0;
  for(int i=0; i<nunit; i++){
    for(int j=0; j<nunit; j++){
      for(int k=0; k<nunit; k++){
	for(int l=0; l<4; l++){
	  pos[ip].x = i*unit_size + unit[l].x + 0.1*ush;
	  pos[ip].y = j*unit_size + unit[l].y + 0.1*ush;
	  pos[ip].z = k*unit_size + unit[l].z + 0.1*ush;

	  angle[ip].x = 1.0;
	  angle[ip].y = 0.0;
	  angle[ip].z = 0.0;
	  angle[ip].w = 0.0;
	  ip++;
	}
      }
    }
  }
  assert(ip == n_tot);

  for(int i=0; i<n_tot; i++){
    mass[i] = 1.0;
    const double v_max = 0.1;
    do {
      vel[i][0] = (2. * mt.genrand_res53() - 1.) * v_max;
      vel[i][1] = (2. * mt.genrand_res53() - 1.) * v_max;
      vel[i][2] = (2. * mt.genrand_res53() - 1.) * v_max;
    }while(vel[i] * vel[i] >= v_max * v_max);
  }

  for(int i=0; i<n_tot; i++){
#if 1
    const double v_max = 0.1;
    do {
      angvel[i].x = (2. * mt.genrand_res53() - 1.) * v_max;
      angvel[i].y = (2. * mt.genrand_res53() - 1.) * v_max;
      angvel[i].z = (2. * mt.genrand_res53() - 1.) * v_max;
    }while(angvel[i] * angvel[i] >= v_max * v_max);
#else
    angvel[i].x = 0.0;
    angvel[i].y = 0.0;
    angvel[i].z = 0.0;
#endif
  }


  for(int i=0; i<n_tot; i++){
    pos[i].x -= 0.5*cell_size;
    pos[i].y -= 0.5*cell_size;
    pos[i].z -= 0.5*cell_size;

    if(pos[i].x >   0.5*cell_size) pos[i].x -= cell_size;
    if(pos[i].y >   0.5*cell_size) pos[i].y -= cell_size;
    if(pos[i].z >   0.5*cell_size) pos[i].z -= cell_size;
    if(pos[i].x <= -0.5*cell_size) pos[i].x += cell_size;
    if(pos[i].y <= -0.5*cell_size) pos[i].y += cell_size;
    if(pos[i].z <= -0.5*cell_size) pos[i].z += cell_size;
  }

  PS::F64vec cm_vel = 0.0;
  double  cm_mass = 0.0;
  for(int i=0; i<n_tot; i++){
    cm_vel += mass[i] * vel[i];
    cm_mass += mass[i];
  }
  cm_vel /= cm_mass;
  for(int i=0; i<n_tot; i++){
    vel[i] -= cm_vel;
  }
}


void MakePlane(const PS::S32 n_tot,
	       PS::F64    * mass,
	       PS::F64vec * pos,
	       PS::F64vec * vel,
	       Quaternion * angle,
	       PS::F64vec * angvel,
	       const double density,
	       const int seed = 0){
  //static const double PI = atan(1.0) * 4.0;
  PS::MTTS mt;
  double cell_size = sqrt((double)n_tot/density);
  fprintf(stderr,"%d boxdh = %lf\n",PS::Comm::getRank(),cell_size*0.5);
  int nunit = 1;
  while(nunit*nunit < n_tot) nunit++;
  if (n_tot != nunit*nunit){
    std::cerr << "MakeFaceCubicCenter: n_tot and nunit^2 must be the same. "
	      << n_tot << "!= " << nunit*nunit <<std::endl;
    PS::Abort();
  }
  mt.init_genrand(PS::Comm::getRank()*PS::Comm::getNumberOfThread()+PS::Comm::getThreadNum());

  double unit_size = cell_size/(double)nunit;
  int ip=0;
  for(int i=0; i<nunit; i++){
    for(int j=0; j<nunit; j++){
      pos[ip].x = i*unit_size;
      pos[ip].y = j*unit_size;
      pos[ip].z = 0.0;

      angle[ip].x = 1.0;
      angle[ip].y = 0.0;
      angle[ip].z = 0.0;
      angle[ip].w = 0.0;
      ip++;
    }
  }
  assert(ip == n_tot);

  for(int i=0; i<n_tot; i++){
    mass[i] = 1.0;
    const double v_max = 0.1;
    do {
      vel[i][0] = (2. * mt.genrand_res53() - 1.) * v_max;
      vel[i][1] = (2. * mt.genrand_res53() - 1.) * v_max;
      vel[i][2] = (2. * mt.genrand_res53() - 1.) * v_max;
      //vel[i][2] = 0.0;
    }while(vel[i] * vel[i] >= v_max * v_max);
  }

  for(int i=0; i<n_tot; i++){
    const double v_max = 0.1;
    do {
      angvel[i].x = (2. * mt.genrand_res53() - 1.) * v_max;
      angvel[i].y = (2. * mt.genrand_res53() - 1.) * v_max;
      angvel[i].z = (2. * mt.genrand_res53() - 1.) * v_max;
    }while(angvel[i] * angvel[i] >= v_max * v_max);
  }


  for(int i=0; i<n_tot; i++){
    pos[i].x -= 0.5*cell_size;
    pos[i].y -= 0.5*cell_size;

    if(pos[i].x >   0.5*cell_size) pos[i].x -= cell_size;
    if(pos[i].y >   0.5*cell_size) pos[i].y -= cell_size;
    if(pos[i].x <= -0.5*cell_size) pos[i].x += cell_size;
    if(pos[i].y <= -0.5*cell_size) pos[i].y += cell_size;

    assert(-0.5*cell_size < pos[i].x && pos[i].x <= 0.5*cell_size);
    assert(-0.5*cell_size < pos[i].y && pos[i].y <= 0.5*cell_size);
  }

  PS::F64vec cm_vel = 0.0;
  double  cm_mass = 0.0;
  for(int i=0; i<n_tot; i++){
    cm_vel += mass[i] * vel[i];
    cm_mass += mass[i];
  }
  cm_vel /= cm_mass;
  for(int i=0; i<n_tot; i++){
    vel[i] -= cm_vel;
  }
}

void MakeChain(const PS::S32 n_tot,
	       PS::F64 *&mass,
	       PS::F64vec *&pos,
	       PS::F64vec *&vel,
	       Quaternion *&angle,
	       PS::F64vec *&angvel,
	       const double density,
	       const int seed = 0){
  //static const double PI = atan(1.0) * 4.0;
  PS::MTTS mt;
  double cell_size = (double)n_tot/density;
  mt.init_genrand(PS::Comm::getRank()*PS::Comm::getNumberOfThread()+PS::Comm::getThreadNum());

  double unit_size = cell_size/(double)n_tot;
  int ip=0;
  for(int i=0; i<n_tot; i++){
    pos[ip].x = 0.0;
    pos[ip].y = 0.0;
    pos[ip].z = unit_size*i;

    angle[ip].x = 1.0;
    angle[ip].y = 0.0;
    angle[ip].z = 0.0;
    angle[ip].w = 0.0;
    ip++;
  }
  assert(ip == n_tot);

  for(int i=0; i<n_tot; i++){
    mass[i] = 1.0;
    const double v_max = 0.1;
    do {
      vel[i][0] = (2. * mt.genrand_res53() - 1.) * v_max;
      vel[i][1] = (2. * mt.genrand_res53() - 1.) * v_max;
      vel[i][2] = (2. * mt.genrand_res53() - 1.) * v_max;
    }while(vel[i] * vel[i] >= v_max * v_max);
  }

  for(int i=0; i<n_tot; i++){
    const double v_max = 0.1;
    do {
      angvel[i].x = (2. * mt.genrand_res53() - 1.) * v_max;
      angvel[i].y = (2. * mt.genrand_res53() - 1.) * v_max;
      angvel[i].z = (2. * mt.genrand_res53() - 1.) * v_max;
    }while(angvel[i] * angvel[i] >= v_max * v_max);
  }


  for(int i=0; i<n_tot; i++){
    pos[i].z -= 0.5*cell_size;
    if(pos[i].z >   0.5*cell_size) pos[i].z -= cell_size;
    if(pos[i].z <= -0.5*cell_size) pos[i].z += cell_size;
  }

  PS::F64vec cm_vel = 0.0;
  double  cm_mass = 0.0;
  for(int i=0; i<n_tot; i++){
    cm_vel += mass[i] * vel[i];
    cm_mass += mass[i];
  }
  cm_vel /= cm_mass;
  for(int i=0; i<n_tot; i++){
    vel[i] -= cm_vel;
  }
}


template<class Tpsys>
void SetParticles(Tpsys & psys,
		  const PS::S32 n_tot,
		  const double density,
		  const double temperature){
#if 0 // particles are generated on each rank (need bug fix)
  PS::F64    *mass   = new PS::F64[n_tot];
  PS::F64vec *pos    = new PS::F64vec[n_tot];
  PS::F64vec *vel    = new PS::F64vec[n_tot];
  Quaternion *angle  = new Quaternion[n_tot];
  PS::F64vec *angvel = new PS::F64vec[n_tot];
#ifdef NANOSLIT
  MakePlane(n_tot, mass, pos, vel, angle, angvel, density);
#elif defined NANOTUBE
  MakeChain(n_tot, mass, pos, vel, angle, angvel, density);
#else
  MakeFaceCubicCenter(n_tot, mass, pos, vel, angle, angvel, density);
#endif

  PS::S32 n_proc = PS::Comm::getNumberOfProc();
  PS::S32 rank = PS::Comm::getRank();

  PS::S32 n_loc = n_tot / n_proc;
  PS::S32 i_h = n_loc*rank;
  if(n_loc%n_proc > rank) i_h += rank;
  if(n_loc%n_proc > rank) n_loc++;

  psys.setNumberOfParticleLocal(n_loc);
  for(int i=0; i<n_loc; i++){
    const int id = i + i_h;
    //assert(id < n_tot);
    psys[i].mass   = mass[id];
    psys[i].pos    = pos[id];
    psys[i].vel    = vel[id];
    psys[i].angle  = angle[id];
    psys[i].angvel = angvel[id];
    psys[i].id     = id;
    psys[i].type   = 0;
    psys[i].search_radius = 3.0;
  }
  const PS::S32 n_solvent = (PS::S32)(n_loc * solvent_ratio);
  int n = 0;
  while(n<n_solvent){
    const PS::S32 pivot = rand()%n_loc;
    if(psys[pivot].type == 0){
      psys[pivot].type = 1;
      psys[pivot].angvel = 0.0;
      n++;

      //std::cout << n << std::endl;
    }
  }

  if(mass   != nullptr) delete [] mass;
  if(pos    != nullptr) delete [] pos;
  if(vel    != nullptr) delete [] vel;
  if(angle  != nullptr) delete [] angle;
  if(angvel != nullptr) delete [] angvel;

#else // all particles are generated on rank 0
  if(PS::Comm::getRank() == 0){
    PS::F64    *mass   = new PS::F64[n_tot];
    PS::F64vec *pos    = new PS::F64vec[n_tot];
    PS::F64vec *vel    = new PS::F64vec[n_tot];
    Quaternion *angle  = new Quaternion[n_tot];
    PS::F64vec *angvel = new PS::F64vec[n_tot];
#ifdef NANOSLIT
    MakePlane(n_tot, mass, pos, vel, angle, angvel, density);
#elif defined NANOTUBE
    MakeChain(n_tot, mass, pos, vel, angle, angvel, density);
#else
    MakeFaceCubicCenter(n_tot, mass, pos, vel, angle, angvel, density);
#endif

    psys.setNumberOfParticleLocal(n_tot);
    for(int i=0; i<n_tot; i++){
      psys[i].mass   = mass[i];
      psys[i].pos    = pos[i];
      psys[i].vel    = vel[i];
      psys[i].angle  = angle[i];
      psys[i].angvel = angvel[i];
      psys[i].id     = i;
      psys[i].type   = 0;
      psys[i].search_radius = 3.0;
    }
    const PS::S32 n_solvent = (PS::S32)(n_tot * solvent_ratio);
    int n = 0;
    while(n<n_solvent){
      const PS::S32 pivot = rand()%n_tot;
      if(psys[pivot].type == 0){
	psys[pivot].type = 1;
	psys[pivot].angvel = 0.0;
	n++;

	//std::cout << n << std::endl;
      }
    }

    if(mass   != nullptr) delete [] mass;
    if(pos    != nullptr) delete [] pos;
    if(vel    != nullptr) delete [] vel;
    if(angle  != nullptr) delete [] angle;
    if(angvel != nullptr) delete [] angvel;

  }else{
    psys.setNumberOfParticleLocal(0);
  }

#endif
  ScaleVelocity(psys,temperature);
}

template<class Tpsys>
void RemoveTotalMomentum(Tpsys &system){
  const PS::S32 n_loc = system.getNumberOfParticleLocal();
  PS::F64vec cm_vel_loc = 0.0;
  PS::F64  cm_mass_loc = 0.0;
  for(int i=0; i<n_loc; i++){
    cm_vel_loc += system[i].mass * system[i].vel;
    cm_mass_loc += system[i].mass;
  }
  PS::F64vec cm_vel=0.0;
  PS::F64 cm_mass=0.0;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.x, &cm_vel.x, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.y, &cm_vel.y, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.z, &cm_vel.z, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_mass_loc, &cm_mass, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
#else
  cm_vel = cm_vel_loc;
  cm_mass = cm_mass_loc;
#endif

  cm_vel /= cm_mass;
  for(int i=0; i<n_loc; i++){
    system[i].vel -= cm_vel;
  }
}

template<class Tpsys>
void ScaleVelocity(Tpsys & system,const PS::F64 T){
  const PS::S32 natom_local = system.getNumberOfParticleLocal();
  PS::F64 etra_loc = 0.0;
  PS::F64vec erot_loc = 0.0;
  for(PS::S32 i=0; i<natom_local; i++){
    etra_loc += system[i].mass * system[i].vel * system[i].vel;
    //const PS::F64vec angvel = Rotate(system[i].angle) * Trans(Rotate(system[i].angle)) * system[i].angvel;
    const PS::F64vec angvel = system[i].angvel;
    erot_loc.x += angvel.x * angvel.x;
    erot_loc.y += angvel.y * angvel.y;
    erot_loc.z += angvel.z * angvel.z;
  }
  etra_loc *= 0.5;
  erot_loc *= 0.5;
  PS::S32 natom;
  PS::F64 etra;
  PS::F64vec erot;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&etra_loc, &etra, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&erot_loc.x, &erot.x, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&erot_loc.y, &erot.y, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&erot_loc.z, &erot.z, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&natom_local, &natom, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
#else
  etra = etra_loc;
  erot = erot_loc;
  natom = natom_local;
#endif
  const PS::F64 s_tra = sqrt(1.5*natom*T / etra);
  for(PS::S32 i=0;i<natom_local;i++) system[i].vel *= s_tra;
#if 1
  const PS::F64 s_rotx = sqrt(0.5*natom*T / erot.x);
  const PS::F64 s_roty = sqrt(0.5*natom*T / erot.y);
  const PS::F64 s_rotz = sqrt(0.5*natom*T / erot.z);
  for(PS::S32 i=0;i<natom_local;i++){
    system[i].angvel.x *= s_rotx;
    system[i].angvel.y *= s_roty;
    system[i].angvel.z *= s_rotz;
  }
#endif
  RemoveTotalMomentum(system);
}

template<class Tpsys>
void CalcKineticEnergy(const Tpsys & system,
		       PS::F64 & etra,
		       PS::F64 & erot,
		       PS::S32 & deg_free){
  PS::F64 etra_loc = 0.0;
  PS::F64 erot_loc = 0.0;
  PS::S32 deg_free_loc = 0;
  const PS::S32 nbody = system.getNumberOfParticleLocal();
  for(PS::S32 i=0; i<nbody; i++){
    etra_loc += 0.5 * system[i].mass * system[i].vel * system[i].vel;
    // intertia is unit {1.0,0.0,0.0; 0.0,1.0,0.0; 0.0,0.0,1.0}
    // A(q) I^-1_b A^T(q) = I (unit)
    /*
    const PS::F64mat3asym A = Rotate(system[i].angle) * Trans(Rotate(system[i].angle));
    printf(" %lf %lf %lf\n %lf %lf %lf\n %lf %lf %lf\n\n",
	   A.xx,A.xy,A.xz,A.yx,A.yy,A.yz,A.zx,A.zy,A.zz);
    //*/
    if(system[i].type == 0){
      const PS::F64vec angvel = system[i].angvel;
      erot_loc += 0.5 * angvel * angvel;
      deg_free_loc += 6;
    }else{
      deg_free_loc += 3;
    }
  }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&etra_loc, &etra, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&erot_loc, &erot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&deg_free_loc, &deg_free, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
#else
  etra = etra_loc;
  erot = erot_loc;
  deg_free = deg_free_loc;
#endif
}

template<class Tpsys>
void CalcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & etra,
		PS::F64 & erot,
                PS::F64 & epot,
                const bool clear=true){
  if(clear){
    etot = etra = erot = epot = 0.0;
  }
  PS::F64 etot_loc = 0.0;
  PS::F64 etra_loc = 0.0;
  PS::F64 erot_loc = 0.0;
  PS::F64 epot_loc = 0.0;
  const PS::S32 nbody = system.getNumberOfParticleLocal();
  for(PS::S32 i=0; i<nbody; i++){
    etra_loc += 0.5 * system[i].mass * system[i].vel * system[i].vel;
    // intertia is unit {1.0,0.0,0.0; 0.0,1.0,0.0; 0.0,0.0,1.0}
    // A(q) I^-1_b A^T(q) = I (unit)
    /*
    const PS::F64mat3asym A = Rotate(system[i].angle) * Trans(Rotate(system[i].angle));
    printf(" %lf %lf %lf\n %lf %lf %lf\n %lf %lf %lf\n\n",
	   A.xx,A.xy,A.xz,A.yx,A.yy,A.yz,A.zx,A.zy,A.zz);
    //*/
    const PS::F64vec angvel = system[i].angvel;
    erot_loc += 0.5 * angvel * angvel;

    epot_loc += system[i].pot;
  }
  epot_loc *= 0.5;
  etot_loc = etra_loc + erot_loc + epot_loc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&etot_loc, &etot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&epot_loc, &epot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&etra_loc, &etra, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&erot_loc, &erot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
#else
  etot = etot_loc;
  epot = epot_loc;
  etra = etra_loc;
  erot = erot_loc;
#endif
}

int main(int argc, char *argv[]){
  PS::Initialize(argc, argv);

  #if (NANOSLIT && NANOTUBE)
  fprintf(stderr,"Do not enable both NANOSLIT and NANOTUBE!\n");
  PS::Abort();
  #endif

  std::cout<<std::setprecision(15);
  std::cerr<<std::setprecision(15);
  PS::F64 theta = 0.5;
  const PS::S32 n_leaf_limit = 8;

  long long int n_tot = 256;
  PS::F64 density     = 1.05;
  PS::F64 temperature = 0.8;

  PS::F64 dt = 0.0001;

  PS::S32 n_group_limit = 64;
  PS::S32 nstep     = 1000;
  PS::S32 nstep_eq  = 1000;
  PS::S32 nstep_snp = 100;
  PS::S32 nstep_diag = 100;

  std::string input_file = "";

  //char sinput[1024];
  char dir_name[1024];
  int c;
  sprintf(dir_name,"./result");
  while((c=getopt(argc,argv,"o:N:d:T:s:e:S:D:t:c:n:i:h")) != -1){
    switch(c){
    case 'o':
      sprintf(dir_name,optarg);
      break;
    case 'N':
      n_tot = atoi(optarg);
      std::cerr<<"n_tot="<<n_tot<<std::endl;
      break;
    case 'd':
      density = atof(optarg);
      std::cerr<<"density="<<density<<std::endl;
      break;
    case 'T':
      temperature = atof(optarg);
      std::cerr<<"temperature="<<temperature<<std::endl;
      break;
    case 's':
      nstep = atoi(optarg);
      std::cerr<<"nstep="<<nstep<<std::endl;
      break;
    case 'e':
      nstep_eq = atoi(optarg);
      std::cerr<<"nstep_eq="<<nstep_eq<<std::endl;
      break;
    case 'S':
      nstep_snp = atoi(optarg);
      std::cerr<<"nstep_snp="<<nstep_snp<<std::endl;
      break;
    case 'D':
      nstep_diag = atoi(optarg);
      std::cerr<<"nstep_diag="<<nstep_diag<<std::endl;
      break;
    case 't':
      dt = atof(optarg);
      std::cerr<<"dt="<<dt<<std::endl;
      break;
    case 'n':
      n_group_limit = atoi(optarg);
      std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
      break;
    case 'i':
      input_file = optarg;
      std::cerr<<"input_file="<<input_file<<std::endl;
      break;
    case 'h':
      std::cerr<<"N: n_tot (default: 1000)"<<std::endl;
      std::cerr<<"d: number density (default: 1.05)"<<std::endl;
      std::cerr<<"T: temperature (default: 0.8)"<<std::endl;
      std::cerr<<"s: number of steps (default: 1000)"<<std::endl;
      std::cerr<<"S: time step for snapshot(default: 100)"<<std::endl;
      std::cerr<<"D: time step for diag(default: 100)"<<std::endl;
      std::cerr<<"e: number of steps for equilibration(default: 1000)"<<std::endl;
      std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
      std::cerr<<"i: checkpoint file name"<<std::endl;
      std::cerr<<"t: time step (default: 0.0001)"<<std::endl;
      std::cerr<<"n: n_group_limit (default: 64.0)"<<std::endl;
      return 0;
    }
  }
#ifdef NANOSLIT
  PS::F64 boxdh = 0.5*sqrt((double)n_tot/density);
#elif defined NANOTUBE
  PS::F64 boxdh = 0.5*n_tot/density;
#else
  PS::F64 boxdh = 0.5*powf((double)n_tot/density,1./3.);
#endif
  if(PS::Comm::getRank()==0) fprintf(stderr, "boxdh = %lf\n",boxdh);
  struct stat st;
  if(stat(dir_name, &st) != 0) {
    PS::S32 rank = PS::Comm::getRank();
    PS::S32 ret_loc, ret=0;
    if(rank == 0)
      ret_loc = mkdir(dir_name, 0777);
    PS::Comm::broadcast(&ret_loc, ret);
    if(ret == 0) {
      if(rank == 0)
	fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
      fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
      PS::Abort();
      exit(0);
    }
  }

  std::ofstream fout_eng;
  std::ofstream fout_tcal;
  if(PS::Comm::getRank() == 0){
    char sout_de[1024];
    char sout_tcal[1024];
    sprintf(sout_de, "%s/t-de.dat", dir_name);
    sprintf(sout_tcal, "%s/t-tcal.dat", dir_name);
    std::cerr<<sout_de<<std::endl;
    std::cerr<<sout_tcal<<std::endl;
    fout_eng.open(sout_de);
    fout_tcal.open(sout_tcal);
  }

  PS::ParticleSystem<FP> system_janus;
  system_janus.initialize();

  PS::S32 n_grav_glb = n_tot;
  if(input_file == ""){
    SetParticles(system_janus, n_tot, density, temperature);
    //if(PS::Comm::getRank()==0) fprintf(stderr,"Particles are generated!\n");
  }else{
    system_janus.setNumberOfParticleLocal(n_tot);
    BinaryHeader header(n_tot);
    system_janus.readParticleBinary(input_file.c_str(),header);
    if(PS::Comm::getRank()==0) fprintf(stderr,"Particles are read from %s!\n",input_file.c_str());
  }

  const PS::F64 coef_ema = 0.3;
  PS::DomainInfo dinfo;
  dinfo.initialize(coef_ema);

#ifdef NANOSLIT
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
  dinfo.setPosRootDomain(PS::F64vec(-boxdh,-boxdh,-Rwall),
			 PS::F64vec( boxdh, boxdh, Rwall));
#elif defined NANOTUBE
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Z);
  dinfo.setPosRootDomain(PS::F64vec(-Rwall,-Rwall,-boxdh),
			 PS::F64vec( Rwall, Rwall, boxdh));
#else
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(-boxdh,-boxdh,-boxdh),
			 PS::F64vec( boxdh, boxdh, boxdh));
#endif
  dinfo.collectSampleParticle(system_janus);
  dinfo.decomposeDomain();
  system_janus.exchangeParticle(dinfo);

  //PS::S32 n_grav_loc = system_janus.getNumberOfParticleLocal();
  PS::TreeForForceShort<Force, EPI, EPJ>::Scatter tree_janus;
  tree_janus.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
  tree_janus.calcForceAllAndWriteBack(CalcForceEpEp(
#ifdef DISSIPATIVE_RANDOM
						    dt,temperature
#endif
						    ),  system_janus, dinfo);

  PS::F64 Epot0, Etra0, Erot0, Etot0, Epot1, Etra1, Erot1, Etot1;
  ScaleVelocity(system_janus,temperature);
  CalcEnergy(system_janus, Etot0, Etra0, Erot0, Epot0);
  //if(PS::Comm::getRank() == 0) printf("Etot = %lf, Epot = %lf, Etra = %lf, Erot = %lf\n",Etot0,Epot0,Etra0,Erot0);

  PS::S32 snp_id = 0;
  PS::S32 time_snp = 0;
  PS::S32 time_diag = 0;
  bool isInitialized = false;
  PS::F64 time_sys = -dt * nstep_eq;

  PS::F64 Epot_ave = 0.0, Ekin_ave = 0.0;
  int n_loc = system_janus.getNumberOfParticleLocal();
#ifdef NOSE_HOOVER
  PS::F64 zeta = 0.0, log_s = 0.0;
#endif
  for(int s=-nstep_eq;s<nstep;s++){
    if(s < 0){
      ScaleVelocity(system_janus,temperature);
#ifdef NOSE_HOOVER
      zeta = 0.0; log_s = 0.0;
#endif
    }
    if(s%1000==0) RemoveTotalMomentum(system_janus);
    PS::Timer timer;
    timer.reset();
    timer.start();
    if(s == time_snp){
      CDVHeader cdvh(PS::F64vec(-boxdh,-boxdh,-boxdh),PS::F64vec(boxdh,boxdh,boxdh));
      char filename[256];
      sprintf(filename, "%s/%05d.cdv", dir_name, snp_id);
      system_janus.writeParticleAscii(filename, cdvh);

      BinaryHeader bh(n_tot);
      sprintf(filename,"%s/checkpoint", dir_name);
      system_janus.writeParticleBinary(filename,bh);

      snp_id++;
      time_snp += nstep_snp;
    }
    if(!isInitialized && s == 0){
      CalcEnergy(system_janus, Etot0, Etra0, Erot0, Epot0);
      //if(PS::Comm::getRank() == 0) fprintf(stderr,"Etot0 = %lf, Epot0 = %lf, Etra0 = %lf, Erot0 = %lf\n",Etot0,Epot0,Etra0,Erot0);
      isInitialized = true;
    }
    for(int i=0;i<n_loc;i++){
#ifdef NOSE_HOOVER
      system_janus[i].ApplyHeatBath(zeta,dt);
#endif
      system_janus[i].IntegrateBeforeForceCalc(dt,boxdh);
    }
    if(s%10 == 0){
      dinfo.collectSampleParticle(system_janus);
      dinfo.decomposeDomain();
    }

    system_janus.exchangeParticle(dinfo);
    n_loc = system_janus.getNumberOfParticleLocal();

    tree_janus.calcForceAllAndWriteBack
      (CalcForceEpEp(
#ifdef DISSIPATIVE_RANDOM
		     dt,temperature
#endif
		     ),  system_janus, dinfo);
#if (NANOSLIT || NANOTUBE)
    for(int i=0;i<n_loc;i++){
      system_janus[i].CalcWallForce();
    }
#endif
#ifdef NOSE_HOOVER
    PS::F64 etra = 0.0, erot = 0.0;
    PS::S32 g = 0;
    CalcKineticEnergy(system_janus,etra,erot,g);

    zeta += (2.0 * (etra + erot) - g*temperature)/Q * dt;
    log_s += zeta*dt;
#endif   
    for(int i=0;i<n_loc;i++){
      system_janus[i].IntegrateAfterForceCalc(dt);
#ifdef NOSE_HOOVER
      system_janus[i].ApplyHeatBath(zeta,dt);
#endif
    }

    CalcEnergy(system_janus, Etot1, Etra1, Erot1, Epot1);
    if(s>=0){
      Epot_ave += Epot1;
      Ekin_ave += Etra1 + Erot1;
    }
    if(s == time_diag) {
      if(PS::Comm::getRank() == 0){
#ifdef NOSE_HOOVER
	const PS::F64 Etstat = 0.5*zeta*zeta*Q + g*temperature*log_s;
	Etot1 += Etstat;

	fout_eng<<time_sys<<"   "<< " " << Epot1 << " " << Etra1 << " " << Erot1 << " " << Etstat << " " << (Etot1-Etot0)/Etot0<<std::endl;
	/*
	fprintf(stderr, "%10.7f %lf %lf %lf %lf %+e\n",
		time_sys, Epot1, Etra1, Erot1, Etstat, (Etot1 - Etot0) / Etot0);
	//*/
#else
	fout_eng<<time_sys<<"   "<< " " << Epot1 << " " << Etra1 << " " << Erot1 << " " <<(Etot1-Etot0)/Etot0<<std::endl;
	/*
	fprintf(stderr, "%10.7f %lf %lf %lf %+e\n",
		time_sys, Epot1, Etra1, Erot1, (Etot1 - Etot0) / Etot0);
	//*/
#endif
	time_diag += nstep_diag;
      }
    }
    time_sys += dt;
  }

  PS::Finalize();
  return 0;
}
