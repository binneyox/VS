/** \file    dust.h
    \brief   Defines a dust model so extinction can be computed down a
    los; also defines obs::los
    \author  James Binney
    \date    2021
*/

#pragma once
#include "coord.h"
#include "obs.h"
#include "units.h"
#include "potential_factory.h"

				   
namespace dust {
/* Class introducing a global spiral pattarn. Distances in internal
 * units
 */
class EXP Spiral{
	double Ampl,phase,alpha,Narms,kz;
	public:
		Spiral(){}
		Spiral(double _Ampl,double _phase,double _alpha,int _Narms,double _kz) :
		    Ampl(_Ampl), phase(_phase), alpha(_alpha), Narms(_Narms), kz(_kz){}			
		double dens(coord::PosCyl Rzphi){
			return Ampl*cos(alpha*log(Rzphi.R+1.e-10)-Narms*(Rzphi.phi-phase))*exp(-kz*Rzphi.z);
		}
};
/* Class creating spiral-shaped dust clouds. (Rc,phic) Galactocentric
 * coords of cloud centre, norm central density
*/
class EXP JTcloud{
	math::LinearInterpolator2d vals;
	double Rc,phic,z0,norm;
	public:
		JTcloud(){}
		JTcloud(double _Rc, double _phic, double _z0, double _norm,const std::string& fname);
		double dens(const coord::PosCyl& Rzphi) const{
			double phi=Rzphi.phi;
			if(fabs(phic)>0.5*M_PI && phi<0) phi+=2*M_PI;
			double x=(Rzphi.R-Rc), y=Rc*(phi-phic);//assume clockwise rotation
			if(fabs(x)>vals.xmax() || fabs(y)>vals.ymax()) return 0;
			return vals.value(x,y)*.5/z0*exp(-fabs(Rzphi.z)/z0);
		}
};
/* Class creating blobs at p
*/
class EXP Blob{
	coord::PosCyl p;
	double Ampl,rad;
	public:
		Blob(){}
		Blob(const coord::PosCyl _p, const double _Ampl, const double _rad) :
		    p(_p), Ampl(_Ampl), rad(_rad){//ensure continuity around 180
			if(fabs(p.phi)>0.5*M_PI && p.phi<0) p.phi+=2*M_PI;
		}
		double dens(const coord::PosCyl& Rzphi) const{
			double phi=Rzphi.phi;
			if(fabs(p.phi)>0.5*M_PI && phi<0) phi+=2*M_PI;
			double x=(Rzphi.R-p.R), y=p.R*(phi-p.phi);//assume clockwise rotation
			//if(fabs(x)>3*rad || fabs(y)>3*rad) return 0;
			return Ampl * exp(-.5*(pow_2(x/rad)+pow_2(y/rad)+pow_2((Rzphi.z-p.z)/rad)));
		}
};

/*  The basic dust model comes with a cavity around the Sun of radius
 *  hole. To this model one can add any number of log spirals or
 *  Julia-Toomre clouds.Input distances in kpc
 *  */
class EXP dustModel{
	private:
		double Rd, zd, log_rho0, Rw, Hw, from_Kpc;
		obs::solarShifter* shifter;
		potential::PtrDensity ptr;
		std::vector<JTcloud>* cl;
		std::vector<Spiral>* sp;
		std::vector<Blob>* bl;
		double Zw(const coord::PosCyl& p) const{
			if(p.R<Rw) return 0;
			else return Hw*(p.R-Rw)/Rw*sin(p.phi);
		}
	public:
		dustModel(const double _Rd,const double _zd,const double _Rw,
			  const double _Hw,const double dAvds,
			  obs::solarShifter* _shifter, const units::InternalUnits& intUnits):
		    Rd(_Rd*intUnits.from_Kpc), zd(_zd*intUnits.from_Kpc), Rw(_Rw*intUnits.from_Kpc), Hw(_Hw*intUnits.from_Kpc), shifter(_shifter){
			from_Kpc = intUnits.from_Kpc;
			ptr = NULL; sp = NULL; cl = NULL;
			log_rho0 = 0;
			log_rho0 = log((dAvds/intUnits.from_Kpc)/dens(coord::toPosCyl(shifter->xyz())));
		}
		dustModel(potential::PtrDensity _ptr, const double dAvds,
			  obs::solarShifter* _shifter, const units::InternalUnits& intUnits):
		    ptr(_ptr), shifter(_shifter) {
			from_Kpc = intUnits.from_Kpc;
			bl = NULL; sp = NULL; cl = NULL;
			log_rho0 = 0;
			log_rho0 = log((dAvds/intUnits.from_Kpc)/ptr->density(coord::toPosCyl(shifter->xyz())));
		}
		~dustModel(void){
			if(bl!=NULL) bl->~vector();
			if(sp!=NULL) sp->~vector();
			if(cl!=NULL) cl->~vector();
		}
		double dens(const coord::PosCyl& p) const;
		void addBlob(const coord::PosCyl p,double Ampl,double rad){
			Blob nbl(p,Ampl,rad);
			if(bl==NULL){
				std::vector<Blob>* blobs = new std::vector<Blob>;
				bl=blobs;
			}
			bl->push_back(nbl);
		}
		void deleteBlob(int n=-1){
			if(bl!=NULL){
				if(n==-1) bl->erase(bl->end());
				else bl->erase(bl->begin()+n);
			}
		}
		void addSpiral(double Ampl,double phase,double alpha,int Narms,double kz){
			Spiral nsp(Ampl,phase,alpha,Narms,kz);
			if(sp==NULL){
				std::vector<Spiral>* spirals = new std::vector<Spiral>;
				sp=spirals;
			}
			sp->push_back(nsp);
		}
		void deleteSpiral(int n=-1){
			if(sp!=NULL){
				if(n==-1) sp->erase(sp->end());
				else sp->erase(sp->begin()+n);
			}
		}
		void addCloud(double Rc, double phic, double z0, double norm,const std::string& fname){
			JTcloud ncl(Rc,phic,z0,norm,fname);
			if(cl==NULL){
				std::vector<JTcloud>* clouds = new std::vector<JTcloud>;
				cl=clouds;
			}
			cl->push_back(ncl);
		}
		void deleteCloud(int n=-1){
			if(cl!=NULL){
				if(n==-1) cl->erase(cl->end());
				else cl->erase(cl->begin()+n);
			}
		}
		double A_V (const coord::PosCar) const;
		double A_B (const coord::PosCar p) const{
			return 1.324 * A_V(p);
		}
		double A_R (const coord::PosCar p) const{
			return 0.748 * A_V(p);
		}
		double A_H (const coord::PosCar p) const{
			return 0.175 * A_V(p);
		}
		double A_K (const coord::PosCar p) const{
			return 0.112 * A_V(p);
		}
};

} // namespace dust

namespace obs {
//class line of sight, specified by (l,b) and location of Sun.
//Extinction in BVRHK down it will be computed if a dustModel isspecified;
//otherwise zero extinction
class EXP los{
	private:
		math::CubicSpline extinct;
		dust::dustModel* dm;
	public: 
		obs::solarShifter sun;
		PosSky pos;// l, b in radians
		double cosl, cosb, sinl, sinb;
		los(const PosSky _pos, const obs::solarShifter _sun,
		    dust::dustModel* _dm=NULL);//angles in radians
		los(const coord::PosCar xyz, const obs::solarShifter _sun,
		    dust::dustModel* _dm=NULL);
		coord::PosCar xyz(const double s) const{// s in internal units
			double x = sun.xyz().x + s*cosb*cosl;
			double y = sun.xyz().y + s*cosb*sinl;
			double z = sun.xyz().z + s*sinb;
			coord::PosCar Car(x, y, z);
			return Car;
		}
		coord::PosCyl Rzphi(double s) const{// s in internal units
			coord::PosCyl Cyl(toPosCyl(xyz(s)));			
			return Cyl;
		}
		double s(const coord::PosCar& xyz) const{// s in internal units
			double s = pow_2(xyz.x-sun.xyz().x) + pow_2(xyz.y-sun.xyz().y) + pow_2(xyz.z-sun.xyz().z);
			return (s>0?  sqrt(s) : 0);
		}
		double s(const coord::PosCyl& p) const{
			return s(coord::toPosCar(p));
		}
		double sKpc(const coord::PosCar& xyz) const{// s in Kpc
			return s(xyz)/from_Kpc;
		}
		double sKpc(const coord::PosCyl& p) const{
			return sKpc(coord::toPosCar(p));
		}
		double A_V(const double sKpc) const{
			if(!dm) return 0;
			else{
				double A,dA;
				extinct.evalDeriv(sKpc,&A,&dA);
				return A;
			}
		}
		double A_B(const double sKpc) const{
			return 1.324 * A_V(sKpc);
		}
		double A_R(const double sKpc) const{
			return 0.748 * A_V(sKpc);
		}
		double A_H(const double sKpc) const{
			return 0.175 * A_V(sKpc);
		}
		double A_K(const double sKpc) const{
			return 0.112 * A_V(sKpc);
		}
		math::CubicSpline tab_extinct(void);
};
}// namespace obs
