#include "math_ODE.h"
#include "obs.h"
#include "dust.h"
#include <cmath>


namespace obs{


//Three constants defining Galactic coordinates
#define	ap (192.85948*torad)
#define dp (27.12825*torad)
#define l0 (122.931918*torad)

EXP PosSky from_RAdec(double RA,double dec){//in and out in degrees
	const double torad=M_PI/180.;
	const double alpha=RA*torad, delta=dec*torad;
	double b=asin(sin(dp)*sin(delta)+cos(dp)*cos(delta)*cos(alpha-ap));
	double l=l0-atan2(cos(delta)*sin(alpha-ap),
			  cos(dp)*sin(delta)-sin(dp)*cos(delta)*cos(alpha-ap));
	PosSky p(l/torad, b/torad);
	return p;
}
EXP PosSky to_RAdec(double l,double b){//in and out in degrees
	const double torad=M_PI/180.;
	l*=torad; b*=torad;
	double delta=asin(sin(dp)*sin(b)+cos(dp)*cos(b)*cos(l0-l));
	double alpha=ap+atan2(cos(b)*sin(l0-l),
			      cos(dp)*sin(b)-sin(dp)*cos(b)*cos(l0-l));
	PosSky p(alpha/torad, delta/torad);
	return p;
}
//muRA =RA dot cos(delta) & similarly for mul
EXP PosVelSky from_muRAdec(double RA,double dec,double muRA,double mudec){
	const double torad=M_PI/180.;
	PosSky poslb(from_RAdec(RA,dec));
	RA*=torad; dec*=torad;
	double lCP=l0,RA_GP=ap,decGP=dp;
	double sb,cb,l=poslb.l*torad,sd,cd;
	math::sincos(poslb.b*torad,sb,cb);
	math::sincos(dec,sd,cd);
	double mub=((sin(decGP)*cd-cos(decGP)*sd*cos(RA-RA_GP))*mudec
		    -cos(decGP)*sin(RA-RA_GP)*muRA)/cb;
	double mul;
	if(fabs(cos(lCP-l))>fabs(sin(lCP-l)))
		mul=(sd*sin(RA-RA_GP)*mudec-cos(RA-RA_GP)*muRA
		     -sb*sin(lCP-l)*mub)/cos(lCP-l);
	else
		mul=((cos(decGP)*cd+sin(decGP)*sd*cos(RA-RA_GP))*mudec
		     +sin(decGP)*sin(RA-RA_GP)*muRA+sb*cos(lCP-l)*mub)/sin(lCP-l);
	VelSky pm(mul,mub);
	PosVelSky p(poslb,pm);
	return p;
}

EXP solarShifter::solarShifter(const units::InternalUnits& intUnits,
				    coord::PosVelCar* _Sun){
	torad = M_PI/180.;
	from_Kpc = intUnits.from_Kpc;
	from_kms = intUnits.from_kms;
	from_mas_per_yr = intUnits.from_mas_per_yr;
	if(_Sun) Sun=*_Sun;
	else Sun=coord::PosVelCar(-8.27 * from_Kpc,0,0.025 * from_Kpc,
			   12 * from_kms,249 * from_kms,7 * from_kms);
}
EXP coord::PosCar solarShifter::toCar(const PosSky pos, double sKpc) const{
	double cb,sb,cl,sl;
	math::sincos(torad*pos.b,sb,cb);
	math::sincos(torad*pos.l,sl,cl);
	double s = sKpc*from_Kpc;
	return coord::PosCar(Sun.x+s*cb*cl, Sun.y+s*cb*sl, Sun.z+s*sb);
}	
EXP coord::PosVelCar solarShifter::toCar(const PosSky pos, double sKpc,
					 const VelSky pm, double Vlos_kms) const{
	double cb,sb,cl,sl;
	math::sincos(torad*pos.b,sb,cb);
	math::sincos(torad*pos.l,sl,cl);
	double s = sKpc*from_Kpc, Vlos = Vlos_kms*from_kms;
	double X = Sun.x+s*cb*cl;
	double Y = Sun.y+s*cb*sl;
	double Z = Sun.z+s*sb;
	double U = Sun.vx+Vlos*cb*cl-s*(sb*cl*pm.mub+sl*pm.mul)*from_mas_per_yr;
	double V = Sun.vy+Vlos*cb*sl-s*(sb*sl*pm.mub-cl*pm.mul)*from_mas_per_yr;
	double W = Sun.vz+Vlos*sb+s*cb*pm.mub*from_mas_per_yr;
	coord::PosVelCar O(X,Y,Z,U,V,W);//xv in internal units
	return O;
}
EXP coord::PosVelCyl solarShifter::toCyl(const PosSky pos, const double sKpc, const VelSky pm, double Vlos_kms) const{
	return toPosVelCyl(toCar(pos,sKpc,pm,Vlos_kms));
}
EXP PosSky solarShifter::toSky(const coord::PosCar p, double& sKpc) const{
	double xhel = p.x-Sun.x, yhel=p.y-Sun.y, zhel=p.z-Sun.z;//position wrt Sun
	double s = sqrt(xhel*xhel+yhel*yhel+zhel*zhel);
	sKpc = s/from_Kpc; //return in kpc
	double b=asin(zhel/s), l=atan2(yhel,xhel);
	return PosSky(l/torad, b/torad);
}
EXP PosVelSky solarShifter::toSky(const coord::PosVelCar pv, double& sKpc, double& Vlos_kms) const{
	double xhel = pv.x-Sun.x, yhel=pv.y-Sun.y, zhel=pv.z-Sun.z;//position wrt Sun
	double s = sqrt(xhel*xhel+yhel*yhel+zhel*zhel);
	sKpc = s/from_Kpc; //return in kpc
	double Vhelx = pv.vx-Sun.vx, Vhely = pv.vy-Sun.vy, Vhelz = pv.vz-Sun.vz;//V wrt Sun
	double Vlos = (xhel*Vhelx + yhel*Vhely + zhel*Vhelz)/s;
	Vhelx -= Vlos*xhel/s; Vhely -= Vlos*yhel/s; Vhelz -= Vlos*zhel/s; //sky velocity
	Vlos_kms = Vlos/from_kms; //return Vlos in kms
	double b=asin(zhel/s), l=atan2(yhel,xhel);
	double sb,cb,sl,cl;
	math::sincos(b,sb,cb); math::sincos(l,sl,cl);
	double mub = Vhelz/(s*cb);
	double mul;
	if(fabs(sl)>fabs(cl))
		mul=-(Vhelx+Vhelz*sb*cl/cb)/(s*sl);
	else
		mul=(Vhely+Vhelz*sb*sl/cb)/(s*cl);
	PosVelSky p(l/torad, b/torad, mul/from_mas_per_yr, mub/from_mas_per_yr);
	return p;
}

EXP VelSky solarShifter::toPM(const coord::PosVelCar pv, double& Vlos_kms) const{
	double sKpc;
	PosVelSky p = toSky(pv, sKpc, Vlos_kms);
	return p.pm;
}


EXP los::los(const PosSky _pos, const solarShifter _sun, dust::dustModel* _dm)//angles in radians// added for sKpc
: pos(_pos), cosl(cos(pos.l)), cosb(cos(pos.b)), sinl(sin(pos.l)), sinb(sin(pos.b)),
sun(_sun), from_Kpc(sun.fromKpc()), dm(_dm)
{
	if(dm) extinct = tab_extinct();
}

EXP los::los(const coord::PosCar xyz, const solarShifter _sun, dust::dustModel* _dm)
: sun(_sun), from_Kpc(Sun.fromKpc()), dm(_dm)// added for sKpc
{
	double xhel = xyz.x-sun.xyz().x, yhel=xyz.y-sun.xyz().y, zhel=xyz.z-sun.xyz().z;//position wrt Sun
	double s = sqrt(xhel*xhel+yhel*yhel+zhel*zhel);
	pos.b = asin(zhel/s); pos.l = atan2(yhel,xhel);
	cosl = cos(pos.l); sinl = sin(pos.l);
	cosb = cos(pos.b); sinb = sin(pos.b);
	if(dm) extinct = tab_extinct(from_Kpc); //sKpc
}

//helper class for tabulating extinctions: we integrate dN/ds = rho
//along LOS. Distances in Kpc
class A_syst: public math::IOdeSystem{
	private:
		const obs::los Los;
		const dust::dustModel* dm; 
	public:
		A_syst(const obs::los _los, dust::dustModel* _dm)
				: Los(_los), dm(_dm){}
		virtual void eval(const double sKpc,const double A[], double dAds[]) const{
			dAds[0] = dm->dens(Los.Rzphi(sKpc*from_Kpc));
		}
		virtual unsigned int size() const {return 1;}
};


EXP math::CubicSpline los::tab_extinct(void){
	std::vector<double> dists, Ns;
	const int ns=3;
	A_syst syst(los(pos, sun), dm);
	math::OdeSolverDOP853 solver(syst);
	double sKpc=0, rho=fmax(0,dm->dens(coord::toPosCyl(sun.xyz()))), N=0, rhoPeak=0; int iz=0;
	solver.init(&N);
	dists.push_back(sKpc); Ns.push_back(N);
	while(iz<200 && rho >= 0.01*rhoPeak){
		solver.doStep();
		double s_last=sKpc, ds=(solver.getTime()-s_last)/(double)ns;
		for(int i=0; i<ns; i++){// store ns points over step
			sKpc+=ds;
			dists.push_back(sKpc); Ns.push_back(solver.getSol(sKpc,0));
			//printf("%f %f %f\n",s,0.175*solver.getSol(sKpc,0),rho);
		}
		rho=dm->dens(Rzphi(sKpc*from_Kpc)); rhoPeak= rho>rhoPeak? rho : rhoPeak;
		iz++;
	}
	return math::CubicSpline(dists,Ns);
}

} // namespace obs