#include "coord_jjb.h"

#define	ap (192.85948*torad)
#define dp (27.12825*torad)
#define l0 (122.931918*torad)

namespace coord {
EXP solarShifter::solarShifter(const units::InternalUnits intUnits, PosVelCar* _Sun){
	if(_Sun) Sun=*_Sun;
	else Sun=PosVelCar(-8.27,0,0.025,12,249,7);
	torad = M_PI/180.;
	from_Kpc = intUnits.from_Kpc;
	from_kms = intUnits.from_kms;
	from_mas_per_yr = intUnits.from_mas_per_yr;
	Sun.x *= from_Kpc;
	Sun.y *= from_Kpc;
	Sun.z *= from_Kpc;
	Sun.vx *= from_kms;
	Sun.vy *= from_kms;
	Sun.vz *= from_kms;
}
EXP PosVelCar solarShifter::toCar(const PosSky pos,double s,const VelSky pm,double vlos){
	double cb,sb,cl,sl;
	math::sincos(torad*pos.b,sb,cb);
	math::sincos(torad*pos.l,sl,cl);
	s *= from_Kpc; vlos *= from_kms;
	double X = Sun.x+s*cb*cl;
	double Y = Sun.y+s*cb*sl;
	double Z = Sun.z+s*sb;
	double U = Sun.vx+vlos*cb*cl-s*(sb*cl*pm.mub+cb*sl*pm.mul)*from_mas_per_yr;
	double V = Sun.vy+vlos*cb*sl-s*(sb*sl*pm.mub-cb*cl*pm.mul)*from_mas_per_yr;
	double W = Sun.vz+vlos*sb+s*cb*pm.mub*from_mas_per_yr;
	PosVelCar O(X,Y,Z,U,V,W);//xv in internal units
	return O;
}
EXP PosVelCyl solarShifter::toCyl(const PosSky pos,const double s,const VelSky pm,double vlos){
	return toPosVelCyl(toCar(pos,s,pm,vlos));
}
EXP VelSky solarShifter::toPM(const PosVelCar pv,double& Vlos){
	double xhel = pv.x-Sun.x, yhel=pv.y-Sun.y, zhel=pv.z-Sun.z;//position wrt Sun
	double s = sqrt(xhel*xhel+yhel*yhel+zhel*zhel);
	double Vhelx = pv.vx-Sun.vx, Vhely = pv.vy-Sun.vy, Vhelz = pv.vz-Sun.vz;//V wrt Sun
	Vlos = (xhel*Vhelx + yhel*Vhely + zhel*Vhelz)/s;
	Vhelx -= Vlos*xhel/s; Vhely -= Vlos*yhel/s; Vhelz -= Vlos*zhel/s;//sky velocity
	double b=asin(zhel/s), l=atan2(yhel,xhel);
	double sb,cb,sl,cl;
	math::sincos(b,sb,cb); math::sincos(l,sl,cl);
	double mub = Vhelz/(s*cb);
	double mul;
	if(fabs(sl)>fabs(cl))
		mul=-(Vhelx+Vhelz*sb*cl/cb)/(s*cb*sl);
	else
		mul=(Vhely+Vhelz*sb*sl/cb)/(s*cb*cl);
	VelSky p(mul/from_mas_per_yr,mub/from_mas_per_yr);
	return p;
}

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

}// namespace