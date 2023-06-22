#include <stdio.h>
#include "coord.h"
#include "units.h"

const units::InternalUnits intUnits(2.7183*units::Kpc, 3.1416*units::Myr);

int main(void){
	double RA=182.8756,dec=63.54138,s=0.49,l,b,mul,mub,U,V,W,torad=acos(-1)/180;
	coord::PosSky pos(coord::from_RAdec(RA,dec));
	printf("(RA,dec) = (%3.0f,%3.0f) => (l,b) = (%3.0f,%3.0f)\n",RA,dec,pos.l,pos.b);
	coord::PosSky ra(coord::to_RAdec(pos.l,pos.b));
	printf("(l,b) = (%3.0f,%3.0f) => (RA,dec) = (%3.0f,%3.0f)\n",pos.l,pos.b,ra.l,ra.b);

	coord::PosVelSky pv(coord::from_muRAdec(RA,dec,3.23463,-6.2911));
	printf("l b: (%f %f) mu_l mu_b: (%f %f)\n",pv.pos.l, pv.pos.b, pv.pm.mul,pv.pm.mub);
	coord::PosVelCar Sun(-8.27,0,0.025,12,249,7);
	coord::solarShifter sun(intUnits);
	double vlos=-25.57;
	coord::PosVelCar XV(sun.toCar(pv.pos,s,pv.pm,vlos));
	double sp=sqrt(pow(XV.x*intUnits.to_Kpc-Sun.x,2)+pow(XV.y*intUnits.to_Kpc,2)+pow(XV.z*intUnits.to_Kpc-Sun.z,2));
	printf("s, sp: %f %f\n",s,sp);
	double Vlos;
	coord::VelSky mu(sun.toPM(XV,Vlos));
	printf("mu, mp: (%f %f) (%f %f)\n",pv.pm.mul,mu.mul,pv.pm.mub,mu.mub);
	printf("Vlos vlos: %f %f\n",vlos,Vlos);
/*	double pmRAcosd,PMdec;
	while(RA<400){
		scanf("%lf %lf %lf %lf",&RA,&dec,&pmRAcosd,&PMdec);
		pv=coord::from_muRAdec(RA,dec,pmRAcosd,PMdec);
		printf("%f %f %f %f\n",pv.pos.l,pv.pos.b,pv.pm.mul,pv.pm.mub);
	}*/
}
