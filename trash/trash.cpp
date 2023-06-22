const double dL=.05;
double Vxy2=pow_2(posvel.vR)+pow_2(posvel.vphi);
double E=.5*(Vxy2+pow_2(posvel.vz))+model.potential.value(posvel);
double DL=dL*potential::L_circ(model.potential,E);//uncertainty in Lz
if(posvel.R*fabs(posvel.vphi)<3*DL){
	double dV=DL/posvel.R;//change in vphi to give DL
	double Vxy=sqrt(Vxy2);;
	double psi=2*(math::random()-.5)*fmin(M_PI,dV/Vxy);
	double cs=cos(psi),sn=sin(psi);
	double sav=posvel.vR;
	posvel.vR  =sav*cs-posvel.vphi*sn;
	posvel.vphi=sav*sn+posvel.vphi*cs;
}
