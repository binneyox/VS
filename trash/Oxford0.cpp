
/// A class and a function required by Oxford & Bologna
void OB_syst::eval(const double t,const double s[],double dsdt[]) const{
	dsdt[0]=-1;//t=const-el
	double a=.5*(k+1), b=.5*(k-1);
	double L=.5*s[0]-b*s[2]+.5*sqrt(fmax(0,4*b*s[2]*(a*s[2]-s[0])+s[0]*s[0]));
	double c=L/(L+s[1]);
//	double c=s[0]/(s[0]+s[1]);
	if(OB==0)
		dsdt[1]=Oxford_rat(s,c,gamma0,gamma1)*exp(coefL*sin(.5*M_PI*c));
	else
		dsdt[1]=Bologna_rat(s,c,gamma0)*exp(coefL*sin(.5*M_PI*c));
	dsdt[2]=0;
}

double Oxford_rat(const double s[],double c,double gamma0,double gamma1){
	double jt=s[0]+1.5*s[1] + gamma0;
	if(jt<=0) return .5;
	double powj=pow(jt,gamma1);
	return .5*(1+c*powj/(powj+1));
}

double Bologna_rat(const double s[],double c,double gamma){
	return .5*(1+c*tanh(s[0]+1.5*s[1]+gamma));
}

double get_h(int OB,double ic[],double coefL,double k,double gamma0,double gamma1){
	OB_syst system(OB,coefL,k,gamma0,gamma1);
	math::OdeSolverDOP853 solver(system);
	solver.init(ic);
	double t=0;
	int j=0;
	while(j<20 && solver.getSol(t,0)>0){
		solver.doStep(0); t=solver.getTime();
		j++;
	}
	if(j==20){
		double s[2]={solver.getSol(t,0),solver.getSol(t,1)};
		printf("%g (%g %g) (%g %g) %g\n",t,ic[0],ic[0]/ic[1],
		       s[0],s[1],Oxford_rat(s,0,gamma0,gamma1)); exit(0);
	}
	return solver.getSol(ic[0],1);
}

Oxford::Oxford(const OxfordParam &inparams) :  par(inparams)
{
	norm=1;
    // sanity checks on parameters
	if(par.mass == 0)
		throw std::invalid_argument("Oxford: mass must be non-zero");
	if(!(par.J0>0))
		throw std::invalid_argument("Oxford: break action J0 must be positive");
	if(par.Jcutoff<0)
		throw std::invalid_argument("Oxford: truncation action Jcutoff must be non-negative");
	if(!(par.slopeOut>3) && par.Jcutoff==0)
		throw std::invalid_argument(
					    "Oxford: mass diverges at large J (outer slope must be > 3)");
	if(!(par.slopeIn<3))
		throw std::invalid_argument(
					    "Oxford: mass diverges at J->0 (inner slope must be < 3)");
	if(par.steepness<=0)
		throw std::invalid_argument(
					    "Oxford: transition steepness parameter must be positive");
	if(par.cutoffStrength<=0)
		throw std::invalid_argument("Oxford: cutoff strength parameter must be positive");
	if( par.L0 <=0 )
		throw std::invalid_argument("Oxford: L0 must be >0");
	if( par.kIn < 1 || par.kOut < 1  )
		throw std::invalid_argument("Oxford: Must have k >= 1");
	if( par.coefLin > 1.5 || par.coefLout > 1.5  )
		throw std::invalid_argument("Oxford: Must have coefL < 1.5");
	if(fabs(par.rotFrac)>1)
		throw std::invalid_argument(
					    "Oxford: amplitude of odd-Jphi component must be between -1 and 1");
	set_beta();
	printf("beta: %f\n",beta);
	/*printf("beta,mass: %f %f %g\n",beta,norm,totalMass());*/
	if(par.mass>0)
		norm*=par.mass/totalMass();
	else
		norm = -par.mass/pow_3(par.J0);

}
void Oxford::set_norm(double fac){
	norm *= fac;
}
double Oxford::g(const double hJ) const{
	double rat=par.Jcore/hJ;
	return pow(rat*rat-beta*rat+1,-.5*par.slopeIn);
}
double Oxford::value(const actions::Actions &J) const
{
	double fJphi=fabs(J.Jphi), L=J.Jz+fJphi;
	double Jt=6*(1.5*J.Jr+L)/par.L0, zeta=Jt/(1+Jt);
	double k=par.kIn + zeta*(par.kOut-par.kIn);
	//Now that we have good estimate of k compute cL and xi
	double a=.5*(k+1), b=.5*(k-1);
	double cL=(a*J.Jz+b*J.Jz*fJphi/L+fJphi);
#ifdef DE
	double ic[3]={cL/par.L0,J.Jr/par.L0,J.Jz/par.L0};
	double hJ=par.L0*get_h(0,ic,par.coefLin,k,par.gamma0,par.gamma1);
	double gJ;
	if(par.coefLout==par.coefLin) gJ=hJ;
	else gJ=par.L0*get_h(0,ic,par.coefLout,k,par.gamma0,par.gamma1);
#else
	double powJt=pow((1.5*J.Jr+L)/par.L0,par.gamma1), xi=powJt/(1+powJt);
	double c=L/(L+J.Jr), stheta=sin(.5*M_PI*c), fac=exp(-par.coefLin*stheta);
	double hJ=J.Jr * fac + .5*(1+xi*c)/fac * cL;
	fac=exp(-par.coefLout*stheta);
	double gJ=J.Jr * fac + .5*(1+xi*c)/fac * cL;
#endif
	double val =  norm *
		      math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		      math::pow(1 + math::pow(gJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	if(par.Jcutoff>0)   // exponential cutoff at large J
		val *= exp(-math::pow((J.Jr+J.Jz+fabs(J.Jphi)) / par.Jcutoff, par.cutoffStrength));
	return g(hJ) * val;
}
void Oxford::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
	strm << "mass\t\t" << par.mass*intUnits.to_Msun << '\n';
	strm << "J0\t\t" << par.J0*intUnits.to_Kpc_kms << '\n';
	strm << "Jcutoff\t\t" << par.Jcutoff*intUnits.to_Kpc_kms << '\n';
	strm << "Jphi0\t\t" << par.Jphi0*intUnits.to_Kpc_kms << '\n';
	strm << "Jcore\t\t" << par.Jcore*intUnits.to_Kpc_kms << '\n';
	strm << "slopeIn\t\t" << par.slopeIn << '\n';
	strm << "slopeOut\t" << par.slopeOut << '\n';
	strm << "steepness\t" << par.steepness << '\n';
	strm << "coefLin\t" << par.coefLin << '\n';
	strm << "coefLout\t" << par.coefLout << '\n';
	strm << "kIn\t" << par.kIn << '\n';
	strm << "kOut\t" << par.kOut << '\n';
	strm << "L0\t" << par.L0 << '\n';
	strm << "gamma0\t" << par.gamma0 << '\n';
	strm << "rotFrac\t\t" << par.rotFrac << '\n';
	strm << "cutoffStrength\t" << par.cutoffStrength << '\n';
}
void Oxford::tab_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
	strm << "mass\t J0\t Jcutoff\t Jphi0\t Jcore\t L0\t slopeIn\t slopeOut\t steepness\t coefLin\n";
	strm << "kIn\t coefLout\t kOut\t rotFrac\t cutoffStrength\t gamma\n";
	strm << par.mass*intUnits.to_Msun << " & "
			<< par.J0*intUnits.to_Kpc_kms << " & "
			<< par.Jcutoff*intUnits.to_Kpc_kms << " & "
			<< par.Jphi0*intUnits.to_Kpc_kms << " & "
			<< par.Jcore*intUnits.to_Kpc_kms << " & "
			<< par.L0*intUnits.to_Kpc_kms << " & "
			<< par.slopeIn << " & "
			<< par.slopeOut << " & "
			<< par.coefLin << " & "
			<< par.kIn << " & "
			<< par.coefLout << " & "
			<< par.kOut << " & "
			<< par.rotFrac << " & "
			<< par.cutoffStrength << " & "
			<< par.gamma0 << " \\cr \n";
}
double Oxford::diff(const double hJ) const{
	double val = 
		    math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		    math::pow(1 + math::pow(hJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
	return hJ*hJ*val*(1-g(hJ));
}
double Oxford::intDiff(const double beta_) const{
	Oxford CB(par,beta_);
	OxfordDiff Diff(CB);
	return math::integrateGL(Diff,1e-5,10*par.Jcore,10);
}
void Oxford::set_beta(void){
	beta=0;
	OxfordInt Int(*this);
	beta=math::findRoot(Int,0,1,1e-5);
}

Bologna::Bologna(const BolognaParam &inparams) :
    par(inparams)
{
	norm=1;
    // sanity checks on parameters
	if(par.mass==0)
		throw std::invalid_argument("Bologna: normalization must be non-zero");
	if(!(par.J0>0))
		throw std::invalid_argument("Bologna: break action J0 must be positive");
	if(par.alpha < 0 || par.kIn < 0 || par.kOut < 0) 
		throw std::invalid_argument("Bologna: steepness and k parameters must be positive");
	if(par.mass>0)
		norm *= par.mass/totalMass();
	else
		norm = -par.mass/pow_3(par.J0);
}

double Bologna::value(const actions::Actions &J) const {
	double fJphi=fabs(J.Jphi), L=J.Jz+fJphi;
	double Jt=2*(1.5*J.Jr+L)/par.L0, zeta=Jt/(1+Jt);
	double k=par.kIn + zeta*(par.kOut-par.kIn);
	//Now that we have good estimate of k compute cL and xi
	double a=.5*(k+1), b=.5*(k-1);
	double cL=a*J.Jz+b*J.Jz*fJphi/(L+par.L1)+fJphi;
	double ic[2]={cL/par.L0,J.Jr/par.L0};
	double gJ=par.L0*get_h(1,ic,par.coefL,k,par.gamma,NAN);
/*	double powJt=pow((1.5*J.Jr+L)/par.L0,par.delta), xi=tanh(powJt);
	double c=L/(L+J.Jr);
	if(par.coefL!=0){
		double fac=exp(-par.coefL*sin(.5*M_PI*c));
		gJ=J.Jr * fac + .5*(1+xi*c)/fac * cL;
	} else {
		gJ=J.Jr + .5*(1+xi*c) * cL;
	}*/
	double val = norm * exp(-pow( gJ / par.J0, par.alpha));
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi);
	return val; 
}
