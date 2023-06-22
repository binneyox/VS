#include "df_halo_RR.h"
#include "math_core.h"
#include <cmath>
#include <fstream>

namespace df{

EXP DoublePowerLaw::DoublePowerLaw(const DoublePowerLawParam &inparams) :
    par(inparams)
{
    // sanity checks on parameters
	if(!(par.norm>0))
		throw std::invalid_argument("DoublePowerLaw: normalization must be positive");
	if(!(par.J0>0))
		throw std::invalid_argument("DoublePowerLaw: break action J0 must be positive");
	if(par.Jcutoff<0)
		throw std::invalid_argument("DoublePowerLaw: truncation action Jcutoff must be non-negative");
	if(!(par.slopeOut>3) && par.Jcutoff==0)
		throw std::invalid_argument(
					    "DoublePowerLaw: mass diverges at large J (outer slope must be > 3)");
	if(!(par.slopeIn<3))
		throw std::invalid_argument(
					    "DoublePowerLaw: mass diverges at J->0 (inner slope must be < 3)");
	if(par.steepness<=0)
		throw std::invalid_argument("DoublePowerLaw: transition steepness parameter must be positive");
	if(par.cutoffStrength<=0)
		throw std::invalid_argument("DoublePowerLaw: cutoff strength parameter must be positive");
	if( par.coefJphiIn <=0 || par.coefJzIn <=0 || par.coefJphiIn +par.coefJzIn >=5 || 
	    par.coefJphiOut<=0 || par.coefJzOut<=0 || par.coefJphiOut+par.coefJzOut>=5 )
		throw std::invalid_argument(
					    "DoublePowerLaw: invalid weights in the linear combination of actions");
	if(fabs(par.rotFrac)>1)
		throw std::invalid_argument(
					    "DoublePowerLaw: amplitude of odd-Jphi component must be between -1 and 1");

}

EXP double DoublePowerLaw::value(const actions::Actions &J) const
{
    // linear combination of actions in the inner part of the model (for J<J0)
	double hJ  = (3-par.coefJphiIn-par.coefJzIn) * J.Jr + 0.7 * (par.coefJzIn * J.Jz +
		     par.coefJphiIn * fabs(J.Jphi));
    // linear combination of actions in the outer part of the model (for J>J0)
	double gJ  = (3-par.coefJphiOut-par.coefJzOut) * J.Jr + 0.7 * (par.coefJzOut* J.Jz +
		     par.coefJphiOut* fabs(J.Jphi));
    double val = par.norm / pow_3(2*M_PI * par.J0) *
		 math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		 math::pow(1 + math::pow(gJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
    if(par.rotFrac!=0)  // add the odd part
	    val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
    if(par.Jcutoff>0)   // exponential cutoff at large J
	    val *= exp(-math::pow(gJ / par.Jcutoff, par.cutoffStrength));
    return val;
}

EXP void DoublePowerLaw::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
//	const units::ExternalUnits& conv);
	strm << "norm\t\t" << par.norm*intUnits.to_Msun << '\n';
	strm << "J0\t\t" << par.J0*intUnits.to_Kpc_kms << '\n';
	strm << "Jcutoff\t\t" << par.Jcutoff*intUnits.to_Kpc_kms << '\n';
	strm << "Jphi0\t\t" << par.Jphi0*intUnits.to_Kpc_kms << '\n';
	strm << "slopeIn\t\t" << par.slopeIn << '\n';
	strm << "slopeOut\t" << par.slopeOut << '\n';
	strm << "steepness\t" << par.steepness << '\n';
	strm << "coefJphiIn\t" << par.coefJphiIn << '\n';
	strm << "coefJzIn\t" << par.coefJzIn << '\n';
	strm << "coefJphiOut\t" << par.coefJphiOut << '\n';
	strm << "coefJzOut\t" << par.coefJzOut << '\n';
	strm << "rotFrac\t\t" << par.rotFrac << '\n';
	strm << "cutoffStrength\t" << par.cutoffStrength << '\n';
}

EXP SoftDoublePowerLaw::SoftDoublePowerLaw(const SoftDoublePowerLawParam &inparams) :
    par(inparams)
{
    // sanity checks on parameters
	if(!(par.norm>0))
		throw std::invalid_argument("DoublePowerLaw: normalization must be positive");
	if(!(par.J0>0))
		throw std::invalid_argument("DoublePowerLaw: break action J0 must be positive");
	if(par.Jcutoff<0)
		throw std::invalid_argument("DoublePowerLaw: truncation action Jcutoff must be non-negative");
	if(!(par.slopeOut>3) && par.Jcutoff==0)
		throw std::invalid_argument(
					    "DoublePowerLaw: mass diverges at large J (outer slope must be > 3)");
	if(!(par.slopeIn<3))
		throw std::invalid_argument(
					    "DoublePowerLaw: mass diverges at J->0 (inner slope must be < 3)");
	if(par.steepness<=0)
		throw std::invalid_argument("DoublePowerLaw: transition steepness parameter must be positive");
	if(par.cutoffStrength<=0)
		throw std::invalid_argument("DoublePowerLaw: cutoff strength parameter must be positive");
	if( par.coefJphiIn <=0 || par.coefJzIn <=0 || par.coefJphiIn +par.coefJzIn >=3 || 
	    par.coefJphiOut<=0 || par.coefJzOut<=0 || par.coefJphiOut+par.coefJzOut>=3 )
		throw std::invalid_argument(
					    "DoublePowerLaw: invalid weights in the linear combination of actions");
	if(fabs(par.rotFrac)>1)
		throw std::invalid_argument(
					    "DoublePowerLaw: amplitude of odd-Jphi component must be between -1 and 1");
	set_beta();
}
double SoftDoublePowerLaw::g(const double hJ) const{
	double rat=par.h0/hJ;
	return pow(rat*rat-beta*rat+1,-.5*par.slopeIn);
}

EXP double SoftDoublePowerLaw::value(const actions::Actions &J) const
{
    // linear combination of actions in the inner part of the model (for J<J0)
	double hJ  = (3-par.coefJphiIn-par.coefJzIn) * J.Jr + 0.7 * (par.coefJzIn * J.Jz +
		     par.coefJphiIn * fabs(J.Jphi));
    // linear combination of actions in the outer part of the model (for J>J0)
	double gJ  = (3-par.coefJphiOut-par.coefJzOut) * J.Jr + 0.7 * (par.coefJzOut* J.Jz +
		     par.coefJphiOut* fabs(J.Jphi));
	double val =  par.norm / pow_3(2*M_PI * par.J0) *
		    math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		    math::pow(1 + math::pow(gJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	if(par.Jcutoff>0)   // exponential cutoff at large J
		val *= exp(-math::pow(gJ / par.Jcutoff, par.cutoffStrength));
	return g(hJ) * val;
}

EXP void SoftDoublePowerLaw::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
	strm << "norm\t\t" << par.norm*intUnits.to_Msun << '\n';
	strm << "J0\t\t" << par.J0*intUnits.to_Kpc_kms << '\n';
	strm << "Jcutoff\t\t" << par.Jcutoff*intUnits.to_Kpc_kms << '\n';
	strm << "Jphi0\t\t" << par.Jphi0*intUnits.to_Kpc_kms << '\n';
	strm << "slopeIn\t\t" << par.slopeIn << '\n';
	strm << "slopeOut\t" << par.slopeOut << '\n';
	strm << "steepness\t" << par.steepness << '\n';
	strm << "coefJphiIn\t" << par.coefJphiIn << '\n';
	strm << "coefJzIn\t" << par.coefJzIn << '\n';
	strm << "coefJphiOut\t" << par.coefJphiOut << '\n';
	strm << "coefJzOut\t" << par.coefJzOut << '\n';
	strm << "rotFrac\t\t" << par.rotFrac << '\n';
	strm << "cutoffStrength\t" << par.cutoffStrength << '\n';
	strm << "coreAction\t" << par.h0 << '\n';
}
double SoftDoublePowerLaw::diff(const double hJ) const{
	double val = 
		    math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		    math::pow(1 + math::pow(hJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
	return hJ*hJ*val*(1-g(hJ));
}
double SoftDoublePowerLaw::intDiff(const double beta_) const{
	SoftDoublePowerLaw CB(par,beta_);
	ColeBinneyDiff Diff(CB);
	double ans=math::integrateGL(Diff,1e-5,10*par.h0,10);
	return ans;
}
void SoftDoublePowerLaw::set_beta(void){
	beta=0;
	ColeBinneyInt Int(*this);
	beta=math::findRoot(Int,0,1,1e-5);
}

}// namespace df
