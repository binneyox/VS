#include "df_disk_RR.h"
#include <cmath>
//#include <io.h>
#include <fstream>

namespace df{

namespace{  // internal

/// compute the average of DF over stellar age:
/// ( \int_0^1 dt B^2(t) \exp[ t/t_0 - A*B(t) ] ) / ( \int_0^1 dt \exp[ t/t_0 ] ),
/// where B(t) = ( (t + t_1) / (1 + t_1) )^{-2\beta}
double averageOverAge(double Ar, double Az, double sigmabirth, const AgeVelocityDispersionParam& par)
{
    if((par.beta_r == 0 && par.beta_z == 0) || par.sigmabirth == 1 || !isFinite(par.Tsfr))
        return exp(-Ar-Az);
    // if we have a non-trivial age-velocity dispersion relation,
    // then we need to integrate over sub-populations convolved with star formation history
    static const int NT = 5;      // number of points in quadrature rule for integration over age
    static const double qx[NT] =  // nodes of quadrature rule
    { 0.04691007703066802, 0.23076534494715845, 0.5, 0.76923465505284155, 0.95308992296933198 };
    static const double qw[NT] =  // weights of quadrature rule
				{ 0.11846344252809454, 0.23931433524968324, 64./225, 0.23931433524968324, 0.11846344252809454 };
    double s_r = std::pow(sigmabirth, 1./par.beta_r), s_z = std::pow(sigmabirth, 1./par.beta_z); 
    double integ=0, norm = 0;
    for(int i=0; i<NT; i++) {
        // t is the lookback time (stellar age) measured in units of galaxy time (ranges from 0 to 1)
        double t = qx[i];
        // star formation rate exponentially increases with look-back time
        double weight = exp(t / par.Tsfr) * qw[i];
        // velocity dispersions {sigma_r, sigma_z} scale as  [ t + s * (1-t) ]^beta
	double alpha_r = std::pow(t + (1-t) * s_r, 2*par.beta_r);
	double alpha_z = std::pow(t + (1-t) * s_z, 2*par.beta_z);
        //double multsq = std::pow(t + (1-t) * s, -2*par.beta);  // multiplied by sigma^-2
        integ += weight * exp(-Ar/alpha_r - Az/alpha_z) / (alpha_r*alpha_z);
        norm  += weight;
    }
    return integ / norm;
}
}

EXP QuasiIsothermal::QuasiIsothermal(const QuasiIsothermalParam &params, const potential::Interpolator& freqs) :
    par(params), freq(freqs)
{
    // sanity checks on parameters
    if(!(par.Sigma0>0))
        throw std::invalid_argument("QuasiIsothermal: surface density Sigma0 must be positive");
    if(!(par.Rdisk>0))
        throw std::invalid_argument("QuasiIsothermal: disk scale radius Rdisk must be positive");
    if(!(par.sigmar0>0))
        throw std::invalid_argument("QuasiIsothermal: velocity dispersion sigmar0 must be positive");
    if(!(par.Rsigmar>0))
        throw std::invalid_argument("QuasiIsothermal: velocity scale radius Rsigmar must be positive");
    if(!( (par.Hdisk>0) ^ (par.sigmaz0>0 && par.Rsigmaz>0) ))
        throw std::invalid_argument("QuasiIsothermal: should have either "
            "Hdisk>0 to assign the vertical velocity dispersion from disk scaleheight, or "
            "Rsigmaz>0, sigmaz0>0 to make it exponential in radius");
    if(par.Hdisk<0 || par.sigmaz0<0 || par.Rsigmaz<0)  // these are optional but non-negative
        throw std::invalid_argument("QuasiIsothermal: parameters cannot be negative");
    if(par.sigmabirth<=0 || par.sigmabirth>1)
        throw std::invalid_argument("QuasiIsothermal: invalid value for velocity dispersion at birth");
}

EXP double QuasiIsothermal::value(const actions::Actions &J) const
{
    // obtain the radius of in-plane motion with the given "characteristic" angular momentum
    double Rcirc = freq.R_from_Lz(sqrt(pow_2(par.Jmin) +
        pow_2(fabs(J.Jphi) + par.coefJr * J.Jr + par.coefJz * J.Jz)) );
    if(Rcirc > 20 * par.Rdisk)
        return 0;   // we're too far out, DF is negligibly small
    double kappa, nu, Omega;   // characteristic epicyclic freqs
    freq.epicycleFreqs(Rcirc, kappa, nu, Omega);
    // surface density follows an exponential profile in radius
    double Sigma = par.Sigma0 * exp( -Rcirc / par.Rdisk );
    // squared radial velocity dispersion is exponential in radius
    double sig_fac = exp ( -Rcirc / par.Rsigmar );
    double sigmabirth= sig_fac * par.sigmabirth;
    double sigmarsq = pow_2(par.sigmar0 * sig_fac ) + pow_2(par.sigmamin);
    // squared vertical velocity dispersion computed by either of the two methods: 
    double sigmazsq = pow_2(par.sigmamin) + (par.Hdisk>0 ?
        2 * pow_2(nu * par.Hdisk) :     // keep the disk thickness approximately equal to Hdisk, or
        pow_2(par.sigmaz0 * exp ( -Rcirc / par.Rsigmaz ) ) );  // make sigmaz exponential in radius
    // suppression factor for counterrotating orbits
    double negJphi = J.Jphi>0 ? 0. : 2*Omega * J.Jphi;
    double result = 1./(2*M_PI*M_PI) * Sigma * nu * Omega / (kappa * sigmarsq * sigmazsq) *
		    averageOverAge( (kappa * J.Jr - negJphi) / sigmarsq, nu * J.Jz / sigmazsq,
				    sigmabirth, par);
    return isFinite(result) ? result : 0;
}

EXP void QuasiIsothermal::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
//	const units::ExternalUnits& conv);
	strm << "Sigma0\t" << par.Sigma0*intUnits.to_Msun_per_pc2 << '\n';
	strm << "Rdisk\t" << par.Rdisk*intUnits.to_Kpc << '\n';
	strm << "Hdisk\t" << par.Hdisk*intUnits.to_Kpc << '\n';
	strm << "sigmar0\t\t" << par.sigmar0*intUnits.to_kms << '\n';
	strm << "sigmaz0\t" << par.sigmaz0*intUnits.to_kms << '\n';
	strm << "Rsigmar\t" << par.Rsigmar*intUnits.to_Kpc << '\n';
	strm << "Rsigmaz\t" << par.Rsigmaz*intUnits.to_Kpc << '\n';
	strm << "coefJr\t" << par.coefJr*intUnits.to_Kpc_kms << '\n';
	strm << "coefJz\t" << par.coefJz*intUnits.to_Kpc_kms << '\n';
	strm << "beta_r\t" << par.beta_r << '\n';
	strm << "beta_z\t" << par.beta_z << '\n';
	strm << "Tsfr\t" << par.Tsfr*intUnits.to_Gyr << '\n';
	strm << "sigmabirth\t" << par.sigmabirth << '\n';
}


EXP Exponential::Exponential(const ExponentialParam& params) :
    par(params)
{
	if(!(par.norm>0))
		throw std::invalid_argument("Exponential: overall normalization must be positive");
	if(!(par.Jr0>0) || !(par.Jz0>0) || !(par.Jphi0>0))
		throw std::invalid_argument("Exponential: scale actions must be positive");
	if(par.sigmabirth<=0 || par.sigmabirth>1)
		throw std::invalid_argument("Exponential: invalid value for velocity dispersion at birth");
	grad0=1/par.addJden-1/par.Jphi0-(par.pr+par.pz)/par.addJvel;
	if(grad0<0)
		throw std::invalid_argument("Exponential: invalid negative value for grad0: f diverges for Jphi<0");
}

EXP double Exponential::value(const actions::Actions &J) const
{
	double Jp =	J.Jphi<=0 ? 0 : J.Jphi; 
	double Jden = fabs(Jp) + par.addJden;
	double Jvel = fabs(Jp) + par.addJvel;
	double Jvelpr=pow(par.Jphi0/Jvel,par.pr);
	double Jvelpz=pow(par.Jphi0/Jvel,par.pz);
	double DF = par.norm / TWO_PI_CUBE / (par.Jr0 * par.Jz0) *
			Jvelpr*Jvelpz/pow_2(par.Jphi0) * Jden * exp(-Jden / par.Jphi0) *
		    exp(- Jvelpr * J.Jr / par.Jr0 - Jvelpz * J.Jz / par.Jz0);
	if(J.Jphi<0) return DF * exp(J.Jphi*grad0);
	else return DF;
}

EXP void Exponential::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
//	const units::ExternalUnits& conv);
	strm << "norm\t" << par.norm*intUnits.to_Msun << '\n';
	strm << "Jr0\t" << par.Jr0*intUnits.to_Kpc_kms << '\n';
	strm << "Jz0\t" << par.Jz0*intUnits.to_Kpc_kms << '\n';
	strm << "Jphi0\t" << par.Jphi0*intUnits.to_Kpc_kms << '\n';
	strm << "p_r\t" << par.pr << '\n';
	strm << "p_z\t" << par.pz << '\n';
	strm << "addJden\t" << par.addJden*intUnits.to_Kpc_kms << '\n';
	strm << "addJvel\t" << par.addJvel*intUnits.to_Kpc_kms << '\n';
}

//added
EXP newExp::newExp(const newExpParam& params) :
    par(params)
{
	if(!(par.norm>0))
		throw std::invalid_argument("newExp: overall normalization must be positive");
	if(!(par.Jr0>0) || !(par.Jz0>0) || !(par.Jphi0>0))
		throw std::invalid_argument("newExp: scale actions must be positive");
	if(par.sigmabirth<=0 || par.sigmabirth>1)
		throw std::invalid_argument("newExp: invalid value for velocity dispersion at birth");
}

EXP double newExp::value(const actions::Actions &J) const
{
	if(J.Jphi<=0) return 0;
    // weighted sum of actions
	double Jsum = fabs(J.Jphi) + par.coefJr * J.Jr + par.coefJz * J.Jz;
	double Jden = sqrt(pow_2(Jsum) + pow_2(par.addJden));
	double Jvel = sqrt(pow_2(Jsum) + pow_2(par.addJvel));
	double Jvelp=pow(par.Jphi0/Jvel,par.power);
	double sigmabirth=par.sigmabirth * sqrt(Jvelp);
	return par.norm / TWO_PI_CUBE / (par.Jr0 * par.Jz0) * pow_2(Jvelp/par.Jphi0) * 
			Jden * tanh(J.Jphi/par.Jcut) * exp(-Jden / par.Jphi0) *
			averageOverAge(Jvelp * J.Jr / par.Jr0, Jvelp * J.Jz / par.Jz0, sigmabirth, par);
}

EXP void newExp::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
//	const units::ExternalUnits& conv);
	strm << "norm\t" << par.norm*intUnits.to_Msun << '\n';
	strm << "Jr0\t" << par.Jr0*intUnits.to_Kpc_kms << '\n';
	strm << "Jz0\t" << par.Jz0*intUnits.to_Kpc_kms << '\n';
	strm << "Jcut\t" << par.Jcut*intUnits.to_Kpc_kms << '\n';
	strm << "Jphi0\t" << par.Jphi0*intUnits.to_Kpc_kms << '\n';
	strm << "Js\t" << par.Js*intUnits.to_Kpc_kms << '\n';
	strm << "power\t" << par.power << '\n';
	strm << "addJden\t" << par.addJden*intUnits.to_Kpc_kms << '\n';
	strm << "addJvel\t" << par.addJvel*intUnits.to_Kpc_kms << '\n';
	strm << "coefJr\t" << par.coefJr << '\n';
	strm << "coefJz\t" << par.coefJz << '\n';
	strm << "beta_r\t" << par.beta_r << '\n';
	strm << "beta_z\t" << par.beta_z << '\n';
	strm << "Tsfr\t" << par.Tsfr << '\n';
	strm << "sigmabirth\t" << par.sigmabirth << '\n';
}

//added
EXP basicExp::basicExp(const basicExpParam& params) :
    par(params)
{
	if(!(par.norm>0))
		throw std::invalid_argument("basicExp: overall normalization must be positive");
	if(!(par.Jr0>0) || !(par.Jz0>0) || !(par.Jphi0>0))
		throw std::invalid_argument("basicExp: scale actions must be positive");
}

EXP double basicExp::value(const actions::Actions &J) const
{
	if(J.Jphi<=0) return 0;
	double Jvel = J.Jphi+par.addJvel;
	double Jden = J.Jphi+par.addJden;
	double Jvelpr=pow(par.Jphi0/Jvel,par.pr);
	double Jvelpz=pow(par.Jphi0/Jvel,par.pz);
	return par.norm / TWO_PI_CUBE / (par.Jr0 * par.Jz0) *
			Jvelpr*Jvelpz/pow_2(par.Jphi0) * Jden * tanh(J.Jphi/par.Jcut) *
			exp(- Jvelpr * J.Jr / par.Jr0 - Jvelpz * J.Jz / par.Jz0 - Jden / par.Jphi0);
}

EXP void basicExp::write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const{
//	const units::ExternalUnits& conv);
	strm << "norm\t" << par.norm*intUnits.to_Msun << '\n';
	strm << "Jr0\t" << par.Jr0*intUnits.to_Kpc_kms << '\n';
	strm << "Jz0\t" << par.Jz0*intUnits.to_Kpc_kms << '\n';
	strm << "Jcut\t" << par.Jcut*intUnits.to_Kpc_kms << '\n';
	strm << "Jphi0\t" << par.Jphi0*intUnits.to_Kpc_kms << '\n';
	strm << "p_r\t" << par.pr << '\n';
	strm << "p_z\t" << par.pz << '\n';
	strm << "addJden\t" << par.addJden*intUnits.to_Kpc_kms << '\n';
	strm << "addJvel\t" << par.addJvel*intUnits.to_Kpc_kms << '\n';
}

}  // namespace df
