/** \file    df_halo.h
    \brief   Distribution functions for the spheroidal component (halo)
    \author  Eugene Vasiliev
    \date    2015-2016
*/
#pragma once
#include <stdexcept>
#include "df_base.h"
#include "utils.h"
#define EXP __declspec(dllexport)

namespace df{

/// \name   Classes for action-based double power-law distribution
/// function (DF), also the same with a core softened per Cole & Binney
///@{

/// Parameters that describe a double power law distribution function.
struct DoublePowerLawParam{
	double
			norm,      ///< normalization factor with the dimension of mass
			J0,        ///< break action (defines the transition between inner and outer regions)
			Jcutoff,   ///< cutoff action (sets exponential suppression at J>Jcutoff, 0 to disable)
			slopeIn,   ///< power-law index for actions below the break action (Gamma)
			slopeOut,  ///< power-law index for actions above the break action (Beta)
			steepness, ///< steepness of the transition between two asymptotic regimes (eta)
			cutoffStrength, ///< steepness of exponential suppression at J>Jcutoff (zeta)
			coefJphiIn,  ///< contribution of azimuthal  action to h(J), controlling anisotropy below J_0 (h_r)
			coefJzIn,  ///< contribution of vertical action to h(J), controlling anisotropy below J_0 (h_z)
			coefJphiOut, ///< contribution of azimuthal action to g(J), controlling anisotropy above J_0 (g_r)
			coefJzOut, ///< contribution of vertical action to g(J), controlling anisotropy above J_0 (g_z)
			rotFrac,   ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
			Jphi0;     ///< controls the steepness of rotation and the size of non-rotating core
	DoublePowerLawParam() :  ///< set default values for all fields (NAN means that it must be set manually)
	    norm(NAN), J0(NAN), Jcutoff(0), slopeIn(NAN), slopeOut(NAN), steepness(1), cutoffStrength(2),
	    coefJphiIn(1), coefJzIn(1), coefJphiOut(1), coefJzOut(1), rotFrac(0), Jphi0(0) {}
};
/// Parameters that describe a double power law distribution function.
struct SoftDoublePowerLawParam{
	double
			norm,      ///< normalization factor with the dimension of mass
			J0,        ///< break action (defines the transition between inner and outer regions)
			Jcutoff,   ///< cutoff action (sets exponential suppression at J>Jcutoff, 0 to disable)
			slopeIn,   ///< power-law index for actions below the break action (Gamma)
			slopeOut,  ///< power-law index for actions above the break action (Beta)
			steepness, ///< steepness of the transition between two asymptotic regimes (eta)
			cutoffStrength, ///< steepness of exponential suppression at J>Jcutoff (zeta)
			coefJphiIn,  ///< contribution of azimuthal action to h(J), controlling anisotropy below J_0 (h_r)
			coefJzIn,  ///< contribution of vertical action to h(J), controlling anisotropy below J_0 (h_z)
			coefJphiOut, ///< contribution of azimuthal action to g(J), controlling anisotropy above J_0 (g_r)
			coefJzOut, ///< contribution of vertical action to g(J), controlling anisotropy above J_0 (g_z)
			rotFrac,   ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
			Jphi0,     ///< controls the steepness of rotation and the size of non-rotating core
			h0;        ///< controls Cole&Binney core size
	SoftDoublePowerLawParam() :  ///< set default values for all fields (NAN means that it must be set manually)
	    norm(NAN), J0(NAN), Jcutoff(0), slopeIn(NAN), slopeOut(NAN), steepness(1), cutoffStrength(2),
	    coefJphiIn(1), coefJzIn(1), coefJphiOut(1), coefJzOut(1), rotFrac(0), Jphi0(0), h0(0) {}
};

/** General double power-law model.
    The distribution function is given by
    \f$  f(J) = norm / (2\pi J_0)^3
         (1 + (J_0 /h(J))^\eta )^{\Gamma / \eta}
         (1 + (g(J)/ J_0)^\eta )^{-B / \eta }
         \exp[ - (g(J) / J_{cutoff})^\zeta ] \f$,  where
    \f$  g(J) = (3-g_phi-g_z) J_r + 0.7(g_z J_z + g_\phi |J_\phi|)  \f$,
    \f$  h(J) = (3-h_phi-h_z) J_r + 0.7(h_z J_z + h_\phi |J_\phi|)  \f$.
    Gamma is the power-law slope of DF at small J (slopeIn), and Beta -- at large J (slopeOut),
    the transition occurs around J=J0, and its steepness is adjusted by the parameter eta.
    h_r, h_z and h_phi control the anisotropy of the DF at small J (their sum is always taken
    to be unity, so that there are two free parameters -- coefJphiIn =
    h_phi, coefJzIn = h_z),
    and g_z, g_phi do the same for large J (coefJphiOut = g_phi, coefJzOut = g_z).
    Jcutoff is the threshold of an optional exponential suppression, and zeta measures its strength.
*/
class EXP DoublePowerLaw: public BaseDistributionFunction{
	const DoublePowerLawParam par;  ///< parameters of DF
	public:
    /** Create an instance of double-power-law distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
		DoublePowerLaw(const DoublePowerLawParam &params);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
		virtual double value(const actions::Actions &J) const;
		void write_params(std::ofstream&,const units::InternalUnits&) const;
};

//Now code that generates DFs as above but modified per Cole & Binney
//to produce a finite central action-space density. The extent in J of
//the region of ~const action-space density is set by h0. The machine,  not the user
//sets beta

class EXP SoftDoublePowerLaw: public BaseDistributionFunction{
	const SoftDoublePowerLawParam par;  ///< parameters of DF
	private:
		double beta;
		double g(const double hJ) const;
		void set_beta(void);
		SoftDoublePowerLaw(const SoftDoublePowerLawParam &params_,double beta_) : par(params_), beta(beta_){}
	public:
    /** Create an instance of ColeBinney double-power-law distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
		SoftDoublePowerLaw(const SoftDoublePowerLawParam&);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
		virtual double value(const actions::Actions &J) const;
		void write_params(std::ofstream&,const units::InternalUnits&) const;
		double diff(const double) const;
		double intDiff(const double) const;
};

//We start with two helper classes used to determine Cole&Binney
//parameter beta

class ColeBinneyDiff: public math::IFunctionNoDeriv{
	public:
		ColeBinneyDiff(SoftDoublePowerLaw &Sft_): Sft(Sft_) {}
		virtual double value(const double hJ) const{
			return Sft.diff(hJ);
		}
		const SoftDoublePowerLaw& Sft;
};
class ColeBinneyInt: public math::IFunctionNoDeriv{
	public:
		ColeBinneyInt(SoftDoublePowerLaw &Sft_): Sft(Sft_) {}
		virtual double value(const double x) const{
			return Sft.intDiff(x);
		}
		const SoftDoublePowerLaw& Sft;
};

}  // namespace df
