#pragma once
#include "actions_newgenfnc.h"
#include "actions_newisochrone.h"
#include "actions_focal_distance_finder.h"
#include "potential_base.h"
#include "potential_utils.h"
#include "math_core.h"
#include "math_specfunc.h"
#include "math_fit.h"
#include "math_spline.h"
#include "math_fourier.h"
#include "coord.h"
#include "particles_base.h"
#include <iostream>
#include <map>
#include <complex>
#include <cassert>

namespace actions{

/* Class to hold the Fourier decomposition of the residual Hamiltonian
 * after torus fitting
*/
class EXP PerturbingHamiltonian{
	private:
		GenFncIndices indices;
		std::vector<std::complex<double> > values;
	public:
		PerturbingHamiltonian(GenFncIndices& _indices,
				      std::vector<std::complex<double> >& _values) :
		    indices(_indices), values(_values){}
		GenFncIndex index(const int i) const {
			return indices[i];
		}
		std::complex<double> value(const int i) const {
			return values[i];
		}
		int numTerms(){
			return indices.size();
		}
		PerturbingHamiltonian& operator *= (const double a);
		PerturbingHamiltonian  operator * (const double a);
		PerturbingHamiltonian& operator += (const PerturbingHamiltonian&);
		PerturbingHamiltonian  operator + (const PerturbingHamiltonian&);
};
/*
 * Base of all torus classes
 */
class EXP BaseTorus {
	public:
		Actions J;
		Frequencies freqs;
		GenFnc GF;
		double E;
		coord::UVSph cs;
		Isochrone TM;
		/* Creator called by TorusGenerator rather than users */
		BaseTorus(const Actions& _J, const Frequencies& _freqs, const GenFnc& _GF,
			  const Isochrone& _TM, const coord::UVSph& _cs, double _E) :
		    J(_J), freqs(_freqs), GF(_GF), TM(_TM), cs(_cs), E(_E) {
			printf("BaseTorus created: %d params\n",_GF.numParams());
		}
		void printGF(void){
			GF.print();
		}
		coord::PosMomCyl from_toy (const Angles&) const;
		coord::PosMomCyl from_true(const Angles&) const;
		coord::PosMomCyl from_aaT(const ActionAngles&) const;//toy inputs
		coord::PosMomCyl from_aa(const Angles&) const;//true input
		coord::PosCyl PosDerivs(const Angles&, actions::DerivAng<coord::Cyl>&,
					double* det=NULL) const;/*dR/dtheta etx*/
		/* compute the surface of section z=0 pz>0 */
		void zSoS(std::vector<double>& R,std::vector<double>& vR,int n,
			  double& Rmin, double& Rmax, double &Vmax) const;
		Frequencies Omega(void) const{
			return freqs;
		}
		/* Compute the orbit from given angles using computed frequencies */
		std::vector<std::pair<coord::PosVelCyl,double> > orbit(const Angles& theta0, double dt, double T) const;
		/* Does the torus pass through given point? If so for what thetas? */
		bool containsPoint(const coord::PosCyl& pt, std::vector<Angles>& As,
				   std::vector<coord::VelCyl>& Vs,
				   std::vector<double>& Jacobs, const double& tol) const;
		/* returns density contributed at location */
		double density(const coord::PosCyl&) const;
		void write(FILE*) const;
		void read(FILE*);
};
/* The class actually used for basic tori. It's set up for interpolation */
class EXP Torus : public BaseTorus {
	public:
		Torus(const Actions& _J, const Frequencies& _freqs, const GenFnc& _GF,
		      const Isochrone& _TM, const coord::UVSph& _cs, double _E) :
		    BaseTorus(_J, _freqs, _GF, _TM, _cs, _E) {}
		Torus(BaseTorus& T) : BaseTorus(T) {}
		Torus&	operator *= (const double);
		Torus&	operator += (const Torus&);
		const Torus operator * (const double);
		const Torus operator + (const Torus&);
};
inline Torus interp(const double x, Torus T0, Torus T1){
	Torus T2(T0*x);
	T2+=T1*(1-x);
	return T2;
}
/* The class of tori that include perturbing Hamiltonians. Can be
 * interpolated */
class EXP eTorus : public BaseTorus {
	private:
		PerturbingHamiltonian pH;
	public:
		eTorus(const Actions& _J, const Frequencies& _freqs, const GenFnc& _GF,
		       const Isochrone& _TM, const coord::UVSph& _cs, double _E,
		       const PerturbingHamiltonian& _pH) :
		    BaseTorus(_J, _freqs, _GF, _TM, _cs, _E), pH(_pH) {
			printf("eTorus created: %d terms in pH\n", pH.numTerms());
		}
		eTorus(const BaseTorus& T,const PerturbingHamiltonian& _pH) :
		    BaseTorus(T), pH(_pH) {
			printf("eTorus created: %d terms in pH\n", pH.numTerms());
		}
		eTorus&	operator *= (const double);
		eTorus&	operator += (const eTorus&);
		const eTorus operator * (const double);
		const eTorus operator + (const eTorus&);
};
/*
 * Class for generators of tori.
*/
class EXP TorusGenerator : public math::IFunctionNdimDeriv{
	private:
		const potential::BasePotential& pot;
		double Js_iso, b_iso, Delta, tol, Rsh, freqScale, invPhi0;
		Isochrone TM;
		double NANfrac;
		GenFncFit* GFF;
		coord::UVSph cs;
		math::LinearInterpolator2d interpD; //for Delta values
		math::LinearInterpolator2d interpR; //for Rshell values
#ifdef TEST
		/* Test_it compares analytic and numerical derivatives */
		void test_it(actions::Actions& J, std::vector<double>&);
#endif
		/* Hamiltonian computes residual H for Phi+ePhi where
		 * ePhi is an additional potential not used in fitting */
		double Hamilton(const BaseTorus&, const potential::BasePotential*, const Angles&);
		PerturbingHamiltonian get_pH (const BaseTorus&,
			int, bool, const potential::BasePotential*);
		void setConsts(actions::Actions, double&);
		double fitAngleMap(double*, double&, Frequencies&, GenFncDerivs&) const;
		int tmax;// Max number of terms retained in residual H
	public:
		/* Creator of tori in given potential. GF deemed ok if
		 * dispersion in H < tol*freqScale*Jtotal */
		TorusGenerator(const potential::BasePotential& _pot,const double _tol = 1e-9);
		BaseTorus fitBaseTorus(const Actions&);// The actual workhorse
		Torus fitTorus(const Actions&);
		eTorus fiteTorus(const Actions&, const potential::BasePotential* _addPhi=NULL);
		double computeHamiltonianDisp(const std::vector<double> &params, double& Hbar);
		virtual unsigned int numVars() const {
			return GFF->numParams();
		}
		virtual unsigned int numValues() const {
			return GFF->numPoints();
		}
		virtual void evalDeriv(const double params[],
				       double* deltaHvalues=NULL, double* dHdParams=NULL) const;
		double computeHamiltonianAtPoint(const double params[], const unsigned int indPoint,
			Actions *dHdJ=NULL, double *derivGenFnc=NULL) const;
		double getDelta(double&,double&);
		double getDelta(Actions&);
		double getRsh(Actions&);
		void test_it(const Actions&, std::vector<double>&);
		std::vector<Torus> constE(const double Jrmin, const Actions& Jstart, const int Nstep);
};

}//namespace actions