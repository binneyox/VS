#pragma once
#include "actions_base.h"
#include "actions_staeckel.h"
#include "actions_newgenfnc.h"
#include "actions_newisochrone.h"
#include "actions_harmonicoscillator.h"
#include "actions_focal_distance_finder.h"
#include "actions_spherical.h"
#include "potential_base.h"
#include "potential_utils.h"
#include "math_core.h"
#include "math_specfunc.h"
#include "math_fit.h"
#include "math_linalg.h"
#include "math_fourier.h"
#include "coord.h"
#include "particles_base.h"
#include <iostream>
#include <map>
#include <complex>
#include <cassert>


namespace actions {
	class scaler {
	public:
		double r0;
		scaler() :r0(1.0) {}
		scaler(double _r0) :r0(_r0) {}
		double s(double r, double* dsdR = NULL, double* d2sdR2 = NULL) const {
			//double sq = sqrt(R * R + 1);
			double tnhx = tanh(r / r0);
			double sechx2 = 1 - tnhx * tnhx;
			double s = .5 * (tnhx + 1);
			if (dsdR) *dsdR = .5 * sechx2 / (r0);
			if (d2sdR2) *d2sdR2 = -sechx2 * tnhx / (r0 * r0);
			return s;
		}
		double R(double s, double* dRds = NULL, double* d2Rds2 = NULL) const {
			double r1 = r0 * atanh(2 * s - 1);
			double inv = 1. / (1 - pow_2(2 * s - 1));
			if (dRds) *dRds = 2 * r0 * inv;
			if (d2Rds2) *d2Rds2 = 8 * r0 * (2 * s - 1) * pow_2(inv);
			return r1;
		}

	};
	//integrates pz^n between endpoints
	class pzf :public math::IFunction {
	private:
		const double E;
		const potential::BasePotential& pot;
	public:
		pzf(double _E, const potential::BasePotential& _pot) :E(_E), pot(_pot) {};
		virtual void evalDeriv(const double z,
			double* value = 0, double* deriv = 0, double* deriv2 = 0) const {
			coord::PosCar X(0, 0, z);
			double Phi;
			coord::GradCar dx;
			pot.eval(X, &Phi, &dx);
			//std::cout << "Value:" << x << " " << pot.value(X) << "\n";
			double p2 = 2 * (E - Phi);
			if (p2 <= 0) {
				if (value)*value = 0;
			}
			else {
				double p = sqrt(p2);
				if (value)*value = p;
				if (deriv)*deriv = -dx.dz / (p);
			}
		};
		virtual unsigned int numDerivs() const { return 1; }
	};
	//finds the z1 s.t. Iz=int pzdz between 0 and z1, pz is the z momentum at that energy.
	class Ifind :public math::IFunction {
	private:
		const potential::BasePotential& pot;
		const double I;
		const double E;
		pzf pzfunc;
	public:
		Ifind(double _Iz, const potential::BasePotential& _pot, double _E)
			:I(_Iz), pot(_pot), E(_E),pzfunc(_E,_pot){
		};
		virtual void evalDeriv(const double z,
			double* value = 0, double* deriv = 0, double* deriv2 = 0)  const {
			if (value)*value = math::integrateGL(math::ScaledIntegrand<math::ScalingCub>
				(math::ScalingCub(0, z), pzfunc), 0, 1, math::MAX_GL_ORDER) - I;
			if (deriv) {
				double pz;
				pzfunc.evalDeriv(z, &pz);
				*deriv = pz;
			}
		}
		virtual unsigned int numDerivs() const { return 1; }
	};
	/*
	 Maps sphere to prolate ellipsoid extended along z
	*/
	class EXP PointTrans {
	private:
	public:
		coord::UVSph cs;
		bool map;
		std::vector<double> paramsF;
		std::vector<double> paramsFr;
		scaler sc;
		scaler scz;
		int N;
		int Nr;
		//scale factor r'=k*r+fourier series in s
		PointTrans();
		PointTrans(double _D);
		PointTrans(coord::UVSph _cs);
		PointTrans(double _D, const std::vector<double>& _p) :
			map(true), cs(_D), paramsF(_p), Nr(0), paramsFr({}), N(_p.size()) {
		}
		PointTrans(coord::UVSph _cs, const std::vector<double>& _p) :
			map(true), cs(_cs), paramsF(_p), Nr(0), paramsFr({}), N(_p.size()) {
		}
		PointTrans(double _D, scaler _sc, std::vector<double>& _p, std::vector<double> _pr) :
			map(true), cs(_D), paramsF(_p), paramsFr(_pr), Nr(_pr.size()), N(_p.size()), sc(_sc) {
		}
		PointTrans(coord::UVSph _cs, scaler _sc, std::vector<double>& _p, std::vector<double> _pr) :
			map(true), cs(_cs), paramsF(_p), paramsFr(_pr), Nr(_pr.size()), N(_p.size()), sc(_sc) {
		}
		PointTrans(double _D, scaler _sc, scaler _scz, std::vector<double>& _p, std::vector<double> _pr) :
			map(true), cs(_D), paramsF(_p), paramsFr(_pr), Nr(_pr.size()), N(_p.size()), sc(_sc), scz(_scz) {
		}
		PointTrans(coord::UVSph _cs, scaler _sc, scaler _scz, std::vector<double>& _p, std::vector<double> _pr) :
			map(true), cs(_cs), paramsF(_p), paramsFr(_pr), Nr(_pr.size()), N(_p.size()), sc(_sc), scz(_scz) {
		}
		double t2v(const double theta, double* dvdt = NULL, double* d2vdt2 = NULL) const;
		double v2t(const double, double* = NULL) const;
		double r2rn(const double rn, double* drndr = NULL, double* d2rndr2 = NULL) const;
		double rn2r(const double, double* = NULL) const;
		double zntov(const double zn, double* dvdt = NULL, double* d2vdt2 = NULL) const;
		double vtozn(const double v, double* dzdv = NULL) const;
		coord::PosMomSph Cyl2Sph(const coord::PosMomCyl Rz) const;
		coord::PosMomCyl Sph2Cyl(const coord::PosMomSph rp,
			coord::PosMomCyl* dRzdDelta = NULL) const;
		coord::PosMomCyl Sph2Cyl(const coord::PosMomSph rp,
			math::Matrix<double>& dRzdrv, coord::PosMomCyl* dRzdDelta = NULL) const;
		coord::PosMomCyl R2Rn(const coord::PosMomCyl Rn, math::Matrix<double>& dRnzndRz) const;
		coord::PosMomCyl R2Rn(const coord::PosMomCyl Rn) const;
		coord::PosMomCyl Rn2R(const coord::PosMomCyl Rz) const;
	};
	EXP PointTrans interpPointTrans(double x, const PointTrans&, const PointTrans&);
	/*
	 Combination of a PointTrans and an Isochrone aa map
	*/
	class EXP ToyMap {
	private:
	public:
		Isochrone Is;
		HarmonicOscilattor HOs;
		PointTrans PT;
		bool useIso;
		ToyMap() {}
		ToyMap(Isochrone _Is, PointTrans _PT) : Is(_Is), PT(_PT), useIso(true) {}
		ToyMap(HarmonicOscilattor _HOs, PointTrans _PT) : HOs(_HOs), PT(_PT), useIso(false) {}
		ToyMap(const ToyMap& TM) : PT(TM.PT), Is(TM.Is), useIso(TM.useIso), HOs(TM.HOs) {
		}
		Actions pq2J(const coord::PosMomCyl Rzp) const {
			if (useIso)return Is.pq2J(PT.Cyl2Sph(Rzp));
			else return HOs.pq2J(PT.Rn2R(Rzp));
		}
		ActionAngles pq2aa(const coord::PosMomCyl& Rz) const {
			if (useIso) {
				coord::PosMomSph rp(PT.Cyl2Sph(Rz));
				if (Is.H(rp) < 0) return Is.pq2aa(rp);
				printf("ToyMap::pq2aa: Iso H = %g >=0!",Is.H(rp));
				exit(0);
			} else {
				coord::PosMomCyl rp(PT.Rn2R(Rz));
				return HOs.pq2aa(rp);
			}
		}
		coord::PosMomCyl from_aaT(const ActionAngles& aaT) const {
			if (useIso) return PT.Sph2Cyl(Is.aa2pq(aaT));
			else return PT.R2Rn(HOs.aa2pq(aaT));
		}
		coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta) const {
			return from_aaT(ActionAngles(J, theta));
		}
		coord::PosMomCyl from_aaT(const ActionAngles& aaT, coord::PosMomCyl* dRzdPs) const;

		coord::PosMomCyl from_aaT(const Actions& J, Angles& theta, coord::PosMomCyl* dRzdPs) const {
			return from_aaT(ActionAngles(J, theta), dRzdPs);
		}
		coord::PosMomCyl from_aaT(const ActionAngles& aaT, DerivAct<coord::Cyl>& dRzdJ) const;

		coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta,
			DerivAct<coord::Cyl>& dRzdJ) const {
			return from_aaT(ActionAngles(J, theta), dRzdJ);
		}
		coord::PosMomCyl from_aaT(const ActionAngles& aaT, DerivAng<coord::Cyl>& dRzdT) const;
		coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta,
			DerivAng<coord::Cyl>& dRzdt) const {
			return from_aaT(ActionAngles(J, theta), dRzdt);
		}

	};

	EXP ToyMap interpToyMap(double, const ToyMap&, const ToyMap&);

	/*
	  Map of equatorial plane of (r,theta,phi) to system in (x,y) plane
	  extended along y
	*/
	class PsiCoord {
	private:
		const double a, b;
	public:
		PsiCoord(const double _a, const double _b) : a(_a), b(_b) {}
		double Psi(const double phi, double& dpsidphi) const {
			dpsidphi = 1 + 2 * a * cos(2 * phi) + 4 * b * cos(4 * phi);
			return phi + a * sin(2 * phi) + b * sin(4 * phi);
		}
		double Phi(const double psi, double& dpsidphi) const {
			double phi = psi, D = Psi(phi, dpsidphi) - psi;
			while (fabs(D) > 1e-5) {
				phi -= D / dpsidphi;
				D = Psi(phi, dpsidphi) - psi;
			}
			return phi;
		}
	};
	class EXP xyPointTrans {
		const double Delta, Delta2;
		PsiCoord PC;
	public:
		xyPointTrans(double D, double a, double b) : Delta(D), Delta2(D* D), PC(a, b) {}
		coord::PosMomCar rp2xp(const coord::PosMomSph rp) const;
		coord::PosMomSph xp2rp(const coord::PosMomCar xp) const;
	};
	class EXP xyToyMap {
	private:
	public:
		Isochrone Is;
		xyPointTrans PT;
		xyToyMap(double _D, double _aPT, double _bPT, double _Js, double _bIso) :
			PT(_D, _aPT, _bPT), Is(_Js, _bIso) {
		}
		double H(const coord::PosMomCar xp) const {
			coord::PosMomSph rp(PT.xp2rp(xp));
			return Is.H(rp);
		}
		Actions pq2J(const coord::PosMomCar xp) const {
			coord::PosMomSph rp(PT.xp2rp(xp));
			if (Is.H(rp) < 0) return Is.pq2J(rp);
			printf("ToyMap::pq2J: H>=0!");
			exit(0);
		}
		ActionAngles pq2aa(const coord::PosMomCar& xp) const {
			coord::PosMomSph rp(PT.xp2rp(xp));
			if (Is.H(rp) < 0) return Is.pq2aa(rp);
			printf("ToyMap::pq2aa: H>=0!");
			exit(0);
		}
		coord::PosMomCar aa2pq(ActionAngles aaT) const {
			return PT.rp2xp(Is.aa2pq(aaT));
		}
		coord::PosMomCar aa2pq(Actions J, Angles thetas) const {
			return aa2pq(ActionAngles(J, thetas));
		}
	};
	/* Class to hold the Fourier decomposition of the residual Hamiltonian
	 * after torus fitting
	*/
	class EXP PerturbingHamiltonian {
	private:
	public:
		GenFncIndices indices;
		std::vector<std::complex<double> > values;
		PerturbingHamiltonian(const GenFncIndices& _indices,
			const std::vector<std::complex<double> >& _values) :
			indices(_indices), values(_values) {
		}
		GenFncIndex index(const int i) const {
			return indices[i];
		}
		std::complex<double> value(const int i) const {
			return values[i];
		}
		int numTerms() {
			return indices.size();
		}
		std::vector<std::complex<double> > get_hn(const GenFncIndex&, std::vector<float>&) const;
	};
	EXP PerturbingHamiltonian interpPerturbingHamiltonian(const double x,
		const PerturbingHamiltonian& H0, const PerturbingHamiltonian& H1);

	/*
	 * Base of all torus classes
	 */
	class EXP Torus {
	public:
		Actions J;
		Frequencies freqs;
		GenFnc GF;
		ToyMap TM;
		double E;
		bool negJr, negJz;
		Torus() {}
		/* Creator called by TorusGenerator rather than users */
		Torus(const Actions& _J, const Frequencies& _freqs, const GenFnc& _GF,
			const ToyMap _TM, double _E, bool _negJr, bool _negJz) :
		    J(_J), freqs(_freqs), GF(_GF), TM(_TM), E(_E),
		    negJr(_negJr), negJz(_negJz) {
		}
		Torus(const Torus& T) :
		    J(T.J), freqs(T.freqs), GF(T.GF), TM(T.TM), E(T.E),
		    negJr(T.negJr), negJz(T.negJz) {
		}
		void printGF(void) {
			GF.print();
		}
		coord::PosMomCyl from_toy(const Angles&) const;
		coord::PosMomCyl from_true(const Angles&) const;
		coord::PosCyl new_PosDerivJ(const Angles&, actions::DerivAct<coord::Cyl>&) const;
		coord::PosCyl new_PosDerivs(const Angles&, actions::DerivAng<coord::Cyl>&,
			double* det = NULL) const;/*dR/dtheta etx*/
		/* compute the surface of section z=0 pz>0 */
		void zSoS(std::vector<double>& R, std::vector<double>& vR, const int N,
			double& Rmin, double& Rmax, double& Vmax, const double z0 = 0) const;
		/* compute the surface of section R=Rbar pR>0 */
		void rSoS(std::vector<double>& z, std::vector<double>& vz, const double Rbar, const int N,
			double& zmax, double& Vmax, double delta = 1e-3) const;
		void SoSthetaz(std::vector<double>& X, std::vector<double>& pX, const double thetaz,
			const int N, double& Xmax, double& pXmax) const;
		void SoSthetar(std::vector<double>& X, std::vector<double>& pX, const double thetar,
			const int N, double& Xmax, double& pXmax) const;
		Frequencies Omega(void) const {
			return freqs;
		}
		/* Compute the orbit from given angles using computed frequencies */
		std::vector<std::pair<coord::PosVelCyl, double> > orbit(const Angles& theta0, double dt, double T) const;
		/* Does the torus pass through given point? If so for
		 * what true thetas? */
		bool containsPoint(const coord::PosCyl& pt, std::vector<Angles>& As,
			std::vector<coord::VelCyl>& Vs,
			std::vector<double>& Jacobs,
			std::vector<actions::DerivAngCyl>* dA = NULL,
			const double tol = 1e-5) const;
		/* returns density contributed at location */
		double density(const coord::PosCyl&) const;
		void write(FILE*) const;
		void read(FILE*);
	};

	//interpolate between 2 tori
	EXP Torus interpTorus(const double x, const Torus& T0, const Torus& T1);

	//interpolae on an indexed array of tori 
	EXP Torus interpTorus(const double x, std::vector<double>&, std::vector<Torus>&);

	class EXP TorusGrid1 {
	private:
		const std::vector<double>& xs;
		const::std::vector<Torus>& Ts;
		const int nx;
		int botX(const double) const;
	public:
		TorusGrid1(std::vector<double>& _xs,
			std::vector<actions::Torus>& _Ts) :
			xs(_xs), Ts(_Ts), nx(xs.size()) {
		}
		Torus T(const double x) const;
	};
	class EXP TorusGrid3 {
	private:
		const std::vector<double>& xs, ys, zs;
		const::std::vector<Torus>& Ts;
		const int nx, ny, nz;
		int botX(const double) const;
		int botY(const double) const;
		int botZ(const double) const;
		Torus Tn(int ix, int iy, int iz) const {
			return Ts[(ix * ny + iy) * ny + iz];
		}
	public:
		TorusGrid3(std::vector<double>& _xs, std::vector<double>& _ys,
			std::vector<double>& _zs, std::vector<actions::Torus>& _Ts) :
			xs(_xs), ys(_ys), zs(_zs), Ts(_Ts), nx(xs.size()), ny(ys.size()), nz(zs.size()) {
		}
		Torus T(const double x, const double y, const double z) const;
		Torus T(const Actions J) const {
			return T(J.Jr, J.Jz, J.Jphi);
		}
	};
	/* The class of tori that include perturbing Hamiltonians. */
	class EXP eTorus : public Torus {
	private:
	public:
		PerturbingHamiltonian pH;
		eTorus(const Actions& _J, const Frequencies& _freqs, const GenFnc& _GF,
			const ToyMap _TM, double _E, bool _negJr, bool _negJz,
			const PerturbingHamiltonian& _pH) :
			Torus(_J, _freqs, _GF, _TM, _E, negJr, negJz), pH(_pH) {
			//printf("eTorus created: %d terms in pH\n", pH.numTerms());
		}
		eTorus(const Torus& T, const PerturbingHamiltonian& _pH) :
			Torus(T), pH(_pH) {
			//printf("eTorus at E = %f created: %d terms in pH\n",T.E, pH.numTerms());
		}
		PerturbingHamiltonian Hns() const {
			return pH;
		}
		std::vector<std::complex<double> > get_hn(const GenFncIndex& Indx, std::vector<float>& multiples) const {
			return pH.get_hn(Indx, multiples);
		}
	};

	//interpolate between 2 tori
	EXP eTorus interpeTorus(const double x, const eTorus& T0, const eTorus& T1);

	//interpolate on an indexed array of tori 
	EXP eTorus interpeTorus(const double x, std::vector<double>&, std::vector<eTorus>&);

	/*
	 * Class for fitting torus to an orbit
	*/
	class EXP TMfitter : public math::IFunctionNoDeriv {
	private:
		double pphi, Delta2, xmin, ymax, xbar, Frat, aPT, bPT;
		std::vector<std::pair<coord::PosMomCyl, double> >& traj;
	public:
		TMfitter(const potential::BasePotential&,
			std::vector<std::pair<coord::PosMomCyl, double> >&, double);
		std::vector<double> fitTM() const;
		virtual double value(double) const;
		virtual unsigned int numVars() const { return 4; }
		virtual unsigned int numValues() const { return 1; }
	};
	/*
	 * Class for generators of tori.
	*/
	class EXP TorusGenerator {
	private:
		const potential::BasePotential& pot;
		const double defaultTol, invPhi0;
		math::LinearInterpolator2d interpD; //for Delta values
		math::LinearInterpolator2d interpR; //for Rshell(L,Xi) values
		math::LinearInterpolator2d interpV; //for Vshell(L,Xi) values
		math::LinearInterpolator2d interpRE; //for Rshell(E,Xi) values
		math::LinearInterpolator interpJz;//Jr(Jr) for box loop orbit transition
		ActionFinderSpherical afs;
		std::string logfname;
#ifdef TEST
		/* Test_it compares analytic and numerical derivatives */
		void test_it(actions::Actions& J, std::vector<double>&);
#endif
		/* Hamilton computes residual H for Phi+ePhi where
		 * ePhi is an optional additional potential not used in fitting */
		double Hamilton(const Torus&, const potential::BasePotential*, const Angles&);
		/* PerturbingHamiltonian computes Fourier decomposition
		 * of the residual H */
		PerturbingHamiltonian get_pH(const Torus&,
			int nf, bool ifp, const potential::BasePotential*);
		//		void setConsts(actions::Actions, double, double&, double&, double&, Isochrone&, coord::UVSph&) const;
				//ToyMap chooseTM(actions::Actions, double&, double&, double&) const;
		int tmax;// Max number of terms retained in residual H
	public:
		/* Creator of tori in given potential. GF deemed ok if
		 * dispersion in H < tol*freqScale*Jtotal */
		TorusGenerator(const potential::BasePotential& _pot,
			const double _tol = 1e-9, std::string _logfname = "TG.log");
		ToyMap chooseTM(GenFncFitSeries&, std::vector<double>&, const Actions&,
			double&, double&, double&, bool = true, FILE* logfile = NULL) const;
		Torus fitTorus(const Actions& J, const double tighten = 1, const bool useIso = true) const;
		/* Fit a torus with all Sn=0 */
		Torus fitBaseTorus(const Actions&, const double tighten = 1) const;
		/* Build a torus with all Sn=0 around the given ToyMap */
		Torus giveBaseTorus(const Actions&, const ToyMap&) const;
		/* fitFrom uses given TM & varies Sn starting from given values
		*/
		Torus fitFrom(const Actions&, const Torus&, const double tighten = 1) const;
		eTorus fiteTorus(const Actions&, const potential::BasePotential* _addPhi = NULL);
		eTorus fiteTorus(const Actions&, const double, const potential::BasePotential* _addPhi = NULL);
		double getDelta(double&, double&);
		double getDelta(Actions&);
		double getRsh(Actions&);
		void test_it(const Actions&, std::vector<double>&);
		std::vector<Torus> constE(const double Jrmin, const Actions& Jstart, const int Nstep);
		void old_getHn(const Torus&, int);
	};

	class EXP ActionFinderTG : public BaseActionFinder {
	private:
		const potential::PtrPotential pot;
		const TorusGenerator& TG;
		const ActionFinderAxisymFudge* AF;

	public:
		ActionFinderTG(const potential::PtrPotential& _pot,
			const TorusGenerator& _TG,
			const ActionFinderAxisymFudge* _AF) :
			pot(_pot), TG(_TG), AF(_AF) {
		}
		ActionFinderTG(const potential::PtrPotential& _pot,
			const TorusGenerator& _TG) :
			pot(_pot), TG(_TG), AF(new ActionFinderAxisymFudge(_pot)) {
		}
		//virtual Actions actions(const coord::PosMomCyl& point) const;
		virtual ActionAngles actionAngles(const coord::PosVelCyl& point,
			Frequencies* freq = NULL) const;
		virtual ActionAngles actionAngles(const coord::PosVelCar(point)) const {
			return actionAngles(coord::toPosVelCyl(point));
		}
		virtual ActionAngles actionAnglesTorus(const coord::PosVelCyl& point,
			Torus& T) const;
		virtual ActionAngles actionAngles(const coord::PosCar& x, const coord::VelCar& v)  const {
			return actionAngles(coord::PosVelCar(x, v));
		}
		virtual Actions actions(const coord::PosVelCyl& point) const {
			return Actions(actionAngles(point));
		}
	};
	EXP void getGridBoxLoop(const potential::BasePotential& pot,
				std::vector<double>& gridE,
				std::vector<double>& gridJr,
				std::vector<double>& gridJz);
	EXP void mapJcrit(const potential::BasePotential& pot,
			  math::CubicSpline& Jrcrit, math::CubicSpline& Jzcrit);
}//namespace actions