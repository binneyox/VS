#include "actions_newtorus.h"
#include "orbit.h"
#include <omp.h>

/// accuracy parameter determining the spacing of the interpolation grid along the energy axis
static const double ACCURACY_INTERP2 = 1e-6;

namespace actions {

	namespace {//internal
		double NANfrac, NANbar;

		struct PVUVSph {
			double u, v, phi, pu, pv, pphi;
		};
		struct DerivActUVSph {
			PVUVSph dbyJr, dbyJz, dbyJphi;
		};
		//Two functions used in containsPoint
		double angleDiff(Angles A1, Angles A2) {

			double dr = A1.thetar - A2.thetar;
			while (dr > M_PI) dr -= 2 * M_PI;
			while (dr < -M_PI) dr += 2 * M_PI;
			double dz = A1.thetaz - A2.thetaz;
			while (dz > M_PI) dz -= 2 * M_PI;
			while (dz < -M_PI) dz += 2 * M_PI;
			double dp = A1.thetaphi - A2.thetaphi;
			while (dp > M_PI) dp -= 2 * M_PI;
			while (dp < -M_PI) dp += 2 * M_PI;
			return sqrt(dr * dr + dz * dz + dp * dp);
		}

		bool is_new(Angles A1, std::vector<Angles> As) {
			bool ok = true; const double tiny = 1e-5;
			for (int i = 0; i < As.size(); i++) {
				double diff = angleDiff(A1, As[i]);
				ok = ok && diff > tiny;
			}
			return ok;
		}
		class locFinder : public math::IFunctionNdimDeriv {
		private:
			const Torus T;
			const coord::PosCyl P0;
		public:
			locFinder(const Torus _T, const coord::PosCyl& _P0) : T(_T), P0(_P0) {}
			virtual void evalDeriv(const double params[], double* values, double* derivs = NULL) const {
				Angles A1(params[0], params[1], 0);
				coord::PosCyl P1;  DerivAng<coord::Cyl> dA;
				if (derivs) {
					P1 = (T.new_PosDerivs(A1, dA));
					derivs[0] = dA.dbythetar.R; derivs[1] = dA.dbythetaz.R;
					derivs[2] = dA.dbythetar.z; derivs[3] = dA.dbythetaz.z;
				}
				if (values && !derivs)
					P1 = T.from_toy(A1);
				if (values) {
					values[0] = P1.R - P0.R; values[1] = P1.z - P0.z;
					double dist = sqrt(pow_2(values[0]) + pow_2(values[1]));
				}
			}
			virtual unsigned int numVars() const {
				return 2;
			}
			virtual unsigned int numValues() const {
				return 2;
			}
		};

		double insertLine(int& ntop, const int tmax, double s, GenFncIndex I,
			std::vector<double>& Hmods, GenFncIndices& Hindices) {
			if (ntop == 0) {
				Hmods.push_back(s);
				Hindices.push_back(I);
				ntop++; return Hmods[0];
			}
			int l = 0;
			while (l < ntop && s <= Hmods[l]) l++;
			if (l == ntop) {//line isn't stronger than any previous line
				if (ntop >= tmax) {
					return Hmods[ntop - 1];//no room for this line
				}
				else {//add line to end of list
					Hmods.push_back(s);
					Hindices.push_back(I);
					ntop++;
					return Hmods.back();
				}
			}
			else {//we should insert current line
				if (ntop < tmax) {//move existing terms down
					Hmods.push_back(Hmods[ntop - 1]);
					Hindices.push_back(Hindices[ntop - 1]);
				}
				for (int m = ntop - 1; m > l; m--) {
					Hmods[m] = Hmods[m - 1];
					Hindices[m] = Hindices[m - 1];
				}//we've created a space at l, so fill it
				Hmods[l] = s;
				Hindices[l] = I;
				if (ntop < tmax) ntop++;//we've added rather than replaced a line
				return Hmods.back();//return weakest retained line
			}
		}
		/* compute the best focal distance at a 2d grid in L, Xi=Jz/L
		 * on input gridL, which is a grid in Jcirc, times gridXi is a uniform grid on (0,1)
		*/
		void createGridFocalDistance(const potential::BasePotential& pot,
			std::vector<double>& gridL,
			std::vector<double>& gridXi,
			std::vector<double>& gridE,
			math::Matrix<double>& grid2dD,
			math::Matrix<double>& grid2dR,
			math::Matrix<double>& grid2dV,
			math::Matrix<double>& grid2dRE
		)
		{
			int sizeL = gridL.size(), sizeXi = gridXi.size();
			math::Matrix<double> grid2dL(sizeL, sizeXi);
			for (int iL = 1; iL < sizeL - 1; iL++) {//omit bdy values
				double Jz, Jc = gridL[iL];
				//Find E, R, V of circular orbit with Jphi=Jc
					//double E = E_circ(pot, Jc, &Rc, &Vc);
				double E = gridE[iL];
				std::vector<double> L_vals(sizeXi);
				std::vector<double> Xi_vals(sizeXi);
				std::vector<double> D_vals(sizeXi);
				std::vector<double> V_vals(sizeXi);
				//				std::vector<double> ThetaT_vals(sizeXi);
				std::vector<double> R_vals(sizeXi);
				//Run over inclinations of orbits with this L
				for (int iXi = 1; iXi < sizeXi - 1; iXi++) {
					double Jphi = fabs(Jc * (1 - gridXi[iXi]));
					double Rsh, FD;
					std::vector<coord::PosVelCyl> shell;
					FD = estimateFocalDistanceShellOrbit(pot, E, Jphi, &Rsh, &Jz, &shell);
					double L = Jphi + Jz;
					L_vals[iXi] = L; Xi_vals[iXi] = Jz / L;
					D_vals[iXi] = FD; R_vals[iXi] = Rsh;
					V_vals[iXi] = shell[0].vz;
					if (iXi > 1 && Xi_vals[iXi] < Xi_vals[iXi - 1])
						printf("createGridFocalDistances: non-monotonic Xi:\n",
							"%d %f %f\n", iXi, Xi_vals[iXi - 1], Xi_vals[iXi]);
				}
				// bdy values
				L_vals[0] = gridL[iL]; L_vals[sizeXi - 1] = L_vals[sizeXi - 2];
				Xi_vals[0] = 0;       Xi_vals[sizeXi - 1] = 1;
				D_vals[0] = D_vals[1]; D_vals[sizeXi - 1] = D_vals[sizeXi - 2];
				R_vals[0] = R_vals[1]; R_vals[sizeXi - 1] = R_vals[sizeXi - 2];
				V_vals[0] = 0;         V_vals[sizeXi - 1] = V_vals[sizeXi - 2];
				// make interpolators L(L,xi), D(L,xi),
				// Rsh(L,xi), Vsh(L,Xi)
				math::LinearInterpolator interpL(Xi_vals, L_vals);
				math::LinearInterpolator interpD(Xi_vals, D_vals);
				math::LinearInterpolator interpR(Xi_vals, R_vals);
				math::LinearInterpolator interpV(Xi_vals, V_vals);
				for (int iXi = 0; iXi < sizeXi; iXi++) {//interpolate L, D onto regular grid in Xi
					interpL.evalDeriv(gridXi[iXi], &grid2dL(iL, iXi));
					interpD.evalDeriv(gridXi[iXi], &grid2dD(iL, iXi));
					interpR.evalDeriv(gridXi[iXi], &grid2dR(iL, iXi));
					interpV.evalDeriv(gridXi[iXi], &grid2dV(iL, iXi));
				}
			}
			//We now need to fill in rows iL=0, iL=sizeL-1
			for (int iXi = 0; iXi < sizeXi; iXi++) {
				grid2dL(0, iXi) = 0; grid2dL(sizeL - 1, iXi) = 1.01 * grid2dL(sizeL - 2, iXi);
				grid2dD(0, iXi) = 0; grid2dD(sizeL - 1, iXi) = grid2dD(sizeL - 2, iXi);
				grid2dR(0, iXi) = 0; grid2dR(sizeL - 1, iXi) = grid2dR(sizeL - 2, iXi);
				grid2dV(0, iXi) = 0; grid2dV(sizeL - 1, iXi) = grid2dV(sizeL - 2, iXi);
			}
			//Now grid2dD contains D on grids in E,Xi
			for (int iL = 0; iL < sizeL; iL++)
				for (int iXi = 0; iXi < sizeXi; iXi++)
					grid2dRE(iL, iXi) = grid2dR(iL, iXi);
			//grid2dD contains D on regular grid in Xi but irregular
			//values of L that are stored in grid2dL
			for (int iXi = 0; iXi < sizeXi; iXi++) {
				std::vector<double> L_vals(sizeL);
				std::vector<double> D_vals(sizeL);
				std::vector<double> R_vals(sizeL);
				std::vector<double> V_vals(sizeL);
				for (int iL = 0; iL < sizeL; iL++) {
					L_vals[iL] = grid2dL(iL, iXi);
					D_vals[iL] = grid2dD(iL, iXi);
					R_vals[iL] = grid2dR(iL, iXi);
					V_vals[iL] = grid2dV(iL, iXi);
					if (iL > 0 && L_vals[iL] <= L_vals[iL - 1])
						printf("createGridFocalDistance: non-monotonic L_vals %g %g\n",
							L_vals[iL - 1], L_vals[iL]);
				}
				math::LinearInterpolator DL(L_vals, D_vals);
				math::LinearInterpolator RL(L_vals, R_vals);
				math::LinearInterpolator VL(L_vals, V_vals);
				for (int iL = 0; iL < sizeL; iL++) {
					DL.evalDeriv(gridL[iL], &grid2dD(iL, iXi));
					RL.evalDeriv(gridL[iL], &grid2dR(iL, iXi));
					VL.evalDeriv(gridL[iL], &grid2dV(iL, iXi));
				}
			}
		}
		// H & its deriv of H wrt R,z,p_R,p_z,p_phi
		double H_dHdX(const potential::BasePotential& pot, const coord::PosMomCyl Rzphi,
			coord::PosMomCyl& dHdX) {
			double Phi; coord::GradCyl grad;
			pot.eval(Rzphi, &Phi, &grad);
			dHdX.R = grad.dR - pow_2(Rzphi.pphi / Rzphi.R) / Rzphi.R;
			dHdX.z = grad.dz; dHdX.phi = grad.dphi;
			dHdX.pR = Rzphi.pR; dHdX.pz = Rzphi.pz; dHdX.pphi = Rzphi.pphi / pow_2(Rzphi.R);
			return .5 * (pow_2(Rzphi.pR) + pow_2(Rzphi.pz) + pow_2(Rzphi.pphi / Rzphi.R)) + Phi;
		}
		//Class used in choice of isochrone
		class Iso {
		private:
			double L, Jr;
		public:
			Iso(double _L, double _Jr) : L(_L), Jr(_Jr) {}
			double b2cE(double Js) const {//general E
				return .5 * pow_2(Js * Js) / pow_2(Jr + .5 * (L + sqrt(L * L + 4 * Js * Js)));
			}
			double b2cEc(double Js) const {//circular E
				return .5 * pow_2(Js * Js) / pow_2(.5 * (L + sqrt(L * L + 4 * Js * Js)));
			}
			double g(double Js) const {//Rsh=b*g(Js)
				double g2 = (pow((L + sqrt(L * L + 4 * Js * Js)) / (2 * Js), 4) - 1);
				if (g2 < 0) printf("g2<0: %f\n", g2);
				return sqrt(g2);
			}
			double cob(double Js) const {//c/b
				return .5 * pow_2(Js) / b2cE(Js) - 1;;
			}
			double ecc(double Js) const {
				double boc = 1 / cob(Js);
				double e2 = 1 - pow_2(L / Js) * boc * (1 + boc);
				return e2 < 1 ? sqrt(e2) : 0;
			}
			double f(double b, double Js, double& e) const {//ratio of forces aopo/peri
				double cb = cob(Js), boc = 1 / cb;
				e = ecc(Js);
				double up = 1 + e, um = 1 - e;
				double rp = b * cb * sqrt(up * (up + 2 * boc)), ap = sqrt(b * b + rp * rp);
				double rm = b * cb * sqrt(um * (um + 2 * boc)), am = sqrt(b * b + rm * rm);
				return pow_2((b + am) / (b + ap)) * am / ap * rp / rm;
			}
		};
		/*
		 * Class to pick TM by minimising Sum (H-E)^2 around ToyMap torus
		*/
		class DeltaFinder : public math::IFunctionNdimDeriv {
		private:
			const potential::BasePotential& pot;
			const Actions J;
			const double E;
			const int nrmax, nzmax, npts;
		public:
			DeltaFinder(const potential::BasePotential& _pot,
				const Actions _J, const double _E) :
				pot(_pot), J(_J), E(_E),
				nrmax(4), nzmax(4), npts(nrmax* nzmax) {
			}
			virtual unsigned int numVar() const {
				return 3;
			}
			virtual unsigned int numValues() const {
				return npts;
			}
			virtual void evalDeriv(const double params[],
				double* dHvalues, double* dHdParams) const {
				ToyMap TM(Isochrone(params[0], params[1]),
					PointTrans(coord::UVSph(params[2])));
				std::vector<double> Hvalues(npts);
				for (int i = 0; i < nrmax; i++) {
					for (int j = 0; j < nzmax; j++) {
						int ind = i * nzmax + j;
						Angles aT(i * M_PI / (double)nrmax, j * .5 * M_PI / (double)nzmax, 0);
						ActionAngles aa(J, aT);
						coord::PosMomCyl dHdR, dRzdPs[3], dHdPs[3];
						coord::PosMomCyl Rz(TM.from_aaT(aa, dRzdPs));
						Hvalues[ind] = H_dHdX(pot, Rz, dHdR) - E;
						if (dHdParams)
							for (int k = 0; k < 3; k++)
								dHdParams[3 * ind + k] = dHdR.R * dRzdPs[k].R + dHdR.z * dRzdPs[k].z
								+ dHdR.pR * dRzdPs[k].pR + dHdR.pz * dRzdPs[k].pz;
					}
				}
			}
			virtual unsigned int numDerivs(void) const {
				return 0;
			}
		};
		/*
		 Class to find Js (and thus b) of isochrone for which the isochrone's
		 force ratio matches that in real Phi
		*/
		class JsFinder : public math::IFunction {
			const potential::BasePotential& pot;
			const Iso Is;
			const double Rsh;
		public:
			JsFinder(const potential::BasePotential& _pot, const Iso& _Is, const double _Rsh) :
				pot(_pot), Is(_Is), Rsh(_Rsh) {
			}
			virtual void evalDeriv(double Js, double* val, double* deriv = 0, double* deriv2 = 0) const {
				double b = Rsh / Is.g(Js);
				double e, f_apo_peri = Is.f(b, Js, e);
				double c = Is.cob(Js) * b;
				double F[2];
				for (int k = -1; k < 2; k += 2) {
					double u = 1 + k * e;
					double Phi, r = c * sqrt(u * (u + 2 * b / c));
					coord::PosCyl Rz(r, 0, 0); coord::GradCyl grad;
					pot.eval(Rz, &Phi, &grad);
					F[(k + 1) / 2] = grad.dR;
				}
				if (val) *val = f_apo_peri - F[1] / F[0];
			}
			virtual unsigned int numDerivs(void) const {
				return 0;
			}

		};

		/// create the array of indices of the generating function with all terms up to the given maximum order
		static GenFncIndices makeGridIndices(int irmax, int izmax)
		{   /// NOTE: here we specialize for the case of axisymmetric systems!
			GenFncIndices indices;
			for (int ir = 0; ir <= irmax; ir++)
				for (int iz = -izmax; iz <= (ir == 0 ? -2 : izmax); iz += 2)
					indices.push_back(GenFncIndex(ir, iz, 0));
			return indices;
		}

		/// return the absolute value of an element in a map, or zero if it doesn't exist
		static inline double absvalue(const std::map< std::pair<int, int>, double >& indPairs, int ir, int iz)
		{
			if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end())
				return fabs(indPairs.find(std::make_pair(ir, iz))->second);
			else
				return 0;
		}
		/** Compute the derivative of Hamiltonian by toy actions:
				dH/dJ = dH/d{x,v} d{x,v}/dJ, where the lhs is a covector of length 3,
				the first term on rhs is a covector of length 6 (the gradient dPhi/dx and the velocity),
				and the second term is a 6x3 matrix of partial derivs provided by the toy map.
		*/
		static inline Actions dHbydJ(const coord::PosMomCyl& dHby,
			const DerivAct<coord::Cyl>& dXdJ)
		{
			return Actions(
				dXdJ.dbyJr.R * dHby.R + dXdJ.dbyJr.z * dHby.z + dXdJ.dbyJr.phi * dHby.phi +
				dXdJ.dbyJr.pR * dHby.pR + dXdJ.dbyJr.pz * dHby.pz + dXdJ.dbyJr.pphi * dHby.pphi,
				dXdJ.dbyJz.R * dHby.R + dXdJ.dbyJz.z * dHby.z + dXdJ.dbyJz.phi * dHby.phi +
				dXdJ.dbyJz.pR * dHby.pR + dXdJ.dbyJz.pz * dHby.pz + dXdJ.dbyJz.pphi * dHby.pphi,
				dXdJ.dbyJphi.R * dHby.R + dXdJ.dbyJphi.z * dHby.z + dXdJ.dbyJphi.phi * dHby.phi +
				dXdJ.dbyJphi.pR * dHby.pR + dXdJ.dbyJphi.pz * dHby.pz + dXdJ.dbyJphi.pphi * dHby.pphi);
		}
		//class to help with zSoS computation
		class CrossingFinder : public math::IFunction {
		private:
			const Torus* T;
			const double thetaT_r, z0;
		public:
			CrossingFinder(const Torus* _T, const double& _thetaT, const double _z0 = 0) :
				T(_T), thetaT_r(_thetaT), z0(_z0) {
			}
			virtual void evalDeriv(const double thetaT_z, double* value, double* deriv = NULL, double* deriv2 = NULL) const {
				Angles thetaT(thetaT_r, thetaT_z, 0);
				coord::PosMomCyl Rz(T->from_toy(thetaT));
				*value = Rz.z - z0;
			}
			virtual unsigned int numDerivs(void) const {
				return 0;
			}
		};
		//class to help with RSoS computation
		class r_crossingFinder : public math::IFunction {
		private:
			const Torus* T;
			const double rbar;
			const double thetaT_z;
		public:
			r_crossingFinder(const Torus* _T, const double _rbar, const double& _thetaT) :
				T(_T), rbar(_rbar), thetaT_z(_thetaT) {
			}
			virtual void evalDeriv(const double thetaT_r, double* value, double* deriv = NULL, double* deriv2 = NULL) const {
				Angles thetaT(thetaT_r, thetaT_z, 0);
				coord::PosMomSph rtheta(coord::toPosMomSph(T->from_toy(thetaT)));
				*value = rtheta.r - rbar;
			}
			virtual unsigned int numDerivs(void) const {
				return 0;
			}
		};
		/*
		 * Function to comput dx,dp/d true angles. It
		 * takes in (dx,dp)/dthetaT and dtheta/dthetaT, inverts the
		 * latter and multiplies on the former
		 */
		void new_assemble(DerivAng<coord::Cyl>& dA, math::Matrix<double>& M) {
			math::Matrix<double> dxdt(3, 3);
			dxdt(0, 0) = dA.dbythetar.R; dxdt(0, 1) = dA.dbythetaz.R;
			dxdt(0, 2) = dA.dbythetaphi.R;
			dxdt(1, 0) = dA.dbythetar.z; dxdt(1, 1) = dA.dbythetaz.z;
			dxdt(1, 2) = dA.dbythetaphi.z;
			dxdt(2, 0) = dA.dbythetar.phi; dxdt(2, 1) = dA.dbythetaz.phi;
			dxdt(2, 2) = dA.dbythetaphi.phi;
			math::Matrix<double> dpdt(3, 3);
			dpdt(0, 0) = dA.dbythetar.pR; dpdt(0, 1) = dA.dbythetaz.pR;
			dpdt(0, 2) = dA.dbythetaphi.pR;
			dpdt(1, 0) = dA.dbythetar.pz; dpdt(1, 1) = dA.dbythetaz.pz;
			dpdt(1, 2) = dA.dbythetaphi.pz;
			dpdt(2, 0) = dA.dbythetar.pphi; dpdt(2, 1) = dA.dbythetaz.pphi;
			dpdt(2, 2) = dA.dbythetaphi.pphi;
			math::LUDecomp L2(M);
			math::Matrix<double> M3 = (L2.inverse(3));
			math::Matrix<double> M2(3, 3);
			math::blas_dgemm(math::CblasNoTrans, math::CblasNoTrans, 1.0,
				dxdt, M3, 0.0, M2);
			dA.dbythetar.R = M2(0, 0);
			dA.dbythetaz.R = M2(0, 1);
			dA.dbythetaphi.R = M2(0, 2);
			dA.dbythetar.z = M2(1, 0);
			dA.dbythetaz.z = M2(1, 1);
			dA.dbythetaphi.z = M2(1, 2);
			dA.dbythetar.phi = M2(2, 0);
			dA.dbythetaz.phi = M2(2, 1);
			dA.dbythetaphi.phi = M2(2, 2);
			math::blas_dgemm(math::CblasNoTrans, math::CblasNoTrans, 1.0,
				dpdt, M3, 0.0, M2);
			dA.dbythetar.pR = M2(0, 0);
			dA.dbythetaz.pR = M2(0, 1);
			dA.dbythetaphi.pR = M2(0, 2);
			dA.dbythetar.pz = M2(1, 0);
			dA.dbythetaz.pz = M2(1, 1);
			dA.dbythetaphi.pz = M2(1, 2);
			dA.dbythetar.pphi = M2(2, 0);
			dA.dbythetaz.pphi = M2(2, 1);
			dA.dbythetaphi.pphi = M2(2, 2);
		}

		double avgH;
		void setH(double Hbar) { avgH = Hbar; }

		/* A helper class used by TorusGenerator. It evaluates H and its
		 * derivatives on a toy-angle grid over a proposed torus using a
		 * proposed GF. Once the GF has been
		 * optimised, it runs over the grid a final time to determine d S_k/d J
		 * so we can recover the true angles
		 */
		class torusFitter : public math::IFunctionNdimDeriv {
		private:
			const Actions& J;
			const potential::BasePotential& pot;
			GenFncFitSeries& GFFS;
			const ToyMap TM;
		public:
			torusFitter(const Actions& _J,
				const potential::BasePotential& _pot,
				const ToyMap _TM, GenFncFitSeries& _GFFS) :
				J(_J), pot(_pot), TM(_TM), GFFS(_GFFS) {
			}
			virtual unsigned int numVars() const {
				return GFFS.numParams();
			}
			virtual unsigned int numValues() const {
				return GFFS.numPoints();
			}
			double new_computeHamiltonianAtPoint(const double params[],
				const unsigned int indPoint, Actions& JT, Actions* dHdJ = NULL,
				double* dHdParams = NULL) const
			{
				// Generating function computes the toy actions from the real actions
				// at the given point in the grid of toy angles grid
				ActionAngles toyAA = GFFS.toyActionAngles(indPoint, params);

				// do not allow to stray into forbidden region of negative actions
				if (toyAA.Jr < 0 || toyAA.Jz < 0) {
					JT.Jr = toyAA.Jr; JT.Jz = toyAA.Jz;
					return 0;
				}
				DerivAct<coord::Cyl> dXdJ;
				coord::PosMomCyl Rzphi(TM.from_aaT(toyAA, dXdJ));
				//printf("R,z.. %f %f %f %f %f %f\n",Rzphi.R,Rzphi.z,Rzphi.pR,Rzphi.pz,Rzphi.pphi);
		// obtain the value of the real Hamiltonian at the given point and its
		// derivatives w.r.t. coordinates/momenta
				coord::PosMomCyl dHdX;
				double H = H_dHdX(pot, Rzphi, dHdX);
				if (dHdParams) {
					// derivative of Hamiltonian by toy actions
					Actions dHby = dHbydJ(dHdX, dXdJ);
					if (dHdJ) *dHdJ = dHby;
					for (unsigned int p = 0; p < GFFS.numParams(); p++) {
						// derivs of toy actions by gen.fnc.params
						Actions dbyP = GFFS.deriv(indPoint, p, &params[0]);
						dHdParams[p] = dHby.Jr * dbyP.Jr + dHby.Jz * dbyP.Jz
							+ dHby.Jphi * dbyP.Jphi;
					}
				}
				return H;
			}
			double computeHamiltonianDisp(const std::vector<double>& params, double& Hbar)
			{
				const unsigned int numParams = GFFS.numParams(), max_threads = 8;
				double Hsm[max_threads] = { 0 }, Hsq[max_threads] = { 0 };
				int N[max_threads] = { 0 };
				int Nnan[max_threads] = { 0 };
#pragma omp parallel for schedule(dynamic)
				for (int indPoint = 0; indPoint < GFFS.numPoints(); indPoint++) {
					int nth = omp_get_thread_num();
					Actions JT(1, 1, 1);
					double H = new_computeHamiltonianAtPoint(&params[0], indPoint, JT);
					if (JT.Jr < 0 || JT.Jz < 0) Nnan[nth]++;
					else {
						Hsm[nth] += H; Hsq[nth] += H * H; N[nth]++;
					}
				}
				for (int i = 1; i < max_threads; i++) {//concatenate sums
					Nnan[0] += Nnan[i]; Hsm[0] += Hsm[i]; Hsq[0] += Hsq[i]; N[0] += N[i];
				}
				Hbar = Hsm[0] / N[0];
				Hsq[0] = Hsq[0] / N[0] - Hbar * Hbar;
				NANfrac = (double)Nnan[0] / (double)GFFS.numPoints();
				NANbar += NANfrac;
				return sqrt(Hsq[0]);
			}
			void evalDeriv(const double params[],
				double* deltaHvalues, double* dHdParams) const
			{
				const unsigned int numPoints = GFFS.numPoints();
				const unsigned int numParams = GFFS.numParams();
				const unsigned int max_threads = 16;

				// we need to store the values of Hamiltonian at grid points even if this is not requested,
				// because they are used to correct the entries of the Jacobian matrix
				// to account for the fact that the mean <H> also depends on the parameters
				std::vector<double> Hvalues(numPoints);
				double Havg[max_threads] = { 0 };  // accumulator for the average Hamiltonian
				int nNAN[max_threads] = { 0 };
				// loop over grid of toy angles
#pragma omp parallel for schedule(dynamic)
				for (int indPoint = 0; indPoint < numPoints; indPoint++) {
					int nth = omp_get_thread_num();
					Actions JT(1, 1, 1);
					double H = new_computeHamiltonianAtPoint(params, indPoint, JT, NULL,
						dHdParams ? dHdParams + indPoint * numParams : NULL);
					if (JT.Jr < 0 || JT.Jz < 0) {
						nNAN[nth]++;
						Havg[nth] += 1e6;
						Hvalues[indPoint] = 1e6; break;
					}
					else {
						// accumulate the average value and store the output
						Havg[nth] += H;
						Hvalues[indPoint] = H;
					}
				}
				for (int i = 1; i < max_threads; i++) {
					nNAN[0] += nNAN[i]; Havg[0] += Havg[i];
				}

				// convert from  H_k  to  deltaH_k = H_k - <H>
				NANfrac = (double)nNAN[0] / numPoints;//we can't put nNAN in NANfrac!
				NANbar += NANfrac;
				Havg[0] /= numPoints; setH(Havg[0]);
				if (deltaHvalues) {
					double disp = 0;
					for (unsigned int indPoint = 0; indPoint < numPoints; indPoint++) {
						deltaHvalues[indPoint] = Hvalues[indPoint] - Havg[0];
						disp += pow_2(deltaHvalues[indPoint]);
					}
				}
				// convert derivatives:  d(deltaH_k) / dP_p = dH_k / dP_p - d<H> / dP_p
				if (dHdParams) {
					std::vector<double> dHavgdP(numPoints);
					for (unsigned int p = 0; p < numParams; p++) dHavgdP[p] = 0;
					for (unsigned int pp = 0; pp < numPoints * numParams; pp++) {
						dHavgdP[pp % numParams] += dHdParams[pp] / numPoints;
						if (std::isnan(dHavgdP[pp % numParams])) {
							printf("nan@ %d %d %g\n", pp % numParams, pp / numParams, dHdParams[pp]);
							exit(0);
						}
					}
					for (unsigned int pp = 0; pp < numPoints * numParams; pp++) {
						unsigned int indPoint = pp / numParams;
						unsigned int indParam = pp % numParams;
						dHdParams[pp] = dHdParams[pp] - dHavgdP[indParam];
					}
				}
			}
			void testit(std::vector<double> params, int N) {
				int nPts = GFFS.numPoints(), nPars = GFFS.numParams();
				double varH;
				double* HmHbar = new double[nPts];
				double* HmX = new double[nPts];
				double* dHdP = new double[nPts * nPars];
				double dA = 1e-5;
				//explore changing Nth parameter 
				params[N] += .5 * dA;
				evalDeriv(&params[0], HmHbar, NULL);
				params[N] -= dA;
				evalDeriv(&params[0], HmX, NULL);
				for (int i = 0;i < nPts;i++)
					HmX[i] = (HmHbar[i] - HmX[i]) / dA;
				params[N] += .5 * dA;
				evalDeriv(&params[0], HmHbar, dHdP);
				for (int i = 0; i < nPts; i += 10)
					printf("(%g %g) ", HmX[i], dHdP[i * nPars + N]);
				printf("\n");
				delete[] HmHbar; delete[] HmX; delete[] dHdP;
			}
			/** Compute the frequencies and the derivatives of generating function by real actions,
				used in angle mapping.
				The three arrays of derivatives dS_i/dJ_{r,z,phi}, i=0..numParamsGenFnc-1,
				together with three frequencies Omega_{r,z,phi}, are the solutions of
				an overdetermined system of linear equations:
				\f$  M_{k,i} X_{i} = RHS_{k}, k=0..numPoints-1  \f$,
				where numPoints is the number of individual triplets of toy angles,
				\f$  X_i  \f$ is the solution vector {Omega, dS_i/dJ} for each direction (r,z,theta),
				\f$  RHS_k = dH/dJ(\theta_k)  \f$, again for three directions independently, and
				\f$  M_{k,i}  \f$ is the matrix of coefficients shared between all three equation systems:
				\f$  M_{k,0} = 1, M_{k,i+1} = -dH/dS_i(\theta_k)  \f$.
				The matrix M and three RHS vectors are filled using the same approach as during
				the Levenberg-Marquardt minimization, from the provided parameters of toy map and
				generating function; then the three linear systems are solved using
				the singular-value decomposition of the shared coefficient matrix,
				and the output frequencies and gen.fnc.derivatives are returned in corresponding arguments.
				The return value of this function is the same as `computeHamiltonianDisp()`.
			*/
			double fitAngleMap(const double params[],
				double& Hbar, Frequencies& freqs, GenFncDerivs& dPdJ) const {
				unsigned int numPoints = GFFS.numPoints();
				unsigned int numParams = GFFS.numParams();
				// the matrix of coefficients shared between three linear systems
				math::Matrix<double> coefsdHdS(numPoints, numParams + 1);
				// tmp storage for dH/dS
				std::vector<double> dHdParams(numParams + 1);
				// derivs of Hamiltonian by toy actions (RHS vectors)
				std::vector<double> dHdJr(numPoints), dHdJz(numPoints), dHdJphi(numPoints);
				// accumulator for computing dispersion in H
				math::Averager Havg;

				// loop over grid of toy angles
				for (unsigned int indPoint = 0; indPoint < numPoints; indPoint++) {
					Actions JT(1, 1, 1), dHby;  // derivative of Hamiltonian by toy actions
					double H = new_computeHamiltonianAtPoint(&params[0], indPoint, JT,
						&dHby, &dHdParams.front());
					if (!(JT.Jr < 0 || JT.Jz < 0)) Havg.add(H);
					else printf("J<0 @ %d ", indPoint);
					// fill the elements of each of three rhs vectors
					dHdJr[indPoint] = dHby.Jr;
					dHdJz[indPoint] = dHby.Jz;
					dHdJphi[indPoint] = dHby.Jphi;
					// fill the matrix row
					coefsdHdS(indPoint, 0) = 1;  // matrix coef for omega
					for (unsigned int p = 0; p < numParams; p++) // matrix coefs for dSk/dJ
						coefsdHdS(indPoint, p + 1) = -dHdParams[p];
				}
				Hbar = Havg.mean();
				// solve the overdetermined linear system in the least-square sense:
				// step 1: prepare the SVD of coefs matrix
				math::SVDecomp SVD(coefsdHdS);

				// step 2: solve three linear systems with the same matrix but different rhs
				std::vector<double> dSdJr(SVD.solve(dHdJr)), dSdJz(SVD.solve(dHdJz)),
					dSdJphi(SVD.solve(dHdJphi));

				// store output
				freqs.Omegar = dSdJr[0];
				freqs.Omegaz = dSdJz[0];
				freqs.Omegaphi = J.Jphi>0? dSdJphi[0] : -dSdJphi[0];
				dPdJ.resize(numParams);
				for (unsigned int p = 0; p < numParams; p++) {
					dPdJ[p].Jr = dSdJr[p + 1];
					dPdJ[p].Jz = dSdJz[p + 1];
					dPdJ[p].Jphi = dSdJphi[p + 1];
				}
				return sqrt(Havg.disp());
			}

		};
		
/* Helper class to fit v(theta) for the PointTrans by minimising dH */
		class fitvtmap : public math::IFunctionNdimDeriv {
			const int N;
			GenFncFitSeries GFFS;
			const int np;
			const Actions J;
			const Isochrone iso;
			double delta;
			const potential::BasePotential& pot;
			std::vector<double> paramsGF;
		public:
			fitvtmap(const int _N, GenFncFitSeries _GFFS, Actions _J,
				 Isochrone _iso, double _delta,
				 const potential::BasePotential& _pot,
				 std::vector<double> _params)
					: N(_N), GFFS(_GFFS), np(GFFS.numPoints()),
			J(_J), iso(_iso), delta(_delta), pot(_pot), paramsGF(_params) {
			}
			virtual unsigned int numVars() const {
				return N;
			}
			virtual unsigned int numValues() const {
				return np;
			}
			double t2v(const double params1[], double X, double* dfdParams = NULL, double* deriv = NULL, double* ddfdparams = NULL, double* d2fdx2 = NULL) const {
				double f1 = X;
				double f2 = 1.0, f3 = 0.0;
				for (int i = 0;i < N;i++) {
					f1 += params1[i] * sin(2 * (i + 1) * X);
					if (dfdParams) {
						dfdParams[i] = sin(2 * (i + 1) * X);
					}
					if (deriv)f2 += 2 * (i + 1) * params1[i] * cos(2 * (i + 1) * X);
					if (ddfdparams)ddfdparams[i] = 2 * (i + 1) * cos(2 * (i + 1) * X);
					if (d2fdx2)f3 += -pow_2(2 * (i + 1)) * params1[i] * sin(2 * (i + 1) * X);
				}
				if (deriv) *deriv = f2;
				if (d2fdx2)*d2fdx2 = f3;
				return f1;
			}
			double H(const double params[], int i, double* dHdParams = NULL, double* dHdtheta0 = NULL, double* dHdptheta0 = NULL) const {
				double dvdt;std::vector<double> dthetadparams(N);std::vector<double> ddthetadparams(N);
				ActionAngles aaT=GFFS.toyActionAngles(i,&paramsGF[0]);
				coord::PosMomSph xv = iso.aa2pq(aaT);
				PointTrans PT(delta);
				double thetan = t2v(&params[0], xv.theta, &dthetadparams[0], &dvdt, &ddthetadparams[0]);
				double theta0 = xv.theta;
				xv.theta = thetan;
				xv.ptheta /= dvdt;
				coord::PosMomCyl Rp = PT.Sph2Cyl(xv);
				coord::PosMomCyl dHdx;

				double E = H_dHdX(pot, Rp, dHdx);
				if (dHdParams) {
					double sq = pow_2(xv.r) + pow_2(delta), rt = sqrt(sq);
					double snt, cst; math::sincos(thetan, snt, cst);
					double R = xv.r * snt, z = rt * cst;
					double bot1 = pow_2(xv.r) + pow_2(delta) * pow_2(snt); //(sq * pow_2(snt) + pow_2(rp.r * cst));
					double dbot1dr = 2 * xv.r;
					double dbot1dtheta = 2 * pow_2(delta) * cst * snt;
					double bot2 = (pow_2(xv.r * cst) / rt + rt * pow_2(snt));
					double dbot2dr = -xv.r / pow_3(rt) * bot1 + 2 * xv.r / rt;
					double dbot2dtheta = dbot1dtheta / rt;
					double pR = (sq * snt * xv.pr + xv.r * cst * xv.ptheta) / bot1;
					double pz = (xv.r * cst * xv.pr - snt * xv.ptheta) / bot2;
					double dRdr = snt; double dRdt = xv.r * cst; 
					double dzdr= xv.r / rt * cst;   double dzdt= -rt * snt;
					double dprdr = (2 * xv.r * snt * xv.pr + cst * xv.ptheta) / bot1
						- pR / bot1 * dbot1dr; //dpR/dr 
					double dpRdt = ((sq * cst * xv.pr - xv.r * snt * xv.ptheta) / bot1
						- pR / bot1 * dbot1dtheta) ;//dpR/dpsi
					double dpzdr = cst * xv.pr / bot2 - pz / bot2 * dbot2dr;// dpz/dr
					double dpzdt = ((-xv.r * snt * xv.pr - cst * xv.ptheta) / bot2
						- pz / bot2 * dbot2dtheta);
					double dpRdpr = sq * snt / bot1;double dpRdpt = xv.r * cst / (bot1);//dpR/dpr ppsi
					double dpzdpr = xv.r * cst / bot2; double dpzdpt = -snt / (bot2);//dpz/dpr ppsi
					double dHdtheta = dpRdt * dHdx.pR + dpzdt * dHdx.pz + dHdx.R * dRdt + dHdx.z * dzdt;
					double dHdptheta = dpRdpt * dHdx.pR + dHdx.pz * dpzdpt;
					if (dHdtheta0)*dHdtheta0 = dHdtheta;
					if (dHdptheta0)*dHdptheta0 = dHdptheta;
					for (int i = 0;i < N;i++) {
						dHdParams[i] = dHdtheta * dthetadparams[i] + dHdptheta / dvdt * xv.ptheta * (-ddthetadparams[i]);
					}
				}
				return E;

			}
			virtual void evalDeriv(const double params[],
				double* dH, double* dfdParams)const
			{
				std::vector<double>Hval(np);
				std::vector<double> dfdparams0(N * np);

				double Hav = 0.0;
				std::vector<double> dfdav(N, 0.0);
				std::vector<double> dfdp(N);
				for (int i = 0;i < np;i++) {
					double E = H(&params[0], i, &dfdp[0]);
					Hav += E;
					Hval[i] = E;
					if (dfdParams) {
						for (int j = 0;j < N;j++) {
							dfdparams0[i * N + j] = dfdp[j];
							dfdav[j] += dfdp[j];
						}
					}
					//}
				}
				Hav /= np;
				for (int i = 0;i < N;i++) {
					dfdav[i] /= np;
				}
				//double disp = 0.0;
				if (dH) {
					for (int i = 0;i < np;i++) {
						dH[i] = Hval[i] - Hav;
						//disp += pow_2(dH[i]) / np;
					}
				}
				if (dfdParams) {
					for (int i = 0;i < np;i++) {
						for (int j = 0;j < N;j++) {
							dfdParams[i * N + j] = (dfdparams0[i * N + j] - dfdav[j]);
						}
					}
				}
			}
		};

		//Finds where x lies in xs[]. returns fractional distance from xs[top]
		double bot_top(const double x, const std::vector<double>& xs,
			       int& top, int& bot){
			if (xs[0] > xs[xs.size() - 1]) {
				top = 0; bot = xs.size() - 1;
			}
			else {
				top = xs.size() - 1; bot = 0;
			}//now top should point to largest x
			double f;
			if ((xs[top] - x) * (x - xs[bot]) < 0.) {//x lies out of grid
				if (x < xs[bot]) {
					top = bot + 1; f = 1;//T = Tgrid[bot];
				}
				else {
					bot = top - 1; f = 0;// T = Tgrid[top];
				}
			}
			else {
				while (abs(top - bot) > 1) {
					int n = (top + bot) / 2;
					if ((xs[top] - x) * (x - xs[n]) >= 0) bot = n;
					else top = n;
				}
				f = (xs[top] - x) / (xs[top] - xs[bot]);//distance from top
			}
			return f;
		}


	}//internal

	//function used to invert v(theta) to get theta(v)
	class vthet : public math::IFunction {
	private:
		std::vector<double> p;
		double v;
	public:
		vthet(std::vector<double> _p, double _v) :v(_v), p(_p) {};
		virtual unsigned int numDerivs()const { return 1; }
		virtual void evalDeriv(const double x,
			double* value, double* deriv, double* deriv2) const {
			double f1 = x - v;
			double f2 = 1.;
			for (int i = 0;i < p.size();i++) {
				if(value)f1 += sin(2 * (i + 1) * x) * p[i];
				if(deriv)f2 += 2 * (i + 1) * cos(2 * (i + 1) * x) * p[i];
			}
			if (value)*value = f1;
			if (deriv)*deriv = f2;

		}
	};
	EXP PointTrans interpPointTrans(double x, const PointTrans& PT0, const PointTrans& PT1) {
		double xp = 1-x;
		coord::UVSph cs1(x * PT0.cs.Delta + xp * PT1.cs.Delta);
		if(!(PT0.map || PT1.map))//neither PT has paramsF
			return PointTrans(cs1);
		std::vector<double> p;
		if(PT0.N == PT1.N)
			for(int i=0; i<PT0.N; i++)
				p.push_back(x * PT0.paramsF[i] + xp * PT1.paramsF[i]);
		else if(PT0.N>PT1.N){
			for(int i=0; i<PT1.N; i++)
				p.push_back(x * PT0.paramsF[i] + xp * PT1.paramsF[i]);
			for(int i=PT1.N; i<PT0.N; i++)
				p.push_back(x * PT0.paramsF[i]);
		} else {
			for(int i=0; i<PT0.N; i++)
				p.push_back(x * PT0.paramsF[i] + xp * PT1.paramsF[i]);
			for(int i=PT0.N; i<PT1.N; i++)
				p.push_back(xp * PT1.paramsF[i]);
		}
		return PointTrans(cs1,p);
	}
	std::vector<double> trivXs = { 0,.5 * M_PI };
	PointTrans::PointTrans() : map(false), paramsF({}) {
		N = 0;
	}
	PointTrans::PointTrans(double _D) :
		map(false), cs(_D), paramsF({}) {
		N = 0;
	}
	PointTrans::PointTrans(coord::UVSph _cs) :
		map(false), cs(_cs),paramsF({}) {
		N = 0;
	}

	double PointTrans::t2v(const double theta, double* dvdt,double* d2tdv2) const {
		double v = theta;
		double deriv = 1.;double deriv2 = 0.0;
		if (map) {
			for (int i = 0; i < N; i++) {
				v += paramsF[i] * sin(2 * (i + 1) * theta);
				deriv += 2*(i+1)*paramsF[i] * cos(2 * (i + 1) * theta);
				deriv2 += -pow_2(2 * (i + 1)) * paramsF[i] * sin(2 * (i + 1) * theta);
			}
		}
		if (dvdt)*dvdt = deriv;
		if (d2tdv2)*d2tdv2 = deriv2;
		return v;
	}

	double PointTrans::v2t(const double v, double* dtdv) const {
		if (!map) {
			if (dtdv) *dtdv = 1; return v;
		}
		vthet V(paramsF, v);
		double t = math::findRoot(V, 0., M_PI, 1e-9);
		if (dtdv) {
			double dvdt;
			t2v(t, &dvdt);
			*dtdv = 1. / dvdt;
		}
		return t;
	}
	coord::PosMomSph PointTrans::Cyl2Sph(const coord::PosMomCyl Rz) const {
		double R2 = pow_2(Rz.R), z2 = pow_2(Rz.z);
		double B = R2 + z2 - cs.Delta2;
		double r2 = .5 * (B + sqrt(B * B + 4 * R2 * cs.Delta2)), r = sqrt(r2);
		double rt = sqrt(r2 + cs.Delta2);
		double v = acos(Rz.z / rt), dtdv;
		double snt, cst; math::sincos(v, snt, cst);
		double pr = snt * Rz.pR + r / rt * cst * Rz.pz;
		double ptheta = r * cst * Rz.pR - rt * snt * Rz.pz;
		double psi = v2t(v, &dtdv);
		return coord::PosMomSph(r, psi, Rz.phi, pr, ptheta/dtdv, Rz.pphi);
	}

	coord::PosMomCyl PointTrans::Sph2Cyl(const coord::PosMomSph rp, coord::PosMomCyl* dRzdDelta) const {
		double sq = pow_2(rp.r) + cs.Delta2, rt = sqrt(sq);
		double dvdt;
		double v = t2v(rp.theta, &dvdt);
		double snt, cst; math::sincos(v, snt, cst);
		double R = rp.r * snt, z = rt * cst;
		double bot1 = (sq * pow_2(snt) + pow_2(rp.r * cst));
		double bot2 = (pow_2(rp.r * cst) / rt + rt * pow_2(snt));
		double pv = rp.ptheta / dvdt;
		double pR = (sq * snt * rp.pr + rp.r * cst * pv) / bot1;
		double pz = (rp.r * cst * rp.pr - snt * pv) / bot2;
		if (dRzdDelta) {
			dRzdDelta->R = 0; dRzdDelta->z = cs.Delta / rt * cst;
			dRzdDelta->pR = 2 * cs.Delta * snt * rp.pr / bot1
				- pR / bot1 * 2 * cs.Delta * pow_2(snt);
			dRzdDelta->pz = -pz / bot2 * cs.Delta / rt * (-pow_2(rp.r * cst) / sq
				+ pow_2(snt));
			dRzdDelta->phi = 0; dRzdDelta->pphi = 0;
		}
		return coord::PosMomCyl(R, z, rp.phi, pR, pz, rp.pphi);
	}

	coord::PosMomCyl PointTrans::Sph2Cyl(const coord::PosMomSph rp, math::Matrix<double>& dRzdrt,
		coord::PosMomCyl* dRzdDelta) const {
		double sq = pow_2(rp.r) + cs.Delta2, rt = sqrt(sq);
		double dvdt,d2vdt2;
		double v = t2v(rp.theta, &dvdt,&d2vdt2);
		double snt, cst; math::sincos(v, snt, cst);
		double R = rp.r * snt, z = rt * cst;
		double bot1 = pow_2(rp.r) + cs.Delta2 * pow_2(snt); //(sq * pow_2(snt) + pow_2(rp.r * cst));
		double dbot1dr = 2 * rp.r;
		double dbot1dtheta = 2 * cs.Delta2 * cst * snt;
		double bot2 = (pow_2(rp.r * cst) / rt + rt * pow_2(snt));
		double dbot2dr = -rp.r / pow_3(rt) * bot1 + 2 * rp.r / rt;
		double dbot2dtheta = dbot1dtheta / rt;
		double pR = (sq * snt * rp.pr + rp.r * cst * rp.ptheta/dvdt) / bot1;
		double pz = (rp.r * cst * rp.pr - snt * rp.ptheta/dvdt) / bot2;
		dRzdrt(0, 0) = snt; dRzdrt(0, 1) = rp.r * cst * dvdt; dRzdrt(0, 2) = 0;//dR/dr dR/dpsi
		dRzdrt(1, 0) = rp.r / rt * cst;   dRzdrt(1, 1) = -rt * snt * dvdt; dRzdrt(1, 2) = 0;//dz/dr dz/dpsi
		dRzdrt(2, 0) = (2 * rp.r * snt * rp.pr + cst * rp.ptheta/dvdt) / bot1
			- pR / bot1 * dbot1dr; //dpR/dr 
		dRzdrt(2, 1) = ((sq * cst * rp.pr - rp.r * snt * rp.ptheta/dvdt) / bot1
			- pR / bot1 * dbot1dtheta) * dvdt-rp.r*cst/bot1*rp.ptheta/pow_2(dvdt)*d2vdt2;//dpR/dpsi
		dRzdrt(3, 0) = cst * rp.pr / bot2 - pz / bot2 * dbot2dr;// dpz/dr
		dRzdrt(3, 1) = ((-rp.r * snt * rp.pr - cst * rp.ptheta/dvdt) / bot2
			- pz / bot2 * dbot2dtheta) * dvdt+snt/bot2*rp.ptheta/pow_2(dvdt)*d2vdt2;// dpz/dpsi
		dRzdrt(0, 2) = 0; dRzdrt(0, 3) = 0;//dR/dpr ppsi
		dRzdrt(1, 2) = 0; dRzdrt(1, 3) = 0;//dz/dpr ppsi
		dRzdrt(2, 2) = sq * snt / bot1; dRzdrt(2, 3) = rp.r * cst / (bot1*dvdt);//dpR/dpr ppsi
		dRzdrt(3, 2) = rp.r * cst / bot2; dRzdrt(3, 3) = -snt / (bot2 * dvdt);//dpz/dpr ppsi
		if (dRzdDelta) {
			dRzdDelta->R = 0; dRzdDelta->z = cs.Delta / rt * cst;
			dRzdDelta->pR = 2 * cs.Delta * snt * rp.pr / bot1
				- pR / bot1 * 2 * cs.Delta * pow_2(snt);
			dRzdDelta->pz = -pz / bot2 * cs.Delta / rt * (-pow_2(rp.r * cst) / sq
				+ pow_2(snt));
			dRzdDelta->phi = 0; dRzdDelta->pphi = 0;
		}
		return coord::PosMomCyl(R, z, rp.phi, pR, pz, rp.pphi);
	}

	coord::PosMomCar xyPointTrans::rp2xp(const coord::PosMomSph rp) const {
		double spsi, cpsi, dpsidphi; math::sincos(PC.Psi(rp.phi, dpsidphi), spsi, cpsi);
		double x = rp.r * cpsi, y = sqrt(pow_2(rp.r) + Delta2) * spsi;
		double r2 = pow_2(rp.r);
		double A = (r2 + Delta2 * pow_2(cpsi)) * dpsidphi, rt = sqrt(r2 + Delta2);
		double px = ((r2 + Delta2) * cpsi * dpsidphi * rp.pr - rp.r * spsi * rp.pphi) / A;
		double py = rt * (rp.r * spsi * dpsidphi * rp.pr + cpsi * rp.pphi) / A;
		return coord::PosMomCar(x, y, 0, px, py, 0);
	}
	coord::PosMomSph xyPointTrans::xp2rp(const coord::PosMomCar xp) const {
		double x2 = pow_2(xp.x), y2 = pow_2(xp.y), B = x2 + y2 - Delta2;
		double r2 = .5 * (B + sqrt(B * B + 4 * x2 * Delta2)), r = sqrt(r2);
		double cpsi = xp.x / r, rt = sqrt(r2 + Delta2);
		double spsi = xp.y > 0 ? sqrt(1 - cpsi * cpsi) : -sqrt(1 - cpsi * cpsi);
		double pr = cpsi * xp.px + r * spsi / rt * xp.py;
		double psi = atan2(spsi, cpsi), dpsidphi, phi = PC.Phi(psi, dpsidphi);
		double pphi = (-r * spsi * xp.px + rt * cpsi * xp.py) * dpsidphi;
		return coord::PosMomSph(r, .5 * M_PI, phi, pr, 0, pphi);
	}

	EXP ToyMap interpToyMap(double x, const ToyMap& TM0, const ToyMap& TM1) {
		const double xp = 1 - x;
		return ToyMap(interpIsochrone(x, TM0.Is, TM1.Is),
			interpPointTrans(x, TM0.PT, TM1.PT));
	}
	//Derivs of RzpR.. wrt parameters Delta, Js, b
	coord::PosMomCyl ToyMap::from_aaT(const ActionAngles& aaT, coord::PosMomCyl* dRzdPs) const {
		coord::PosMomSph drdJs, drdb;
		coord::PosMomSph rp(Is.aa2pq(aaT, drdJs, drdb));
		math::Matrix<double> dRzdrt(4, 4);
		coord::PosMomCyl Rz(PT.Sph2Cyl(rp, dRzdrt, &dRzdPs[0]));
		dRzdPs[1].R = dRzdrt(0, 0) * drdJs.r + dRzdrt(0, 1) * drdJs.theta
			+ dRzdrt(0, 2) * drdJs.pr + dRzdrt(0, 3) * drdJs.ptheta;
		dRzdPs[1].z = dRzdrt(1, 0) * drdJs.r + dRzdrt(1, 1) * drdJs.theta
			+ dRzdrt(1, 2) * drdJs.pr + dRzdrt(1, 3) * drdJs.ptheta;
		dRzdPs[1].pR = dRzdrt(2, 0) * drdJs.r + dRzdrt(2, 1) * drdJs.theta
			+ dRzdrt(2, 2) * drdJs.pr + dRzdrt(2, 3) * drdJs.ptheta;
		dRzdPs[1].pz = dRzdrt(3, 0) * drdJs.r + dRzdrt(3, 1) * drdJs.theta
			+ dRzdrt(3, 2) * drdJs.pr + dRzdrt(3, 3) * drdJs.ptheta;
		dRzdPs[2].R = dRzdrt(0, 0) * drdb.r + dRzdrt(0, 1) * drdb.theta
			+ dRzdrt(0, 2) * drdb.pr + dRzdrt(0, 3) * drdb.ptheta;
		dRzdPs[2].z = dRzdrt(1, 0) * drdb.r + dRzdrt(1, 1) * drdb.theta
			+ dRzdrt(1, 2) * drdb.pr + dRzdrt(1, 3) * drdb.ptheta;
		dRzdPs[2].pR = dRzdrt(2, 0) * drdb.r + dRzdrt(2, 1) * drdb.theta
			+ dRzdrt(2, 2) * drdb.pr + dRzdrt(2, 3) * drdb.ptheta;
		dRzdPs[2].pz = dRzdrt(3, 0) * drdb.r + dRzdrt(3, 1) * drdb.theta
			+ dRzdrt(3, 2) * drdb.pr + dRzdrt(3, 3) * drdb.ptheta;
		return Rz;
	}
	//Derivs of RzpR.. wrt actions
	coord::PosMomCyl ToyMap::from_aaT(const ActionAngles& aaT, DerivAct<coord::Cyl>& dRzdJ) const {
		DerivAct<coord::Sph> drdJ;
		coord::PosMomSph rp(Is.aa2pq(aaT, NULL, &drdJ));
		math::Matrix<double> dRzdrt(4, 4);
		coord::PosMomCyl Rz(PT.Sph2Cyl(rp, dRzdrt));
		dRzdJ.dbyJr.R = dRzdrt(0, 0) * drdJ.dbyJr.r + dRzdrt(0, 1) * drdJ.dbyJr.theta
			+ dRzdrt(0, 2) * drdJ.dbyJr.pr + dRzdrt(0, 3) * drdJ.dbyJr.ptheta;
		dRzdJ.dbyJr.z = dRzdrt(1, 0) * drdJ.dbyJr.r + dRzdrt(1, 1) * drdJ.dbyJr.theta
			+ dRzdrt(1, 2) * drdJ.dbyJr.pr + dRzdrt(1, 3) * drdJ.dbyJr.ptheta;
		dRzdJ.dbyJr.phi = drdJ.dbyJr.phi;
		dRzdJ.dbyJr.pR = dRzdrt(2, 0) * drdJ.dbyJr.r + dRzdrt(2, 1) * drdJ.dbyJr.theta
			+ dRzdrt(2, 2) * drdJ.dbyJr.pr + dRzdrt(2, 3) * drdJ.dbyJr.ptheta;
		dRzdJ.dbyJr.pz = dRzdrt(3, 0) * drdJ.dbyJr.r + dRzdrt(3, 1) * drdJ.dbyJr.theta
			+ dRzdrt(3, 2) * drdJ.dbyJr.pr + dRzdrt(3, 3) * drdJ.dbyJr.ptheta;
		dRzdJ.dbyJr.pphi = 0;
		dRzdJ.dbyJz.R = dRzdrt(0, 0) * drdJ.dbyJz.r + dRzdrt(0, 1) * drdJ.dbyJz.theta
			+ dRzdrt(0, 2) * drdJ.dbyJz.pr + dRzdrt(0, 3) * drdJ.dbyJz.ptheta;
		dRzdJ.dbyJz.z = dRzdrt(1, 0) * drdJ.dbyJz.r + dRzdrt(1, 1) * drdJ.dbyJz.theta
			+ dRzdrt(1, 2) * drdJ.dbyJz.pr + dRzdrt(1, 3) * drdJ.dbyJz.ptheta;
		dRzdJ.dbyJz.phi = drdJ.dbyJz.phi;
		dRzdJ.dbyJz.pR = dRzdrt(2, 0) * drdJ.dbyJz.r + dRzdrt(2, 1) * drdJ.dbyJz.theta
			+ dRzdrt(2, 2) * drdJ.dbyJz.pr + dRzdrt(2, 3) * drdJ.dbyJz.ptheta;
		dRzdJ.dbyJz.pz = dRzdrt(3, 0) * drdJ.dbyJz.r + dRzdrt(3, 1) * drdJ.dbyJz.theta
			+ dRzdrt(3, 2) * drdJ.dbyJz.pr + dRzdrt(3, 3) * drdJ.dbyJz.ptheta;
		dRzdJ.dbyJphi.R = dRzdrt(0, 0) * drdJ.dbyJphi.r + dRzdrt(0, 1) * drdJ.dbyJphi.theta
			+ dRzdrt(0, 2) * drdJ.dbyJphi.pr + dRzdrt(0, 3) * drdJ.dbyJphi.ptheta;
		dRzdJ.dbyJphi.z = dRzdrt(1, 0) * drdJ.dbyJphi.r + dRzdrt(1, 1) * drdJ.dbyJphi.theta
			+ dRzdrt(1, 2) * drdJ.dbyJphi.pr + dRzdrt(1, 3) * drdJ.dbyJphi.ptheta;
		dRzdJ.dbyJphi.phi = drdJ.dbyJphi.phi;
		dRzdJ.dbyJphi.pR = dRzdrt(2, 0) * drdJ.dbyJphi.r + dRzdrt(2, 1) * drdJ.dbyJphi.theta
			+ dRzdrt(2, 2) * drdJ.dbyJphi.pr + dRzdrt(2, 3) * drdJ.dbyJphi.ptheta;
		dRzdJ.dbyJphi.pz = dRzdrt(3, 0) * drdJ.dbyJphi.r + dRzdrt(3, 1) * drdJ.dbyJphi.theta
			+ dRzdrt(3, 2) * drdJ.dbyJphi.pr + dRzdrt(3, 3) * drdJ.dbyJphi.ptheta;
		dRzdJ.dbyJphi.pphi = 1;
		return Rz;
	}
	//Derivs of RzpR.. wrt angles
	coord::PosMomCyl ToyMap::from_aaT(const ActionAngles& aaT, DerivAng<coord::Cyl>& dRzdT) const {
		DerivAng<coord::Sph> drdT;
		coord::PosMomSph rp(Is.aa2pq(aaT, NULL, NULL, &drdT));
		math::Matrix<double> dRzdrt(4, 4);
		coord::PosMomCyl Rz(PT.Sph2Cyl(rp, dRzdrt));
		dRzdT.dbythetar.R = dRzdrt(0, 0) * drdT.dbythetar.r + dRzdrt(0, 1) * drdT.dbythetar.theta
			+ dRzdrt(0, 2) * drdT.dbythetar.pr + dRzdrt(0, 3) * drdT.dbythetar.ptheta;
		dRzdT.dbythetar.z = dRzdrt(1, 0) * drdT.dbythetar.r + dRzdrt(1, 1) * drdT.dbythetar.theta
			+ dRzdrt(1, 2) * drdT.dbythetar.pr + dRzdrt(1, 3) * drdT.dbythetar.ptheta;
		dRzdT.dbythetar.phi = drdT.dbythetar.phi;
		dRzdT.dbythetar.pR = dRzdrt(2, 0) * drdT.dbythetar.r + dRzdrt(2, 1) * drdT.dbythetar.theta
			+ dRzdrt(2, 2) * drdT.dbythetar.pr + dRzdrt(2, 3) * drdT.dbythetar.ptheta;
		dRzdT.dbythetar.pz = dRzdrt(3, 0) * drdT.dbythetar.r + dRzdrt(3, 1) * drdT.dbythetar.theta
			+ dRzdrt(3, 2) * drdT.dbythetar.pr + dRzdrt(3, 3) * drdT.dbythetar.ptheta;
		dRzdT.dbythetar.pphi = 0;
		dRzdT.dbythetaz.R = dRzdrt(0, 0) * drdT.dbythetaz.r + dRzdrt(0, 1) * drdT.dbythetaz.theta
			+ dRzdrt(0, 2) * drdT.dbythetaz.pr + dRzdrt(0, 3) * drdT.dbythetaz.ptheta;
		dRzdT.dbythetaz.z = dRzdrt(1, 0) * drdT.dbythetaz.r + dRzdrt(1, 1) * drdT.dbythetaz.theta
			+ dRzdrt(1, 2) * drdT.dbythetaz.pr + dRzdrt(1, 3) * drdT.dbythetaz.ptheta;
		dRzdT.dbythetaz.phi = drdT.dbythetaz.phi;
		dRzdT.dbythetaz.pR = dRzdrt(2, 0) * drdT.dbythetaz.r + dRzdrt(2, 1) * drdT.dbythetaz.theta
			+ dRzdrt(2, 2) * drdT.dbythetaz.pr + dRzdrt(2, 3) * drdT.dbythetaz.ptheta;
		dRzdT.dbythetaz.pz = dRzdrt(3, 0) * drdT.dbythetaz.r + dRzdrt(3, 1) * drdT.dbythetaz.theta
			+ dRzdrt(3, 2) * drdT.dbythetaz.pr + dRzdrt(3, 3) * drdT.dbythetaz.ptheta;
		dRzdT.dbythetaz.pphi = 0;
		dRzdT.dbythetaphi.R = dRzdrt(0, 0) * drdT.dbythetaphi.r + dRzdrt(0, 1) * drdT.dbythetaphi.theta
			+ dRzdrt(0, 2) * drdT.dbythetaphi.pr + dRzdrt(0, 3) * drdT.dbythetaphi.ptheta;
		dRzdT.dbythetaphi.z = dRzdrt(1, 0) * drdT.dbythetaphi.r + dRzdrt(1, 1) * drdT.dbythetaphi.theta
			+ dRzdrt(1, 2) * drdT.dbythetaphi.pr + dRzdrt(1, 3) * drdT.dbythetaphi.ptheta;
		dRzdT.dbythetaphi.phi = 1;
		dRzdT.dbythetaphi.pR = dRzdrt(2, 0) * drdT.dbythetaphi.r + dRzdrt(2, 1) * drdT.dbythetaphi.theta
			+ dRzdrt(2, 2) * drdT.dbythetaphi.pr + dRzdrt(2, 3) * drdT.dbythetaphi.ptheta;
		dRzdT.dbythetaphi.pz = dRzdrt(3, 0) * drdT.dbythetaphi.r + dRzdrt(3, 1) * drdT.dbythetaphi.theta
			+ dRzdrt(3, 2) * drdT.dbythetaphi.pr + dRzdrt(3, 3) * drdT.dbythetaphi.ptheta;
		dRzdT.dbythetaphi.pphi = 0;
		return Rz;
	}

	//interpolate between 2 tori
	EXP Torus interpTorus(const double x, const Torus& T0, const Torus& T1) {
		if (x == 1) return T0;
		if (x == 0) return T1;
		const double xp = 1 - x;
		Actions J = T0.J * x + T1.J * xp;
		Frequencies freqs = T0.freqs * x + T1.freqs * xp;
		GenFnc newGF(interpGenFnc(x, T0.GF, T1.GF));
		ToyMap newTM(interpToyMap(x, T0.TM, T1.TM));
		double newE = T0.E * x + T1.E * xp;
		return Torus(J, freqs, newGF, newTM, newE);
	}
	EXP Torus interpTorus(const double x, std::vector<double>& xs, std::vector<Torus>& Tgrid) {
		int bot, top;
		double f = bot_top(x,xs,bot,top);
		return interpTorus(f, Tgrid[bot], Tgrid[top]);
	}
	int TorusGrid1::botX(const double x) const{
		double bot=0, top=nx-1;
		while(top-bot>1){
			int n=.5*(bot+top);
			if(x>xs[n]) bot=n; else top=n;
		}
		return bot;
	}

	Torus TorusGrid1::T(double x) const{
		int botx;
		if(x < xs.front()) botx = 0;
		else if(x > xs.back()) botx = nx-2;
		else botx = botX(x);
		double fx = (x - xs[botx]) / (xs[botx+1] - xs[botx]);
		return interpTorus(fx,Ts[botx],  Ts[botx+1]);
	}
		
	int TorusGrid3::botX(const double x) const{
		double bot=0, top=nx-1;
		while(top-bot>1){
			int n=.5*(bot+top);
			if(x>xs[n]) bot=n; else top=n;
		}
		return bot;
	}
	int TorusGrid3::botY(const double x) const{
		double bot=0, top=ny-1;
		while(top-bot>1){
			int n=.5*(bot+top);
			if(x>ys[n]) bot=n; else top=n;
		}
		return bot;
	}
	int TorusGrid3::botZ(const double x) const{
		double bot=0, top=nz-1;
		while(top-bot>1){
			int n=.5*(bot+top);
			if(x>zs[n]) bot=n; else top=n;
		}
		return bot;
	}
	Torus TorusGrid3::T(const double x, const double y, const double z) const{
		int botx, boty, botz;
		if(x < xs.front()) botx = 0;
		else if(x > xs.back()) botx = nx-2;
		else botx = botX(x);
		if(y < ys.front()) boty = 0;
		else if(y > ys.back()) boty = ny-2;
		else boty = botY(y);
		if(z < zs.front()) botz = 0;
		else if(z > zs.back()) botz = nz-2;
		else botz = botZ(z);
		double fx = (x - xs[botx]) / (xs[botx+1] - xs[botx]);
		double fy = (y - ys[boty]) / (ys[boty+1] - ys[boty]);
		double fz = (z - zs[botz]) / (zs[botz+1] - zs[botz]);
		Torus T1(interpTorus(fx,Tn(botx,boty  ,botz),  Tn(botx+1,boty,  botz)));
		Torus T2(interpTorus(fx,Tn(botx,boty+1,botz),  Tn(botx+1,boty+1,botz)));
		Torus T3(interpTorus(fx,Tn(botx,boty  ,botz+1),Tn(botx+1,boty,  botz+1)));
		Torus T4(interpTorus(fx,Tn(botx,boty+1,botz+1),Tn(botx+1,boty+1,botz+1)));
		Torus T5(interpTorus(fy,T1,T2));
		Torus T6(interpTorus(fy,T3,T4));
		return interpTorus(fz,T5,T6);
	}
	EXP PerturbingHamiltonian interpPerturbingHamiltonian(double x,
		const PerturbingHamiltonian& H0, const PerturbingHamiltonian& H1){
		if (x==1) return H0; 
		if (x==0) return H1;
		double xp = 1-x;
		PerturbingHamiltonian H(H0), H2(H1);
		for (int i = 0; i < H0.indices.size(); i++) {//run over H0 indices
			int mr = H0.indices[i].mr, mz = H0.indices[i].mz, mphi = H0.indices[i].mphi;
			std::vector<std::complex<double> >::iterator jt = H2.values.begin();
			for (GenFncIndices::iterator it = H2.indices.begin(); it != H2.indices.end();) {
				if (mr == (*it).mr && mz == (*it).mz && mphi == (*it).mphi) {//this index matches mine
					H.values[i] = x*H.values[i] + xp*(*jt);
					it = H2.indices.erase(it);
					jt = H2.values.erase(jt);
					break;//my values updated for this index
				}
				else {
					it++; jt++;//move on to the next term in H2
				}
			}
		}
		if (H2.indices.size() > 0) {//something unaccounted for
			for (int j = 0; j < H2.indices.size(); j++) {//add them into my list
				H.indices.push_back(H2.indices[j]);
				H.values.push_back(H2.values[j]);
			}
		}
		return H;
	}

	EXP eTorus interpeTorus(const double x, const eTorus& eT0, const eTorus& eT1) {
		if (x == 1) return eT0;
		if (x == 0) return eT1;
		return eTorus(interpTorus(x, Torus(eT0), Torus(eT1)),
			      interpPerturbingHamiltonian(x,eT0.pH,eT1.pH));
	}
	EXP eTorus interpeTorus(const double x, std::vector<double>& xs,
				std::vector<eTorus>& Tgrid) {
		int bot, top;
		double f = bot_top(x,xs,bot,top);
		return interpeTorus(f, Tgrid[bot], Tgrid[top]);
	}

	coord::PosMomCyl Torus::from_toy(const Angles& thetaT) const {
		ActionAngles aaT(GF.toyJ(J, thetaT), thetaT);
		return TM.from_aaT(aaT);
	}
	coord::PosMomCyl Torus::from_true(const Angles& theta) const {//input true angles
		ActionAngles aaT(GF.true2toy(ActionAngles(J, theta)));//toy AAs computed from true
		return TM.from_aaT(aaT);
	}
	coord::PosCyl Torus::new_PosDerivJ(const Angles& thetaT,
		DerivAct<coord::Cyl>& dRJ) const {
		ActionAngles aaT(GF.toyJ(J, thetaT), thetaT);
		// Compute derivs wrt JT
		DerivAct<coord::Cyl> dRJT;
		coord::PosMomCyl Rz(TM.from_aaT(aaT, dRJT));
		//dtheta_i/dthetaT_j=dJT_j/dJ_i
		math::Matrix<double> dthetadthetaT(3, 3);
		GF.dtbydtT_Jacobian(thetaT, dthetadthetaT);

		math::Matrix<double> dxdJT(3, 3);
		dxdJT(0, 0) = dRJT.dbyJr.R; dxdJT(1, 0) = dRJT.dbyJr.z; dxdJT(2, 0) = dRJT.dbyJr.phi;
		dxdJT(0, 1) = dRJT.dbyJz.R; dxdJT(1, 1) = dRJT.dbyJz.z; dxdJT(2, 1) = dRJT.dbyJz.phi;
		dxdJT(0, 2) = dRJT.dbyJphi.R; dxdJT(1, 2) = dRJT.dbyJphi.z; dxdJT(2, 2) = dRJT.dbyJphi.phi;
		//dJTdJ_{ij} = dthetadthetaT_{ji}
		math::Matrix<double> Mat3(3, 3);
		math::blas_dgemm(math::CblasNoTrans, math::CblasTrans, 1.0, dxdJT, dthetadthetaT, 0.0, Mat3);

		math::Matrix<double> dpdJT(3, 3);
		dpdJT(0, 0) = dRJT.dbyJr.pR; dpdJT(1, 0) = dRJT.dbyJr.pz; dpdJT(2, 0) = dRJT.dbyJr.pphi;
		dpdJT(0, 1) = dRJT.dbyJz.pR; dpdJT(1, 1) = dRJT.dbyJz.pz; dpdJT(2, 1) = dRJT.dbyJz.pphi;
		dpdJT(0, 2) = dRJT.dbyJphi.pR; dpdJT(1, 2) = dRJT.dbyJphi.pz; dpdJT(2, 2) = dRJT.dbyJphi.pphi;

		math::Matrix<double> Mat4(3, 3);
		math::blas_dgemm(math::CblasNoTrans, math::CblasTrans, 1.0, dpdJT, dthetadthetaT, 0.0, Mat4);

		dRJ.dbyJr.R = Mat3(0, 0); dRJ.dbyJr.z = Mat3(1, 0); dRJ.dbyJr.phi = Mat3(2, 0);
		dRJ.dbyJz.R = Mat3(0, 1); dRJ.dbyJz.z = Mat3(1, 1); dRJ.dbyJz.phi = Mat3(2, 1);
		dRJ.dbyJphi.R = Mat3(0, 2); dRJ.dbyJphi.z = Mat3(1, 2); dRJ.dbyJphi.phi = Mat3(2, 2);
		dRJ.dbyJr.pR = Mat4(0, 0); dRJ.dbyJr.pz = Mat4(1, 0); dRJ.dbyJr.pphi = Mat4(2, 0);
		dRJ.dbyJz.pR = Mat4(0, 1); dRJ.dbyJz.pz = Mat4(1, 1); dRJ.dbyJz.pphi = Mat4(2, 1);
		dRJ.dbyJphi.pR = Mat4(0, 2); dRJ.dbyJphi.pz = Mat4(1, 2); dRJ.dbyJphi.pphi = Mat4(2, 2);
		return Rz;
	}
	/* Position from toy angle plus dR/dthetaT and dpR/dthetaT at fixed J (which causes JT
	 * to vary with thetaT)
	*/
	coord::PosCyl Torus::new_PosDerivs(const Angles& thetaT,
		DerivAngCyl& dRtT, double* det) const {
		ActionAngles aaT(GF.toyJ(J, thetaT), thetaT);
		coord::PosMomCyl Rz(TM.from_aaT(aaT, dRtT));
		// Compute derivs wrt JT
		DerivAct<coord::Cyl> dRJT;
		TM.from_aaT(aaT, dRJT);
		//dtheta_i/dthetaT_j=dJT_j/dJ_i
		math::Matrix<double> dthetadthetaT(3, 3);
		GF.dtbydtT_Jacobian(thetaT, dthetadthetaT);

		DerivAng<coord::Cyl> dJA = GF.dJdt(thetaT); //dJT/dthetaT

		dRtT.dbythetar.R += dRJT.dbyJr.R * dJA.dbythetar.R + dRJT.dbyJz.R * dJA.dbythetar.z;
		dRtT.dbythetaz.R += dRJT.dbyJr.R * dJA.dbythetaz.R + dRJT.dbyJz.R * dJA.dbythetaz.z;
		dRtT.dbythetar.z += dRJT.dbyJr.z * dJA.dbythetar.R + dRJT.dbyJz.z * dJA.dbythetar.z;
		dRtT.dbythetaz.z += dRJT.dbyJr.z * dJA.dbythetaz.R + dRJT.dbyJz.z * dJA.dbythetaz.z;
		dRtT.dbythetar.phi += dRJT.dbyJr.phi * dJA.dbythetar.R
			+ dRJT.dbyJz.phi * dJA.dbythetar.z
			+ dRJT.dbyJphi.phi * dJA.dbythetar.phi;
		dRtT.dbythetaz.phi += dRJT.dbyJr.phi * dJA.dbythetaz.R
			+ dRJT.dbyJz.phi * dJA.dbythetaz.z
			+ dRJT.dbyJphi.phi * dJA.dbythetaz.phi;;

		dRtT.dbythetar.pR += dRJT.dbyJr.pR * dJA.dbythetar.R + dRJT.dbyJz.pR * dJA.dbythetar.z;
		dRtT.dbythetaz.pR += dRJT.dbyJr.pR * dJA.dbythetaz.R + dRJT.dbyJz.pR * dJA.dbythetaz.z;
		dRtT.dbythetar.pz += dRJT.dbyJr.pz * dJA.dbythetar.R + dRJT.dbyJz.pz * dJA.dbythetar.z;
		dRtT.dbythetaz.pz += dRJT.dbyJr.pz * dJA.dbythetaz.R + dRJT.dbyJz.pz * dJA.dbythetaz.z;
		dRtT.dbythetar.pphi = dJA.dbythetar.phi;
		dRtT.dbythetaz.pphi = dJA.dbythetaz.phi;
		dRtT.dbythetaphi.pphi = dJA.dbythetaphi.phi;

		if (det) {//assume only index.mphi=0 non-zero
			(*det) = dRtT.dbythetar.R * dRtT.dbythetaz.z - dRtT.dbythetar.z * dRtT.dbythetaz.R;
		}
		return Rz;
	}

	void Torus::zSoS(std::vector<double>& Rs, std::vector<double>& vRs, const int N,
		double& Rmin, double& Rmax, double& Vmax, const double z0) const {
		const double tol = 1e-5;
		Rmin = 1e10, Rmax = 0, Vmax = 0;
		for (int i = 0; i < N; i++) {
			double thetaT_r = 2 * M_PI / (double)N * (-N / 2 + i);
			CrossingFinder CF(this, thetaT_r, z0);
			double thetaT_z, dtheta = .1, th_min = -.5 * M_PI,
				th_max = th_min + dtheta, z_min, z_max;
			CF.evalDeriv(th_min, &z_min); CF.evalDeriv(th_max, &z_max);
			while (z_min * z_max > 0) {// plod round looking for sign change
				th_min = th_max; z_min = z_max; th_max += dtheta;
				CF.evalDeriv(th_max, &z_max);
			}
			thetaT_z = math::findRoot(CF, th_min, th_max, tol);
			coord::PosMomCyl Rz(from_toy(Angles(thetaT_r, thetaT_z, 0)));
			if (Rz.pz < -.0001) {// keep going round
				do {
					th_min = th_max; z_min = z_max;
					th_max += dtheta; CF.evalDeriv(th_max, &z_max);
				} while (z_min * z_max > 0);
				thetaT_z = math::findRoot(CF, th_min, th_max, tol) - 2 * M_PI;
				Rz = from_toy(Angles(thetaT_r, thetaT_z, 0));
			}
			Rs.push_back(Rz.R); vRs.push_back(Rz.pR);
			Rmin = fmin(Rmin, Rz.R); Rmax = fmax(Rmax, Rz.R); Vmax = fmax(Vmax, fabs(Rz.pR));
		}
		Rs.push_back(Rs[0]); vRs.push_back(vRs[0]);
	}
	void Torus::rSoS(std::vector<double>& thetas, std::vector<double>& pthetas, const double rbar,
		const int N, double& thetamax, double& pmax) const {
		const double tol = 1e-5;
		thetamax = 0, pmax = 0;
		for (int i = 0; i < N; i++) {
			double thetaT_z = 2 * M_PI / (double)N * (-N / 2 + i);
			r_crossingFinder RCF(this, rbar, thetaT_z);
			double thetaT_r, dtheta = .1, th_min = 0,
				th_max = th_min + dtheta, dr_min, dr_max;
			RCF.evalDeriv(th_min, &dr_min); RCF.evalDeriv(th_max, &dr_max);
			while (dr_min * dr_max > 0 && th_min < 2 * M_PI) {// plod round looking for sign change
				th_min = th_max; dr_min = dr_max; th_max += dtheta;
				RCF.evalDeriv(th_max, &dr_max);
			}
			if (th_min >= 2 * M_PI) continue;
			thetaT_r = math::findRoot(RCF, th_min, th_max, tol);
			coord::PosMomSph rtheta(coord::toPosMomSph(from_toy(Angles(thetaT_r, thetaT_z, 0))));
			if (rtheta.pr < -.0001) {// keep going round
				do {
					th_min = th_max; dr_min = dr_max;
					th_max += dtheta; RCF.evalDeriv(th_max, &dr_max);
				} while (dr_min * dr_max > 0 && th_min < 2 * M_PI);
				thetaT_z = math::findRoot(RCF, th_min, th_max, tol) - 2 * M_PI;
				rtheta = coord::toPosMomSph(from_toy(Angles(thetaT_r, thetaT_z, 0)));
			}
			thetas.push_back(rtheta.theta); pthetas.push_back(rtheta.ptheta);
			thetamax = fmax(thetamax, rtheta.theta); pmax = fmax(pmax, fabs(rtheta.ptheta));
		}
		thetas.push_back(thetas[0]); pthetas.push_back(pthetas[0]);
	}
	void Torus::SoSthetaz(std::vector<double>& X, std::vector<double>& pX,
		const double thetaz,
		const int N, double& Xmax, double& pXmax) const {
		X.clear(); pX.clear();
		Xmax = 0; pXmax = 0;
		for (int i = 0; i < N; i++) {
			double thetar = i * 2 * M_PI / (double)(N - 1);
			Angles thetaT(thetar, thetaz, 0);
			coord::PosMomCyl Rz(from_toy(thetaT));
			ActionAngles aaT(TM.pq2aa(Rz));
			X.push_back(sqrt(aaT.Jr) * cos(aaT.thetar));
			pX.push_back(sqrt(aaT.Jr) * sin(aaT.thetar));
			Xmax = fmax(Xmax, fabs(X.back()));
			pXmax = fmax(pXmax, fabs(pX.back()));
		}
	}
	void Torus::SoSthetar(std::vector<double>& X, std::vector<double>& pX,
		const double thetar,
		const int N, double& Xmax, double& pXmax) const {
		X.clear(); pX.clear();
		Xmax = 0; pXmax = 0;
		for (int i = 0; i < N; i++) {
			double thetaz = i * 2 * M_PI / (double)(N - 1);
			Angles thetaT(thetar, thetaz, 0);
			coord::PosMomCyl Rz(from_toy(thetaT));
			ActionAngles aaT(TM.pq2aa(Rz));
			X.push_back(sqrt(aaT.Jz) * cos(aaT.thetaz));
			pX.push_back(sqrt(aaT.Jz) * sin(aaT.thetaz));
			Xmax = fmax(Xmax, fabs(X.back()));
			pXmax = fmax(pXmax, fabs(pX.back()));
		}
	}

	std::vector<std::pair<coord::PosVelCyl, double> > Torus::orbit(const Angles& theta0, double dt, double T) const {
		std::vector<std::pair<coord::PosVelCyl, double> > traj;
		double t = 0;
		while (t < T) {
			Angles theta(theta0 + (Angles)(freqs * t));
//			Angles theta(theta0.thetar + freqs.Omegar * t, theta0.thetaz + freqs.Omegaz * t,
//				theta0.thetaphi + freqs.Omegaphi * t);
			traj.push_back(std::pair<coord::PosVelCyl, double>(coord::toPosVelCyl(from_true(theta)), t));
			t += dt;
		}
		return traj;
	}

	// Returns true if (R,z,phi) is ever hit by the orbit, and false otherwise. If the 
	// torus passes through the point given, this happens four times, in each case
	// with a different velocity, but only two of these are independent:
	// theta_r -> -theta_r with theta_z -> Pi-theta_z leaves (R,z) and J^T fixed
	// but changes the sign of both velocities. |d(x,y,z)/d(theta_r,theta_z,theta_phi)|
	// is returned. The latter vanishes on the edge of the
	// orbit, such that its inverse, the density of the orbit, diverges there
	// (that's the reason why the density itself is not returned).

	bool Torus::containsPoint(const coord::PosCyl& p, std::vector<Angles>& As,
		std::vector<coord::VelCyl>& Vs,
		std::vector<double>& Jacobs,
		std::vector<DerivAngCyl>* dRdtheta,
		const double tol) const {
		coord::PosMomCyl peri(from_true(Angles(0, .5 * M_PI, 0))), apo(from_true(Angles(M_PI, 0, 0))),
			top(from_true(Angles(M_PI, .5 * M_PI, 0)));
		double Rmin = .95 * peri.R, Rmax = 1.05 * apo.R, zmax = 1.05 * fabs(top.z);
		if (p.R<Rmin || p.R>Rmax || fabs(p.z) > zmax) return false;
		locFinder LF(*this, p);
		double tolerance = 1e-8;
		double params[2] = { 1,1 }, result[2], dist, det;
		int maxNumIter = 200;
		coord::PosMomCyl P1; Angles A1, Atrue;
		DerivAngCyl dA;
		int done, kmax = 30, nfail = 0;

		while (As.size() < 4 && nfail < kmax) {
			double kount = 0;
			for (int k = 0; k < kmax; k++) {
				done = math::nonlinearMultiFit(LF, params, tolerance, maxNumIter, result);
				A1 = Angles(math::wrapAngle(result[0]), math::wrapAngle(result[1]), 0.0);
				P1 = from_toy(A1); P1.phi = p.phi;
				//from_toy set phi=0
				A1.thetaphi = TM.pq2aa(P1).thetaphi;

				ActionAngles aaT = TM.pq2aa(P1);
				dist = sqrt(pow_2(p.R - P1.R) + pow_2(p.z - P1.z));
				if (dist < 2 * tol) {
					Atrue = GF.trueA(A1);
					if (is_new(Atrue, As)) {
						break;
					}
				}
				params[0] += .3; params[0] = math::wrapAngle(params[0]);
				params[1] += .7;   params[1] = math::wrapAngle(params[1]);
				kount++;
				if (kount == kmax && As.size() == 0) return false;
			}
			if (kount < kmax) {
				math::Matrix<double> M(3, 3);//to hold dtheta/dthetaT
				As.push_back(Atrue); new_PosDerivs(A1, dA, &det);
				Vs.push_back(coord::VelCyl(P1.pR, P1.pz, P1.pphi / P1.R));
				Jacobs.push_back(fabs(det / GF.dtbydtT_Jacobian(A1, M)));
				if (dRdtheta) {
					new_assemble(dA, M); dRdtheta->push_back(dA);
				}
				A1.thetar = -A1.thetar; A1.thetaz = M_PI - A1.thetaz;
				P1 = from_toy(A1); P1.phi = p.phi;
				A1.thetaphi = TM.pq2aa(P1).thetaphi;
				Atrue = GF.trueA(A1);
				As.push_back(Atrue); new_PosDerivs(A1, dA, &det);
				Vs.push_back(coord::VelCyl(P1.pR, P1.pz, P1.pphi / P1.R));
				Jacobs.push_back(fabs(det / GF.dtbydtT_Jacobian(A1, M)));
				if (dRdtheta) {
					new_assemble(dA, M); dRdtheta->push_back(dA);
				}
			}
			else nfail++;
		}
		if (nfail >= kmax) printf("containsPoint error at Rz (%f %f) - %d angles \n",
			p.R, p.z, As.size());
		return true;
	}
	double Torus::density(const coord::PosCyl& Rz) const {
		std::vector<Angles> As; std::vector<coord::VelCyl> Vs;
		std::vector<double> Jacobs;
		//const double tol = 1e-6;
		if (!containsPoint(Rz, As, Vs, Jacobs)) return 0;
		double rho = 0;
		for (int i = 0; i < As.size(); i++)
			rho += 1 / Jacobs[i];
		return rho;
	}
	void Torus::write(FILE* ofile) const {
		fprintf(ofile, "%g %g %g %g %g %g %g %g %g %g\n",
			J.Jr, J.Jz, J.Jphi, freqs.Omegar, freqs.Omegaz, freqs.Omegaphi,
			E, TM.PT.cs.Delta, TM.Is.Js, TM.Is.b);
		GF.write(ofile);
	}
	void Torus::read(FILE* ifile) {
		double Delta, Js, b;
		fscanf_s(ifile, "%g %g %g %g %g %g %g %g %g %g\n",
			&J.Jr, &J.Jz, &J.Jphi, &freqs.Omegar, &freqs.Omegaz, &freqs.Omegaphi,
			&E, &Delta, &Js, &b);
		PointTrans newPT(Delta);
		Isochrone newIs(Js, b);
		TM = ToyMap(newIs, newPT);
		GF.read(ifile);
	}

	/* Gather all terms in the pH that are harmonics of the specified line
	 */
	std::vector<std::complex<double> > PerturbingHamiltonian::get_hn(const GenFncIndex& I,
		std::vector<float>& multiples) const {
		std::vector<std::complex<double> > Hs;
		printf("Looking for (%d %d %d)\n", I.mr, I.mz, I.mphi);
		for (int i = 0; i < indices.size(); i++) {
			GenFncIndex In(indices[i]);
			if ((I.mr == 0 && In.mr != 0) || (I.mr != 0 && In.mr == 0)) continue;
			if ((I.mz == 0 && In.mz != 0) || (I.mz != 0 && In.mz == 0)) continue;
			if ((I.mphi == 0 && In.mphi != 0) || (I.mphi != 0 && In.mphi == 0)) continue;
			//now either matching zero indices or both non-zero
			std::vector<double> Rats;
			if (I.mr != 0) Rats.push_back((double)I.mr / (double)In.mr);
			if (I.mz != 0) Rats.push_back((double)I.mz / (double)In.mz);
			if (I.mphi != 0) Rats.push_back((double)I.mphi / (double)In.mphi);
			if (Rats[1] != Rats[0]) continue;
			if ((Rats.size() == 3) && (Rats[2] != Rats[0])) continue;
			printf("%3d %2d %2d %g %f at rank %d\n", In.mr, In.mz, In.mphi, math::modulus(values[i]),
				math::arg(values[i]) / M_PI, i);
			Hs.push_back(values[i]);
			if (Rats[0] >= 1) multiples.push_back(Rats[0]);
			else multiples.push_back(1 / Rats[0]);
		}
		return Hs;
	}
	TMfitter::TMfitter(const potential::BasePotential& pot,
		std::vector<std::pair<coord::PosMomCyl, double> >& _traj,
		double _pphi) : traj(_traj), pphi(_pphi) {
		xmin = 1e6; ymax = 0;
		double ymin = 1e6, xmax = 0, pxmin = 1e6, pxmax = 0, pymin = 1e6, pymax = 0;
		for (int i = 1; i < traj.size(); i++) {//Find where crosses axes
			if (traj[i].first.R * traj[i - 1].first.R <= 0) {//x axis
				double f = traj[i].first.R / (traj[i].first.R - traj[i - 1].first.R);
				double y = (1 - f) * traj[i].first.z + f * traj[i - 1].first.z;
				double px = (1 - f) * traj[i].first.pR + f * traj[i - 1].first.pR;
				ymin = fmin(ymin, fabs(y)); ymax = fmax(ymax, fabs(y));
				pxmin = fmin(pxmin, fabs(px)); pxmax = fmax(pxmax, fabs(px));
			}
			if (traj[i].first.z * traj[i - 1].first.z <= 0) {//z axis
				double f = traj[i].first.z / (traj[i].first.z - traj[i - 1].first.z);
				double x = (1 - f) * traj[i].first.R + f * traj[i - 1].first.R;
				double py = (1 - f) * traj[i].first.pz + f * traj[i - 1].first.pz;
				xmin = fmin(xmin, fabs(x)); xmax = fmax(xmax, fabs(x));
				pymin = fmin(pymin, fabs(py)); pymax = fmax(pymax, fabs(py));
			}
		}
		//xbar, ybar estimated axes of underlying loop orbit
		xbar = .5 * (xmin + xmax); double ybar = .5 * (ymin + ymax);
		Delta2 = (ybar * ybar - xbar * xbar);
		double Phi; coord::GradCyl grad;
		pot.eval(coord::PosCyl(0, ymax, 0), &Phi, &grad);
		Frat = grad.dz;
		pot.eval(coord::PosCyl(xmin, 0, 0), &Phi, &grad);
		Frat /= grad.dR;
		double pxbar = pphi > 0 ? -.5 * (pxmin + pxmax) : .5 * (pxmin + pxmax),
			pybar = pphi > 0 ? .5 * (pymin + pymax) : -.5 * (pymin + pymax);
		aPT = 0.25 * pphi * (1 / (ybar * pybar) + 1 / (xbar * pxbar));
		bPT = 0.125 * pphi * (1 / (ybar * pybar) - 1 / (xbar * pxbar)) - 0.25;
		printf("xbar %f ybar %f ", xbar, ybar);
		printf("pxbar %f pybar %f\n", pxbar, pybar);
		printf("Delta2 %f Frat %f a %f b %f\n", Delta2, Frat, aPT, bPT);
	}
	//Function with root where we have the right force ratio
	double TMfitter::value(double Js) const {
		double g2 = pow((fabs(pphi) + sqrt(pphi * pphi + 4 * Js * Js)) / (2 * Js), 4) - 1;
		double bIso = xbar / sqrt(g2);
		double ap = sqrt(ymax * ymax + bIso * bIso), am = sqrt(xmin * xmin + bIso * bIso);
		//	printf("Js pphi g2 b %f %g %g %g\n",Js,pphi,g2,b);
		return Frat - (ymax * am * pow_2(bIso + am)) / (xmin * ap * pow_2(bIso + ap));
	}
	//Solve for Js and bIso
	std::vector<double> TMfitter::fitTM() const {
		const double Jsmin = .01, Jsmax = 10;
		double Js = math::findRoot(*this, Jsmin, Jsmax, 1e-5);
		double g2 = (pow((fabs(pphi) + sqrt(pphi * pphi + 4 * Js * Js)) / (2 * Js), 4) - 1);
		double bIso = xbar / sqrt(g2);
		std::vector<double> ans(5);
		ans[0] = sqrt(Delta2); ans[1] = aPT; ans[2] = bPT;
		ans[3] = Js; ans[4] = bIso;
		return ans;
	}

	TorusGenerator::TorusGenerator(const potential::BasePotential& _pot, const double _tol) :
		pot(_pot), defaultTol(_tol), invPhi0(1. / _pot.value(coord::PosCyl(0, 0, 0))), tmax(250) {
		std::vector<double> gridR = potential::createInterpolationGrid(pot, ACCURACY_INTERP2);
		int sizeL = gridR.size(), sizeXi = 20;
		printf("Preparing TorusGenerator...");
		std::vector<double> gridL(sizeL);
		std::vector<double> gridE(sizeL);
		std::vector<double> gridXi(sizeXi);
		for (int i = 0; i < sizeL; i++) {
			gridL[i] = gridR[i] * potential::v_circ(pot, gridR[i]);
			gridE[i] = E_circ(pot, gridL[i]);
		}
		for (int i = 0; i < sizeXi; i++)
			gridXi[i] = i / (double)(sizeXi - 1);//EV wld call this Xiscaled
		math::Matrix<double> grid2dD(sizeL, sizeXi);
		math::Matrix<double> grid2dR(sizeL, sizeXi);
		math::Matrix<double> grid2dV(sizeL, sizeXi);
		math::Matrix<double> grid2dRE(sizeL, sizeXi);
		createGridFocalDistance(pot, gridL, gridXi, gridE, grid2dD,
			grid2dR, grid2dV, grid2dRE);
		interpD = math::LinearInterpolator2d(gridL, gridXi, grid2dD);
		interpR = math::LinearInterpolator2d(gridL, gridXi, grid2dR);
		interpV = math::LinearInterpolator2d(gridL, gridXi, grid2dV);
		interpRE = math::LinearInterpolator2d(gridE, gridXi, grid2dRE);
		printf("done\n");
	}
	double TorusGenerator::Hamilton(const Torus& T, const potential::BasePotential* ePot, const Angles& theta)
	{
		coord::PosMomCyl Rz(T.from_true(theta));
		double H = .5 * (pow_2(Rz.pR) + pow_2(Rz.pz) + pow_2(Rz.pphi / Rz.R)) + pot.value(Rz);
		coord::PosCyl pos(Rz.R, Rz.z, Rz.phi);
		return ePot ? H + ePot->value(pos) : H;
	}
	PerturbingHamiltonian TorusGenerator::get_pH(const Torus& T, int nf, bool ifp,
		const potential::BasePotential* ePot) {//Fourier analyses H
		int nfr = nf, nfz = nf;
		int nfp = ePot ? nf / 4 : 1;
		double N = (nfp * nfr * nfz);
		double* h = new double[nfp * nfr * nfz];
		Angles thetas;
		double dtr = 2 * M_PI / (double)nfr, dtz = 2 * M_PI / (double)nfz, dtp = M_PI / (double)nfp;
#pragma omp parallel for schedule(dynamic) 
		for (int k = 0; k < nfp; k++) {
			thetas.thetaphi = k * dtp;//if !ePot sticks at 0
			for (int i = 0; i < nfr; i++) {
				int i1 = nfr - i;
				thetas.thetar = i * dtr;
				for (int j = 0; j < nfz; j++) {
					int j1 = nfz - j, j2 = nfz / 2 - j, j3 = nfz / 2 + j;
					thetas.thetaz = j * dtz;
					double a = Hamilton(T, ePot, thetas);
					h[nfr * nfz * k + nfz * i + j] = a;
					/*
					if(j3<nfz) h[nfr*nfz*k+i*nfz+j3]=a;//N-S symmetry
					if(i1!=i && i1<nfr){
						if(j1<nfz) h[nfr*nfz*k+i1*nfz+j1]=a;//time-reverse symmetry
						if(j2<nfz) h[nfr*nfz*k+i1*nfz+j2]=a;//both symmetries
					}*/
				}
			}
		}
		//need speq to return cpts at Nyquist vals m_phi, which we won't use 
		math::Matrix<double> speq(nfr, 2 * nfz);
		rlft3(h, speq, nfp, nfr, nfz, 1);
		std::vector<double> Hmods;
		GenFncIndices Hindices;
		int ntop = 0;
		double hmax = 0;
		for (int k = 0; k < nfp; k++) {//find largest perturbing terms
			for (int i = 0; i < nfr; i++) {
				for (int j = 0; j < nfz / 2; j++) {
					double s = sqrt(pow_2(h[nfr * nfz * k + nfz * i + 2 * j])
						+ pow_2(h[nfr * nfz * k + nfz * i + 2 * j + 1])) / N;
					if (ntop<tmax || s>hmax) {
						hmax = insertLine(ntop, tmax, s, GenFncIndex(i, j, k),
							Hmods, Hindices);
					}
				}
			}
		}
		ifp = false;
		std::vector<std::complex<double> > Hvalues;
		for (int i = 0; i < ntop; i++) {
			Hvalues.push_back(std::complex<double>(h[nfr * nfz * Hindices[i].mphi + nfz * Hindices[i].mr + 2 * Hindices[i].mz] / N,
				h[nfr * nfz * Hindices[i].mphi + nfz * Hindices[i].mr + 2 * Hindices[i].mz + 1] / N));
		}
		if (ifp) printf("Terms in perturbing H (*100)\n");
		for (int i = 0; i < ntop; i++) {
			//		printf("%d %d %d ",Hindices[i].mr, Hindices[i].mz, Hindices[i].mphi);
			if (Hindices[i].mphi > nfp / 2) Hindices[i].mphi -= nfp;
			if (Hindices[i].mr > nfr / 2)   Hindices[i].mr -= nfr;
			//		Hindices[i].mz*=2;
			//		if(Hindices[i].mz>nfz/2)   Hindices[i].mz   -= nfz;//won't happen!
			if (ifp) printf("(%3d %3d %3d) (%g %g)\n",
				Hindices[i].mr, Hindices[i].mz, Hindices[i].mphi,
				100 * math::modulus(Hvalues[i]), math::arg(Hvalues[i]));
		}
		delete[] h;
		return PerturbingHamiltonian(Hindices, Hvalues);
	}

	double TorusGenerator::getRsh(Actions& J) {
		const double L = J.Jr + J.Jz + fabs(J.Jphi);
		return interpR.value(L, (J.Jr + J.Jz) / L);
	}
	void writeClosed(std::vector<coord::PosVelCyl>& shell, Actions J, ToyMap& TM, ToyMap& baldTM) {
		FILE* ofile; fopen_s(&ofile, "writeClosed.dat", "w");
		int npt = 200;
		fprintf(ofile, "%d\n", npt);
		std::vector<coord::PosMomCyl> Rzr, Rzt;
		for (int i = 0; i < npt; i++) {//draw images of shell isochrone 
			double xr, yr, xt, yt, theta = 2 * M_PI * i / (double)npt;
			coord::PosMomCyl Rzr(TM.from_aaT(J, Angles(M_PI * theta, theta, 0)));
			xr = Rzr.R; yr = Rzr.pR;
			coord::PosMomCyl Rzt(TM.from_aaT(J, Angles(1, theta, 0)));
			xt = Rzt.z; yt = Rzt.pz;
			fprintf(ofile, "%f %f %f %f\n", xr, yr, xt, yt);
		};
		fprintf(ofile, "%zd\n", shell.size());
		for (int i = 0; i < shell.size(); i++) {
			fprintf(ofile, "%f %f %f %f\n",
				shell[i].R, shell[i].vR, shell[i].z, shell[i].vz);
		}
		fclose(ofile);
	}
	ToyMap TorusGenerator::chooseTM(GenFncFitSeries& GFFS0, std::vector<double>& params,
					const Actions& J, double& Jscale,
					double& freqScale, double& Rsh) const {
		const double L = fabs(J.Jphi) + J.Jz, Xi = J.Jz / L;
		const double Jtot = L + J.Jr;
		Jscale = J.Jr + J.Jz;
		double Delta = interpD.value(L, Xi);// Delta *= exp(-3 * J.Jr / J.Jz);
		double fac = exp(-pow(J.Jr / J.Jz, 1.5));//fac -> 1 for shell
		double E=1;
		if (E < 0) {
			Rsh = (1 - fac) * .8 * interpR.value(L, Xi) + fac * interpRE.value(E, Xi);
		} else {
			Rsh = interpR.value(L, Xi);
			//Rsh = (.9*(1 - fac) + .1* fac) * interpR.value(L, Xi);
		}
		freqScale = potential::v_circ(pot, Rsh) / Rsh; //frequency scale set
		//Now choose isochrone
			//For any Js b=Rsh/ISO.g(Js) f_apo_peri=ISO.f(b.Js,e) & where
			//this matches F_apo_peri in pot we pick Js and b
		Iso ISO(L, J.Jr);
		double Jsmax = 1.1 * Jtot, Jsmin = .1 * Jsmax;
		JsFinder JF(pot, ISO, Rsh);
		double val1, val2;
		JF.evalDeriv(Jsmin, &val1); JF.evalDeriv(Jsmax, &val2);
		while (val1 > 0) {
			Jsmin *= .75; JF.evalDeriv(Jsmin, &val1);
		}
		while (val2 < 0) {
			Jsmax *= 1.5; JF.evalDeriv(Jsmax, &val2);
		}
		const double relToler = 1e-4;
		double Js_iso = math::findRoot(JF, Jsmin, Jsmax, relToler);
		double b_iso = Rsh / ISO.g(Js_iso);
		Isochrone Is(Js_iso, b_iso);
		int Nn = 5;
		std::vector<double> params2(Nn, 0.0);
		fitvtmap vtmapf(Nn, GFFS0, J, Is, Delta, pot, params);
		double tolerance = 1e-9;//controls optimisation
		math::nonlinearMultiFit(vtmapf, &params2[0], tolerance, 20, &params2[0]);
		PointTrans PT(Delta, params2);
		return ToyMap(Is, PT);
	}
	Torus TorusGenerator::giveBaseTorus(const Actions& J, const ToyMap& TM) const {
		std::vector<double> params;
		int nrmax = 0, nzmax = 0;// nzmax must be even
		double Hbar, Hdisp = 1e20;
		GenFncIndices indices = makeGridIndices(nrmax, nzmax);
		GenFncFitFracs fracs;
		GenFncFitSeries GFFS(indices, fracs, J);
		torusFitter TF(J, pot, TM, GFFS);
		Frequencies freqs;
		GenFncDerivs dPdJ;
		Hdisp = TF.fitAngleMap(&params[0], Hbar, freqs, dPdJ);
		GenFncFracs Fracs;
		GenFnc G(indices, params, dPdJ, Fracs);
		return Torus(J, freqs, G, TM, Hbar);
	}
	
	Torus TorusGenerator::fitTorus(const Actions& J, const double tighten) const {
		int nrmax = 2, nzmax = 6;// nzmax must be even
		GenFncIndices indices = makeGridIndices(nrmax, nzmax);
		std::vector<double> params(indices.size(),0);
		GenFncFitFracs fracs;
		GenFncFitSeries GFFS0(indices, fracs, J);
		double Jscale, freqScale, Rsh;
		ToyMap TM(chooseTM(GFFS0, params, J, Jscale, freqScale, Rsh));
		double tolerance = 1e-9;//controls optimisation of the given Sn
		double tol = defaultTol * tighten;
		double Hbar, Hdisp = 1e20, Htarget = tol * freqScale * Jscale;
		bool converged = false;
		int Loop = 0, MaxLoop = 12, maxNumIter = 10;
		std::vector<double> rep;
		do {
			double Hdisp_old=Hdisp;
			GenFncFitSeries GFFS(indices, fracs, J);
			torusFitter TF(J, pot,TM, GFFS);
			//if(fracs.size()>0) TF.testit(params,GFFS.numParams()-1);
			if (std::isnan(Hdisp)) exit(0);
			//			printf("TF: %d angles, %d terms, %d fracs, %d parameters\n",
			//			       TF.numValues(),GFFS.numTerms(),GFFS.numFracs(),TF.numVars());
			try {
				NANbar = 0;
				int numIter = math::nonlinearMultiFit(TF, &params[0], tolerance, maxNumIter, &params[0], &Hdisp);
				Hdisp = sqrt(Hdisp);
				rep.push_back(Hdisp);
				//if (NANbar > 0) printf("\ndH1 :%g NANfrac %f %f\n", Hdisp, NANfrac, NANbar);
				converged = (Hdisp < tol * freqScale * Jscale);
			}
			catch (std::exception& e) {
				std::cout << "Exception in fitTorus: " << e.what() << '\n';
			}
			if (converged || (Loop > 1 && fabs(1 - rep[Loop] / rep[Loop - 1]) < 1e-2)) break;
			Loop++;
			if (Loop < MaxLoop)
				indices = GFFS.expand(params, fracs);
		} while (Loop < MaxLoop);
		GenFncFitSeries GFFS(indices, fracs, J);
		if (!converged) {
			printf("\nfitTorus failed to converge: %e against target %e  ",
				Hdisp, Htarget);
			printf("%zd indices %zd fracs", indices.size(), fracs.size());
			printf("NANfracs: %f %f, resids:\n", NANfrac, NANbar);
			for (int i = 0; i < rep.size(); i++) printf("%7.2e ", rep[i]); printf("\n");
		}
		torusFitter TF(J, pot, TM, GFFS);
		Frequencies freqs;
		GenFncDerivs dPdJ;
		Hdisp = TF.fitAngleMap(&params[0], Hbar, freqs, dPdJ);
		GenFncFracs Fracs;
		for (unsigned int j = 0; j < fracs.size(); j++)
			Fracs.push_back(GenFncFrac(fracs[j],
				&params[indices.size() + 2 * j], &dPdJ[indices.size() + 2 * j]));
		params.resize(indices.size()); dPdJ.resize(indices.size());
		GenFnc G(indices, params, dPdJ, Fracs);
		return Torus(J, freqs, G, TM, Hbar);
	}

	Torus TorusGenerator::fitBaseTorus(const Actions& J, const double tighten) const {
		int nrmax = 2, nzmax = 6;// nzmax must be even
		GenFncIndices indices = makeGridIndices(nrmax, nzmax);
		std::vector<double> params(indices.size(),0);
		GenFncFitFracs fracs;
		GenFncFitSeries GFFS0(indices, fracs, J);
		double Jscale, freqScale, Rsh;
		ToyMap TM(chooseTM(GFFS0, params, J, Jscale, freqScale, Rsh));
		return giveBaseTorus(J, TM);
	}

	Torus TorusGenerator::fitFrom(const Actions& J, const Torus& T, const double tighten) const {
		GenFncIndices indices(T.GF.indices);
		std::vector<double> params(T.GF.values);
		GenFncFitFracs fracs;
		for (int i = 0; i < T.GF.fracs.size(); i++) {
			fracs.push_back(GenFncFitFrac(T.GF.fracs[i].mz, T.GF.fracs[i].mphi));
			params.push_back(T.GF.fracs[i].B);
			params.push_back(atanh(T.GF.fracs[i].b));
		}
		double Jscale = J.Jr + J.Jz, freqScale = T.freqs.Omegaz;
		ToyMap TM(T.TM);
		double tolerance = 1e-9;//controls optimisation of the given Sn
		double Hbar, Hdisp = 1e20;
		int maxNumIter = 10;
		double tol = defaultTol * tighten;
		std::vector<double> rep;
		GenFncFitSeries GFF(indices, fracs, J);
		torusFitter TF(J, pot, TM, GFF);
		try {
			int numIter = math::nonlinearMultiFit(TF, &params[0], tolerance, maxNumIter, &params[0], &Hdisp);
			Hdisp = sqrt(Hdisp);
			rep.push_back(Hdisp);
		}
		catch (std::exception& e) {
			std::cout << "Exception in fitFrom: " << e.what() << '\n';
		}
		Frequencies freqs;
		GenFncDerivs dPdJ;
		Hdisp = TF.fitAngleMap(&params[0], Hbar, freqs, dPdJ);
		GenFncFracs Fracs;
		for (unsigned int j = 0; j < fracs.size(); j++)
			Fracs.push_back(GenFncFrac(fracs[j],
				&params[indices.size() + 2 * j], &dPdJ[indices.size() + 2 * j]));
		params.resize(indices.size()); dPdJ.resize(indices.size());
		GenFnc G(indices, params, dPdJ, Fracs);
		return Torus(J, freqs, G, TM, Hbar);
	}
	eTorus TorusGenerator::fiteTorus(const Actions& J, const double tighten,
		const potential::BasePotential* ePhi) {
		Torus T = fitTorus(J, tighten);
		int nf = 128;
		PerturbingHamiltonian pH(get_pH(T, nf, true, ePhi));
		return eTorus(T, pH);
	}
	eTorus TorusGenerator::fiteTorus(const Actions& J, const potential::BasePotential* ePhi) {
		Torus T = fitTorus(J);
		int nf = 128;
		PerturbingHamiltonian pH(get_pH(T, nf, true, ePhi));
		return eTorus(T, pH);
	}
	double TorusGenerator::getDelta(double& L, double& Xi) {
		return interpD.value(L, Xi);
	}
	double TorusGenerator::getDelta(Actions& J) {
		double L = fabs(J.Jphi) + J.Jz, Xi = J.Jz / L;
		return interpD.value(L, Xi);
	}
	std::vector<Torus> TorusGenerator::constE(const double Jrmin, const Actions& Jstart, const int Nsteps) {
		double fac = exp(-log(Jstart.Jr / Jrmin) / (Nsteps - 1));
		Actions Jnext(Jstart);
		std::vector<Torus> Tgrd;
		Torus T(fitTorus(Jnext));
		Tgrd.push_back(T);
		for (int i = 1; i < Nsteps; i++) {
			double frat = Tgrd.back().freqs.Omegar / Tgrd.back().freqs.Omegaz;
			double dJr = (1 - fac) * Jnext.Jr;
			Jnext.Jr -= dJr; Jnext.Jz += dJr * frat;
			Torus T1(fitTorus(Jnext));
			Torus T2 = interpTorus(0.5, Tgrd.back(), T1);
			double frat1 = T2.freqs.Omegar / T2.freqs.Omegaz;
			Jnext.Jz = Jnext.Jz + dJr * frat - dJr * frat1;
			printf("Next Jr/Jz: %f Rsh: %f\n", Jnext.Jr / Jnext.Jz, getRsh(Jnext));
			Tgrd.push_back(fitTorus(Jnext));
		}
		return Tgrd;
	}

	EXP ActionAngles ActionFinderTG::actionAngles(const coord::PosVelCyl& point, Frequencies* freq) const {
		Torus T;
		ActionAngles aa(actionAnglesTorus(point, T));
		if (freq) *freq = T.freqs;
		return aa;
	}
	EXP ActionAngles ActionFinderTG::actionAnglesTorus(
		const coord::PosVelCyl& point, Torus& T) const {
		double E = potential::totalEnergy(*pot, point);
		if (E > 0) {
			printf("Energy is: %f and positive so no orbit\n", E);
			return ActionAngles(Actions(NAN, NAN, NAN), Angles(NAN, NAN, NAN));
		}
		const double tol = 1e-5;
		double phi0 = point.phi, last_diff = 1e6;
		while (fabs(phi0) > M_PI) phi0 += phi0 > M_PI ? -M_PI : M_PI;
		coord::PosMomCyl P0(point.R, point.z, phi0, point.vR, point.vz, point.vphi * point.R);
		Angles trueA;
		Actions J(AF->actions(point));
		if(std::isnan(J.Jr) || std::isnan(J.Jz)){
			printf("Fudge has returned a NaN: Jr,Jz,Jphi (%f %f %f)\n",
			       J.Jr,J.Jz,J.Jphi);
			exit(0);
		}
		T = TG.fitTorus(J);
		std::vector<double> p1 = T.TM.PT.paramsF;
		std::vector<double> dJt_old;
		int kount = 0;
		ActionAngles aaT = T.TM.pq2aa(P0);
		Angles tT(aaT);
		Actions JT(aaT);
		Actions JT2 = T.GF.toyJ(J, tT);
		double diff = 1e3;
		while (kount < 10) {
			if (kount > 0) 
				T = Torus(TG.fitFrom(J, T));
			Actions JT1 = T.GF.toyJ(J, tT);
			std::vector<double> df = { JT.Jr - JT1.Jr,JT.Jz - JT1.Jz,0.0 };
			math::Matrix<double> dthetadthetaT(3, 3);
			T.GF.dtbydtT_Jacobian(tT, dthetadthetaT);
			math::Matrix<double> Mat3(3, 3);
			math::LUDecomp LUM(dthetadthetaT);
			math::Matrix<double> inv = LUM.inverse(3);
			std::vector<double> dJt(3, 0.0);
			math::blas_dgemv(math::CblasTrans, 1.0, inv, df, 0.0, dJt);
			double old_diff = diff;
			diff = sqrt(pow_2(dJt[0]) + pow_2(dJt[1]));
			if(diff >= old_diff){//We've gone backwards!
				J.Jr -= dJt_old[0]; J.Jz -= dJt_old[1];
				break;
			}
			dJt_old = dJt;
			J.Jr += dJt[0]; J.Jz += dJt[1];
			printf("kount:%d (%f,%f) %f\n", kount, J.Jr, J.Jz,diff);
			if (sqrt(pow_2(dJt[0]) + pow_2(dJt[1])) < tol) {
				if (kount > 0) {
					trueA = T.GF.trueA(tT);
					break;
				}
			}
			kount++;
		}
		return ActionAngles(J, trueA);
	}
}//namespace
