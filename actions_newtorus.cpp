#include "actions_newtorus.h"
//#define TESTIT

/// accuracy parameter determining the spacing of the interpolation grid along the energy axis
static const double ACCURACY_INTERP2 = 1e-6;

namespace actions {

namespace{//internal

struct PVUVSph {
	double u,v,phi,pu,pv,pphi;
};
struct DerivActUVSph {
	PVUVSph dbyJr, dbyJz, dbyJphi;
};
//Two functions used in containsPoint
double angleDiff(Angles A1,Angles A2){
	double dr=A1.thetar-A2.thetar;
	while(dr>M_PI) dr-=2*M_PI;
	while(dr<-M_PI) dr+=2*M_PI;
	double dz=A1.thetaz-A2.thetaz;
	while(dz>M_PI) dz-=2*M_PI;
	while(dz<-M_PI) dz+=2*M_PI;
	double dp=A1.thetaphi-A2.thetaphi;
	while(dp>M_PI) dp-=2*M_PI;
	while(dp<-M_PI) dp+=2*M_PI;
	return sqrt(dr*dr+dz*dz+dp*dp);
}
	
bool is_new(Angles A1, std::vector<Angles> As){
	bool ok=true; const double tiny=1e-5;
	for(int i=0; i<As.size(); i++){
		double diff=angleDiff(A1,As[i]);
		ok = ok && diff>tiny;
	}
	return ok;
}			     
class locFinder : public math::IFunctionNdimDeriv {
	private:
		const BaseTorus T;
		const coord::PosCyl P0;
	public:
		locFinder(const BaseTorus _T,const coord::PosCyl& _P0) : T(_T), P0(_P0){}
		virtual void evalDeriv(const double params[], double* values, double* derivs=NULL) const {
			Angles A1(params[0],params[1],0);
			coord::PosCyl P1;  DerivAng<coord::Cyl> dA;
			if(derivs) {
				P1 = (T.PosDerivs(A1, dA));
				derivs[0]=dA.dbythetar.R; derivs[1]=dA.dbythetaz.R;
				derivs[2]=dA.dbythetar.z; derivs[3]=dA.dbythetaz.z;
			}
			if(values && !derivs)
				P1=T.from_toy(A1);
			if(values) {
				values[0]=P1.R-P0.R; values[1]=P1.z-P0.z;
				double dist=sqrt(pow_2(values[0])+pow_2(values[1]));
			}
		}
		virtual unsigned int numVars() const {
			return 2;
		}
		virtual unsigned int numValues() const {
			return 2;
		}
		void print(void) const{
			printf("LF: %f %f %f",T.TM.Js,T.TM.b,T.cs.Delta);
		}
};
				
double insertLine(int& ntop, const int tmax, double s, int i, int j, int k,
		  std::vector<double>& Hmods, GenFncIndices& Hindices){
	if(ntop==0){
		Hmods.push_back(s);
		Hindices.push_back(GenFncIndex(i,j,k));
		ntop++; return Hmods[0];
	}
	int l=0;
	while(l<ntop && s<=Hmods[l]) l++;
	if(l==ntop){//line isn't stronger than any previous line
		if(ntop==tmax){
			return Hmods[ntop-1];//no room for this line
		} else {//add line to end of list
			Hmods.push_back(s);
			Hindices.push_back(GenFncIndex(i,j,k));
			ntop++;
			return Hmods.back();
		}
	} else {//we should insert current line
		if(ntop<tmax){//move existing terms down
			Hmods.push_back(Hmods[ntop-1]);
			Hindices.push_back(Hindices[ntop-1]);
		}
		for(int m = ntop-2; m > l; m--){
			Hmods[m] = Hmods[m-1];
			Hindices[m] = Hindices[m-1];
		}
		Hmods[l]=s;
		Hindices[l] = GenFncIndex(i,j,k);
		if(ntop<tmax) ntop++;//we are adding rather than replacing a line
		return Hmods.back();
	}
}
/* compute the best focal distance at a 2d grid in L, Xi=Jz/L
 * on input gridL, which is a grid in Jcirc, times gridXi is a uniform grid on (0,1)
*/
void createGridFocalDistance(const potential::BasePotential& pot,
			     std::vector<double>& gridL, std::vector<double>& gridXi,
			     math::Matrix<double>& grid2dD, math::Matrix<double>& grid2dR)     
{
	int sizeL = gridL.size(), sizeXi = gridXi.size();
	math::Matrix<double> grid2dL(sizeL, sizeXi);
	for(int iL=1; iL<sizeL-1; iL++){//omit bdy values
		double Rc, Vc, Jz, Jc  = gridL[iL];
		double E = E_circ(pot, Jc, &Rc, &Vc);
		std::vector<double> L_vals(sizeXi);
		std::vector<double> Xi_vals(sizeXi);
		std::vector<double> D_vals(sizeXi);
		std::vector<double> R_vals(sizeXi);
		for(int iXi=1; iXi<sizeXi-1; iXi++){
			double Jphi = Jc * (1-gridXi[iXi]);
			double Rsh, FD;
			FD  = estimateFocalDistanceShellOrbit(pot, E, Jphi, &Rsh, &Jz);
			double L=fabs(Jphi)+Jz;
			L_vals[iXi]=L; Xi_vals[iXi]=Jz/L;
			D_vals[iXi]=FD; R_vals[iXi]=Rsh;
			if(iXi>0 && Xi_vals[iXi]<Xi_vals[iXi-1])
				printf("createGridFocalDistances: non-monotonic Xi:\n",
				       "%d %f %f\n",iXi,Xi_vals[iXi-1],Xi_vals[iXi]);
		}
		// bdy values
		L_vals[0]=gridL[iL]; L_vals[sizeXi-1]=L_vals[sizeXi-2];
		Xi_vals[0]=0; Xi_vals[sizeXi-1]=1;
		D_vals[0]=D_vals[1]; D_vals[sizeXi-1]=D_vals[sizeXi-2];
		R_vals[0]=R_vals[1]; R_vals[sizeXi-1]=R_vals[sizeXi-2];		
		math::LinearInterpolator interpL(Xi_vals,L_vals);
		math::LinearInterpolator interpD(Xi_vals,D_vals);
		math::LinearInterpolator interpR(Xi_vals,R_vals);
		for(int iXi=0; iXi<sizeXi; iXi++){//interpolae L, D onto regular grid in Xi
			interpL.evalDeriv(gridXi[iXi],&grid2dL(iL,iXi));
			interpD.evalDeriv(gridXi[iXi],&grid2dD(iL,iXi));
			interpR.evalDeriv(gridXi[iXi],&grid2dR(iL,iXi));
		}
	}
	//We now need to fill in rows iL=0, iL=sizeL-1
	for(int iXi=0; iXi<sizeXi; iXi++){
		grid2dL(0,iXi) = 0; grid2dL(sizeL-1,iXi) = 1.01*grid2dL(sizeL-2,iXi);
		grid2dD(0,iXi) = 0; grid2dD(sizeL-1,iXi) = grid2dD(sizeL-2,iXi);
		grid2dR(0,iXi) = 0; grid2dR(sizeL-1,iXi) = grid2dR(sizeL-2,iXi);
	}
	//Now grid2dD contains D on regular grid in Xi but irregular
	//values of L that are stored in grid2dL
	for(int iXi=0; iXi<sizeXi; iXi++){
		std::vector<double> L_vals(sizeL);
		std::vector<double> D_vals(sizeL);
		std::vector<double> R_vals(sizeL);
		for(int iL=0; iL<sizeL; iL++){
			L_vals[iL]=grid2dL(iL,iXi);
			D_vals[iL]=grid2dD(iL,iXi);
			R_vals[iL]=grid2dR(iL,iXi);
			if(iL>0 && L_vals[iL]<=L_vals[iL-1])
				printf("createGridFocalDistance: non-monotonic L_vals %g %g\n",
				       L_vals[iL-1],L_vals[iL]);
		}
		math::LinearInterpolator DL(L_vals,D_vals);
		math::LinearInterpolator RL(L_vals,R_vals);
		for(int iL=0; iL<sizeL; iL++){
			DL.evalDeriv(gridL[iL],&grid2dD(iL,iXi));
			RL.evalDeriv(gridL[iL],&grid2dR(iL,iXi));
		}
	}
}

double H_dHdX(const potential::BasePotential& pot, const coord::PosMomCyl Rzphi,
	      coord::PosMomCyl& dHdX){
	double Phi; coord::GradCyl grad;
	pot.eval(Rzphi, &Phi, &grad); 
	dHdX.R = grad.dR - pow_2(Rzphi.pphi/Rzphi.R)/Rzphi.R;
	dHdX.z = grad.dz; dHdX.phi = grad.dphi; 
	dHdX.pR = Rzphi.pR; dHdX.pz = Rzphi.pz; dHdX.pphi = Rzphi.pphi/pow_2(Rzphi.R);
	return .5*(pow_2(Rzphi.pR) + pow_2(Rzphi.pz) + pow_2(Rzphi.pphi/Rzphi.R)) + Phi;
}

coord::PosMomCyl PosMomDerivs(const coord::PosMomSph& rtheta, const coord::UVSph& cs,
			      const actions::DerivAct<coord::Sph>* dJ,//input from TM
			      actions::DerivAct<coord::Cyl>& dJC) {
	double snv,csv;
	math::sincos(rtheta.theta,snv,csv);
	coord::PosUVSph uv(asinh(rtheta.r/cs.Delta),rtheta.theta,rtheta.phi,cs);
	double chu=cosh(uv.u), shu=sinh(uv.u);
	double dudr=1/(cs.Delta*chu);
	coord::MomUVSph pp(rtheta.pr/dudr,rtheta.ptheta,rtheta.pphi);
	coord::PosMomUVSph uvp(uv,pp);
	coord::PosMomCyl Rzphi = coord::toPosMomCyl(uvp);
	if( !dJ ) return Rzphi;
	double dpudr=cs.Delta*shu*dudr*rtheta.pr;
	double dpvdtheta=0;
	double dpudpr=1/dudr;
	double dpvdptheta=1;
	DerivActUVSph dJUV;
	dJUV.dbyJr.u   = dudr*dJ->dbyJr.r;
	dJUV.dbyJz.u   = dudr*dJ->dbyJz.r;
	dJUV.dbyJphi.u = dudr*dJ->dbyJphi.r;
	dJUV.dbyJr.v   = dJ->dbyJr.theta;
	dJUV.dbyJz.v   = dJ->dbyJz.theta;
	dJUV.dbyJphi.v = dJ->dbyJphi.theta;
	dJUV.dbyJr.phi = dJ->dbyJr.phi;//all 3 non-zero
	dJUV.dbyJz.phi = dJ->dbyJz.phi;
	dJUV.dbyJphi.phi = dJ->dbyJphi.phi;

	dJUV.dbyJr.pu   = dpudr*dJ->dbyJr.r + dpudpr*dJ->dbyJr.pr;
	dJUV.dbyJz.pu   = dpudr*dJ->dbyJz.r + dpudpr*dJ->dbyJz.pr;
	dJUV.dbyJphi.pu = dpudr*dJ->dbyJphi.r+dpudpr*dJ->dbyJphi.pr;
	dJUV.dbyJr.pv   = dJ->dbyJr.ptheta;
	dJUV.dbyJz.pv   = dJ->dbyJz.ptheta;
	dJUV.dbyJphi.pv = dJ->dbyJphi.ptheta;
	dJUV.dbyJr.pphi = dJ->dbyJr.pphi;//both these vanish!
	dJUV.dbyJz.pphi = dJ->dbyJz.pphi;
	dJUV.dbyJphi.pphi=dJ->dbyJphi.pphi;//unity

	double dRdu = cs.Delta*chu*snv, dRdv = cs.Delta*shu*csv;
	double dzdu = cs.Delta*shu*csv, dzdv =-cs.Delta*chu*snv;
	double det = cs.Delta*(pow_2(shu)+pow_2(snv));
	double ddetdu=2*shu*chu*cs.Delta, ddetdv=2*snv*csv*cs.Delta;
	double dpRdu = ( shu*snv*uvp.pu+chu*csv*uvp.pv - Rzphi.pR*ddetdu)/det;
	double dpRdv = ( chu*csv*uvp.pu-shu*snv*uvp.pv - Rzphi.pR*ddetdv)/det;
	double dpzdu = ( chu*csv*uvp.pu-shu*snv*uvp.pv - Rzphi.pz*ddetdu)/det;
	double dpzdv = (-shu*snv*uvp.pu-chu*csv*uvp.pv - Rzphi.pz*ddetdv)/det;
	double dpRdpu = chu*snv/det, dpRdpv = shu*csv/det;
	double dpzdpu = shu*csv/det, dpzdpv = -chu*snv/det;

	dJC.dbyJr.R  = dRdu*dJUV.dbyJr.u  + dRdv*dJUV.dbyJr.v;
	dJC.dbyJz.R  = dRdu*dJUV.dbyJz.u  + dRdv*dJUV.dbyJz.v;
	dJC.dbyJphi.R= dRdu*dJUV.dbyJphi.u+ dRdv*dJUV.dbyJphi.v;
	dJC.dbyJr.z  = dzdu*dJUV.dbyJr.u  + dzdv*dJUV.dbyJr.v;
	dJC.dbyJz.z  = dzdu*dJUV.dbyJz.u  + dzdv*dJUV.dbyJz.v;
	dJC.dbyJphi.z= dzdu*dJUV.dbyJphi.u+ dzdv*dJUV.dbyJphi.v;
	dJC.dbyJr.phi  = dJUV.dbyJr.phi;  //dphidu*dJUV.dbyJr.u  + dphidv*dJUV.dbyJr.v;
	dJC.dbyJz.phi  = dJUV.dbyJz.phi;  //dphidu*dJUV.dbyJz.u  + dphidv*dJUV.dbyJz.v;
	dJC.dbyJphi.phi= dJUV.dbyJphi.phi;//dphidu*dJUV.dbyJphi.u+ dphidv*dJUV.dbyJphi.v;

	dJC.dbyJr.pR = dpRdu*dJUV.dbyJr.u + dpRdv*dJUV.dbyJr.v + dpRdpu*dJUV.dbyJr.pu + dpRdpv*dJUV.dbyJr.pv;
	dJC.dbyJz.pR = dpRdu*dJUV.dbyJz.u + dpRdv*dJUV.dbyJz.v + dpRdpu*dJUV.dbyJz.pu + dpRdpv*dJUV.dbyJz.pv;
	dJC.dbyJphi.pR = dpRdu*dJUV.dbyJphi.u + dpRdv*dJUV.dbyJphi.v + dpRdpu*dJUV.dbyJphi.pu + dpRdpv*dJUV.dbyJphi.pv;
	dJC.dbyJr.pz = dpzdu*dJUV.dbyJr.u + dpzdv*dJUV.dbyJr.v + dpzdpu*dJUV.dbyJr.pu + dpzdpv*dJUV.dbyJr.pv;
	dJC.dbyJz.pz = dpzdu*dJUV.dbyJz.u + dpzdv*dJUV.dbyJz.v + dpzdpu*dJUV.dbyJz.pu + dpzdpv*dJUV.dbyJz.pv;
	dJC.dbyJphi.pz = dpzdu*dJUV.dbyJphi.u + dpzdv*dJUV.dbyJphi.v + dpzdpu*dJUV.dbyJphi.pu + dpzdpv*dJUV.dbyJphi.pv;
	dJC.dbyJr.pphi = 0;
	dJC.dbyJz.pphi = 0;
	dJC.dbyJphi.pphi = 1;
	if(std::isnan(dJC.dbyJr.pR)){
		printf("PosMomD\n%f %f %f %f %f %f %f %f\n",dpRdu, dJUV.dbyJr.u, dpRdv, dJUV.dbyJr.v,
		       dpRdpu, dJUV.dbyJr.pu, dpRdpv, dJUV.dbyJr.pv);
		exit(0);
	}
	return Rzphi;
}

class Iso{
	private:
		double L,Jr;
	public:
		Iso(double _L,double _Jr) : L(_L), Jr(_Jr) {}
		double b2cE(double Js) const{//general E
			return .5*pow_2(Js*Js)/pow_2(Jr+.5*(L+sqrt(L*L+4*Js*Js)));
		}
		double b2cEc(double Js) const{//circular E
			return .5*pow_2(Js*Js)/pow_2(.5*(L+sqrt(L*L+4*Js*Js)));
		}
		double g(double Js) const{//r=b*g(Js)
			double g2=(pow((L+sqrt(L*L+4*Js*Js))/(2*Js),4)-1);
			if(g2<0) printf("g2<0: %f\n",g2);
			return sqrt(g2);
		}
		double cob(double Js) const{//c/b
			return .5*pow_2(Js)/b2cE(Js)-1;;
		}
		double ecc(double Js) const{
			double boc=1/cob(Js);
			double e2=1-pow_2(L/Js)*boc*(1+boc);
			return e2<1? sqrt(e2) : 0;
		}
		double f(double b,double Js,double& e) const{//ratio of forces aopo/peri
			double cb=cob(Js), boc=1/cb;
			e=ecc(Js);
			double up=1+e, um=1-e;
			double rp=b*cb*sqrt(up*(up+2*boc)), ap=sqrt(b*b+rp*rp);
			double rm=b*cb*sqrt(um*(um+2*boc)), am=sqrt(b*b+rm*rm);
			return pow_2((b+am)/(b+ap))*am/ap*rp/rm;
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
		    pot(_pot), Is(_Is), Rsh(_Rsh) {}
		virtual void evalDeriv(double Js,double* val, double *deriv=0, double* deriv2=0) const{
			double b = Rsh/Is.g(Js);
			double e, f_apo_peri = Is.f(b,Js,e);
			double c = Is.cob(Js)*b;
			double F[2];
			for(int k=-1; k<2; k+=2){
				double u=1+k*e;
				double Phi, r=c*sqrt(u*(u+2*b/c));
				coord::PosCyl Rz(r,0,0); coord::GradCyl grad;
				pot.eval(Rz,&Phi,&grad);
				F[(k+1)/2]=grad.dR;
			}
			*val = f_apo_peri -  F[1]/F[0];
		}
		virtual unsigned int numDerivs(void) const{
			return 0;
		}

};

/// create the array of indices of the generating function with all terms up to the given maximum order
static GenFncIndices makeGridIndices(int irmax, int izmax)
{   /// NOTE: here we specialize for the case of axisymmetric systems!
	GenFncIndices indices;
	for(int ir=0; ir<=irmax; ir++)
		for(int iz=-izmax; iz<=(ir==0?-2:izmax); iz+=2)
			indices.push_back(GenFncIndex(ir, iz, 0));
	return indices;
}

/// return the absolute value of an element in a map, or zero if it doesn't exist
static inline double absvalue(const std::map< std::pair<int,int>, double >& indPairs, int ir, int iz)
{
	if(indPairs.find(std::make_pair(ir, iz)) != indPairs.end())
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
		       dXdJ.dbyJr.pR* dHby.pR+ dXdJ.dbyJr.pz* dHby.pz+ dXdJ.dbyJr.pphi* dHby.pphi,
		       dXdJ.dbyJz.R * dHby.R + dXdJ.dbyJz.z * dHby.z + dXdJ.dbyJz.phi * dHby.phi +
		       dXdJ.dbyJz.pR* dHby.pR+ dXdJ.dbyJz.pz* dHby.pz+ dXdJ.dbyJz.pphi* dHby.pphi,
		       dXdJ.dbyJphi.R * dHby.R + dXdJ.dbyJphi.z * dHby.z + dXdJ.dbyJphi.phi * dHby.phi +
		       dXdJ.dbyJphi.pR* dHby.pR+ dXdJ.dbyJphi.pz* dHby.pz+ dXdJ.dbyJphi.pphi* dHby.pphi);
}
//class to help with SoS computation
class CrossingFinder : public math::IFunction{
	private:
		const BaseTorus* T;
		const double thetaT_r;
	public:
		CrossingFinder(const BaseTorus* _T,const double& _thetaT) : T(_T), thetaT_r(_thetaT) {}
		virtual void evalDeriv(const double thetaT_z, double* value, double* deriv=NULL, double* deriv2=NULL) const{
			Angles thetaT(thetaT_r,thetaT_z,0);
			coord::PosMomCyl Rz(T->from_toy(thetaT));
			*value = Rz.z;
		}
		virtual unsigned int numDerivs(void) const{
			return 0;
		}
};

}//internal

PerturbingHamiltonian& PerturbingHamiltonian::operator *= (const double a){
	for(int i=0; i<values.size(); i++)
		values[i]*=a;
	return *this;
}
PerturbingHamiltonian PerturbingHamiltonian::operator * (const double a){
	PerturbingHamiltonian H2(indices, values);
	H2 *=a;
	return H2;
}
PerturbingHamiltonian& PerturbingHamiltonian::operator += (const PerturbingHamiltonian& H){
	PerturbingHamiltonian H2(H);
	for(int i=0; i<indices.size(); i++){
		int mr=indices[i].mr, mz=indices[i].mz, mphi=indices[i].mphi;
		std::vector<std::complex<double> >::iterator jt = H2.values.begin();
		for(GenFncIndices::iterator it =  H2.indices.begin(); it != H2.indices.end();){
			if(mr==(*it).mr && mz==(*it).mz && mphi==(*it).mphi){
				values[i]+=(*jt);
				it=H2.indices.erase(it);
				jt=H2.values.erase(jt);
				break;
			} else {
				it++; jt++;
			}
		}
	}
	if(H2.indices.size() > 0){
		for(int j=0; j<H2.indices.size(); j++){
			indices.push_back(H2.indices[j]);
			values.push_back(H2.values[j]);
		}
	}
	return *this;
}

PerturbingHamiltonian PerturbingHamiltonian::operator + (const  PerturbingHamiltonian& H){
	PerturbingHamiltonian H2 = *this;
	H2 += H;
	return H2;
}

#ifdef TESTIT
#include "actions_test_torus.cpp"
#endif

coord::PosMomCyl BaseTorus::from_toy(const Angles& thetaT) const{
	ActionAngles aaT(GF.toyJ(J, thetaT), thetaT);
	return from_aaT(aaT);
}
		
	
	
	
coord::PosMomCyl BaseTorus::from_aaT(const ActionAngles& aaT) const{//input toy J & theta
	coord::PosMomSph rtheta(TM.aa2pq(aaT));//image from isochrone
	coord::PosUVSph uv(asinh(rtheta.r/cs.Delta), rtheta.theta, rtheta.phi, cs);
	double chu=cosh(uv.u), shu=sinh(uv.u);
	double dudr=1/(cs.Delta*chu);
	coord::MomUVSph pp(rtheta.pr/dudr, rtheta.ptheta, rtheta.pphi);
	return coord::toPosMomCyl(coord::PosMomUVSph(uv, pp));
}	
coord::PosMomCyl BaseTorus::from_true(const Angles& theta) const{//input true angles
	ActionAngles aaT(GF.true2toy(ActionAngles(J,theta)));//toy AAs computed from true
	return from_aaT(aaT);
}

/* Position from toy angle plus dR/dthetaT at fixed J (which causes JT
 * to vary with thetaT)
*/
coord::PosCyl BaseTorus::PosDerivs(const Angles& thetaT,
				   actions::DerivAng<coord::Cyl>& dRA, double* det) const {
	ActionAngles aaT(GF.toyJ(J,thetaT), thetaT);
	actions::DerivAct<coord::Sph> dJ;
	actions::DerivAng<coord::Sph> dA;
	const coord::PosMomSph rtheta(TM.aa2pq(aaT, NULL, &dJ, &dA));
	double snv,csv;
	math::sincos(rtheta.theta,snv,csv);
	coord::PosUVSph uv(asinh(rtheta.r/cs.Delta), rtheta.theta, rtheta.phi, cs);
	double chu=cosh(uv.u), shu=sinh(uv.u);
	double dudr=1/(cs.Delta*chu);
	double dRdu = cs.Delta*chu*snv, dRdv = cs.Delta*shu*csv;
	double dzdu = cs.Delta*shu*csv, dzdv =-cs.Delta*chu*snv;

	dRA.dbythetar.R  = dRdu*dudr*dA.dbythetar.r  + dRdv*dA.dbythetar.theta;
	dRA.dbythetaz.R  = dRdu*dudr*dA.dbythetaz.r  + dRdv*dA.dbythetaz.theta;
	dRA.dbythetaphi.R= dRdu*dudr*dA.dbythetaphi.r+ dRdv*dA.dbythetaphi.theta;
	dRA.dbythetar.z  = dzdu*dudr*dA.dbythetar.r  + dzdv*dA.dbythetar.theta;
	dRA.dbythetaz.z  = dzdu*dudr*dA.dbythetaz.r  + dzdv*dA.dbythetaz.theta;
	dRA.dbythetaphi.z= dzdu*dudr*dA.dbythetaphi.r+ dzdv*dA.dbythetaphi.theta;
	dRA.dbythetar.phi  = dA.dbythetar.phi;
	dRA.dbythetaz.phi  = dA.dbythetaz.phi;
	dRA.dbythetaphi.phi= 1;

	actions::DerivAct<coord::Cyl> dRJ;
	dRJ.dbyJr.R  = dRdu*dudr*dJ.dbyJr.r  + dRdv*dJ.dbyJr.theta;
	dRJ.dbyJz.R  = dRdu*dudr*dJ.dbyJz.r  + dRdv*dJ.dbyJz.theta;
	dRJ.dbyJphi.R= dRdu*dudr*dJ.dbyJphi.r+ dRdv*dJ.dbyJphi.theta;
	dRJ.dbyJr.z  = dzdu*dudr*dJ.dbyJr.r  + dzdv*dJ.dbyJr.theta;
	dRJ.dbyJz.z  = dzdu*dudr*dJ.dbyJz.r  + dzdv*dJ.dbyJz.theta;
	dRJ.dbyJphi.z= dzdu*dudr*dJ.dbyJphi.r+ dzdv*dJ.dbyJphi.theta;
	dRJ.dbyJr.phi  = dJ.dbyJr.phi;
	dRJ.dbyJz.phi  = dJ.dbyJz.phi;
	dRJ.dbyJphi.phi= 1;
	
	actions::DerivAng<coord::Cyl> dJA = GF.dJdt(thetaT); //dJT/dthetaT
	dRA.dbythetar.R   += dRJ.dbyJr.R   * dJA.dbythetar.R + dRJ.dbyJz.R  * dJA.dbythetaz.z;
	dRA.dbythetaz.R   += dRJ.dbyJr.R   * dJA.dbythetaz.R + dRJ.dbyJz.R  * dJA.dbythetaz.z;
	dRA.dbythetar.z   += dRJ.dbyJr.z   * dJA.dbythetar.R + dRJ.dbyJz.z  * dJA.dbythetar.z;
	dRA.dbythetaz.z   += dRJ.dbyJr.z   * dJA.dbythetaz.R + dRJ.dbyJz.z  * dJA.dbythetaz.z;
	dRA.dbythetar.phi += dRJ.dbyJr.phi * dJA.dbythetar.R + dRJ.dbyJz.phi* dJA.dbythetar.z;
	dRA.dbythetaz.phi += dRJ.dbyJr.phi * dJA.dbythetaz.R + dRJ.dbyJz.phi* dJA.dbythetaz.z;

	if(det){//assume only index.mphi=0 non-zero
		(*det) = dRA.dbythetar.R*dRA.dbythetaz.z-dRA.dbythetar.z*dRA.dbythetaz.R;
	}
	return coord::toPosCyl(uv);
}
void BaseTorus::zSoS(std::vector<double>& R,std::vector<double>& vR,int N,
		 double& Rmin, double& Rmax, double& Vmax) const{
	const double tol=1e-3;
	Rmin=1e10, Rmax=0, Vmax=0;
	for(int i=0; i<N; i++){
		double thetaT_r = 2*M_PI/(double)N*(-N/2+i);
		CrossingFinder CF(this,thetaT_r);
		double thetaT_z, dtheta=.1, th_min=-.5*M_PI,
		th_max = th_min+dtheta, z_min, z_max;
		CF.evalDeriv(th_min, &z_min); CF.evalDeriv(th_max,&z_max);
		while(z_min*z_max>0){// plod round looking for sign change
			th_min = th_max; z_min=z_max; th_max += dtheta;
			CF.evalDeriv(th_max,&z_max);
		}
		thetaT_z = math::findRoot(CF,th_min, th_max,tol);
		coord::PosMomCyl Rz(from_toy(Angles(thetaT_r,thetaT_z,0)));
		if(Rz.pz<-.0001){// keep going round
			do {
				th_min=th_max; z_min=z_max;
				th_max += dtheta; CF.evalDeriv(th_max,&z_max);
			} while(z_min*z_max>0);
			thetaT_z = math::findRoot(CF,th_min, th_max,tol);
			Rz = coord::PosMomCyl(from_toy(Angles(thetaT_r,thetaT_z,0)));
		}
		R.push_back(Rz.R); vR.push_back(Rz.pR);
		Rmin=fmin(Rmin,Rz.R); Rmax=fmax(Rmax,Rz.R); Vmax=fmax(Vmax,fabs(Rz.pR));
	}
	R.push_back(R[0]); vR.push_back(vR[0]);
}
std::vector<std::pair<coord::PosVelCyl,double> > BaseTorus::orbit(const Angles& theta0, double dt, double T) const{
	std::vector<std::pair<coord::PosVelCyl,double> > traj;
	double t=0;
	while(t<T){
		actions::Angles theta(theta0.thetar+freqs.Omegar*t,theta0.thetaz+freqs.Omegaz*t,
				      theta0.thetaphi+freqs.Omegaphi*t);
		traj.push_back(std::pair<coord::PosVelCyl,double>(coord::toPosVelCyl(from_true(theta)),t));
		t+=dt;
	}
	return traj;
}

// Returns true if (R,z,phi) is ever hit by the orbit, and false otherwise. If the 
// torus passes through the point given, this happens four times, in each case
// with a different velocity, but only two of these are independent:
// theta_r -> -theta_r with theta_z -> Pi-theta_z leaves (R,z) and J^T fixed
// but changes the sign of both velocities.
// | d(x,y,z)/d(Tr,Tl,phi) | is returned. The latter vanishes on the edge of the
// orbit, such that its inverse, the density of the orbit, diverges there
// (that's the reason why the density itself is not returned).

bool BaseTorus::containsPoint(const coord::PosCyl& p, std::vector<Angles>& As,
			      std::vector<coord::VelCyl>& Vs,
			      std::vector<double>& Jacobs, const double& tol) const{
	coord::PosMomCyl peri(from_true(Angles(0,.5*M_PI,0))), apo(from_true(Angles(M_PI,0,0))),
	top(from_true(Angles(M_PI,.5*M_PI,0)));
	double Rmin=.95*peri.R, Rmax=1.05*apo.R, zmax=1.05*fabs(top.z);
	//printf("containsPoint: Rmin %f, Rmax %f, zmax %f\n",Rmin,Rmax,zmax);
	if(p.R<Rmin || p.R>Rmax || fabs(p.z)>zmax) return false;
	locFinder LF(*this, p);
	double tolerance=1e-8;
	double params[2]={1,1}, result[2], dist, det;
	int maxNumIter=200;
	coord::PosMomCyl P1; Angles A1, Atrue;
	DerivAng<coord::Cyl> dA;
	int done, kmax=30, nfail=0;
	while(As.size()<4 && nfail<kmax){
		double kount=0;
		for(int k=0; k<kmax; k++) {
			done=math::nonlinearMultiFit(LF, params, tolerance, maxNumIter, result);
			A1=Angles(math::wrapAngle(result[0]), math::wrapAngle(result[1]), 0);
			P1=from_toy(A1);
			dist=sqrt(pow_2(p.R-P1.R)+pow_2(p.z-P1.z));
			if(dist<2*tol) {
				Atrue=GF.trueA(A1);
				if(is_new(Atrue,As)) break;
			} //else	printf("Max in containsPoint: %d %g\n", done, dist);
			params[0]+=.3; params[0]=math::wrapAngle(params[0]);
			params[1]+=.7;   params[1]=math::wrapAngle(params[1]);
			kount++;
			if(kount==kmax && As.size()==0) return false;
		}
		if(kount<kmax){
			As.push_back(Atrue); PosDerivs(A1, dA, &det);
			Vs.push_back(coord::VelCyl(P1.pR, P1.pz, P1.pphi/P1.R));
			Jacobs.push_back(fabs(det/GF.dtbydtT_Jacobian(A1)));
			A1.thetar = -A1.thetar; A1.thetaz = M_PI-A1.thetaz;
			P1=from_toy(A1); Atrue=GF.trueA(A1);			
			As.push_back(Atrue); PosDerivs(A1, dA, &det);
			Vs.push_back(coord::VelCyl(P1.pR, P1.pz, P1.pphi/P1.R));
			Jacobs.push_back(fabs(det/GF.dtbydtT_Jacobian(A1)));
		} else nfail++;
	}
	if(nfail>=kmax) printf("containsPoint error at Rz (%f %f) - %d angles \n",
			       p.R,p.z,As.size());
	return true;
}
double BaseTorus::density(const coord::PosCyl& Rz) const{
	std::vector<Angles> As; std::vector<coord::VelCyl> Vs;
	std::vector<double> Jacobs;
	const double tol=1e-6;
	if(!containsPoint(Rz,As,Vs,Jacobs,tol)) return 0;
	double rho=0;
	for(int i=0; i<As.size(); i++)
		rho+=1/Jacobs[i];
	return rho;
}
void BaseTorus::write(FILE* ofile) const{
	fprintf(ofile, "%g %g %g %g %g %g %g %g %g %g\n",
		  J.Jr, J.Jz, J.Jphi, freqs.Omegar, freqs.Omegaz, freqs.Omegaphi,
		  E, cs.Delta, TM.Js, TM.b);
	GF.write(ofile);
}
void BaseTorus::read(FILE* ifile) {
	double Delta, Js, b;
	fscanf_s(ifile, "%g %g %g %g %g %g %g %g %g %g\n",
		 &J.Jr, &J.Jz, &J.Jphi, &freqs.Omegar, &freqs.Omegaz, &freqs.Omegaphi,
		 &E, &Delta, &Js, &b);
	cs=coord::UVSph(Delta);
	TM=Isochrone(Js,b);
	GF.read(ifile);
}
Torus& Torus::operator *= (const double a){
	J*=a; freqs*=a; GF*=a; TM*=a; cs=coord::UVSph(cs.Delta*a); E*a;
	return *this;
}
Torus& Torus::operator += (const Torus& T){
	J+=T.J; freqs+=T.freqs; GF+=T.GF; TM+=T.TM;
	cs=coord::UVSph(cs.Delta+=T.cs.Delta); E+T.E;
	return *this;
}
const Torus Torus::operator * (const double a){
	Torus T2(J*a, freqs*a, GF*a, TM*a, coord::UVSph(cs.Delta*a), E*a);
	return T2;
}
const Torus Torus::operator + (const Torus& T){
	Torus T2(J+T.J, freqs+T.freqs, GF+T.GF, TM+T.TM, coord::UVSph(cs.Delta+T.cs.Delta), E+T.E);
	return T2;
}
const eTorus eTorus::operator + (const eTorus& T){
	eTorus T2(J+T.J, freqs+T.freqs, GF+T.GF, TM+T.TM, coord::UVSph(cs.Delta+T.cs.Delta),
		  E+T.E, pH+T.pH);
	return T2;
}

TorusGenerator::TorusGenerator(const potential::BasePotential& _pot, const double _tol) :
    pot(_pot), tol(_tol), invPhi0(1./_pot.value(coord::PosCyl(0,0,0))), tmax(10) {
	std::vector<double> gridR = potential::createInterpolationGrid(pot, ACCURACY_INTERP2);
	int sizeL = gridR.size(), sizeXi=20;
	//printf("L grid: %d\n",sizeL);
	math::ScalingSemiInf Lscale;
	std::vector<double> gridL(sizeL);
	std::vector<double> gridLscaled(sizeL);
	std::vector<double> gridXi(sizeXi);
	for(int i=0; i<sizeL; i++){
		gridL[i] = gridR[i] * potential::v_circ(pot, gridR[i]);
//		gridLscaled[i]=i/(double)(sizeL-1);
//		gridL[i]=math::unscale(Lscale,gridLscaled[i]);
	}
	for(int i=0; i<sizeXi; i++)
		gridXi[i]=i/(double)(sizeXi-1);//EV wld call this Xiscaled
	math::Matrix<double> grid2dD(sizeL,sizeXi);
	math::Matrix<double> grid2dR(sizeL,sizeXi);
	createGridFocalDistance(pot, gridL, gridXi, grid2dD, grid2dR);
	interpD = math::LinearInterpolator2d(gridL, gridXi, grid2dD);
	interpR = math::LinearInterpolator2d(gridL, gridXi, grid2dR);
}

double TorusGenerator::computeHamiltonianAtPoint(const double params[],
	const unsigned int indPoint, Actions *dHdJ, double *derivGenFnc) const
{
    // Generating function computes the toy actions from the real actions
    // at the given point in the grid of toy angles grid
	ActionAngles toyAA = GFF->toyActionAngles(indPoint, params);

    // do not allow to stray into forbidden region of negative actions
	if(toyAA.Jr<0 || toyAA.Jz<0){
		printf("J<0: %f %f ",toyAA.Jr,toyAA.Jz);
		return NAN;
	}
    // aa2pq computes the position and velocity from the toy actions and angles,
    // and optionally their derivatives w.r.t. toy actions and toy map parameters,
	DerivAct<coord::Sph> derivAct;//derivs of coords wrt Js
	coord::PosMomSph rtheta = TM.aa2pq(toyAA, NULL, derivGenFnc!=NULL ? &derivAct : NULL);
	actions::DerivAct<coord::Cyl> dXdJ;
	coord::PosMomCyl Rzphi = PosMomDerivs(rtheta, cs,
					      derivGenFnc!=NULL ? &derivAct : NULL, dXdJ);
    // obtain the value of the real Hamiltonian at the given point and its
    // derivatives w.r.t. coordinates/momenta
	coord::PosMomCyl dHdX;
	double H = H_dHdX(pot, Rzphi, dHdX);

    // derivatives of Hamiltonian w.r.t. parameters of gen.fnc.
	if(derivGenFnc) {
		Actions dHby = dHbydJ(dHdX, dXdJ);// derivative of Hamiltonian by toy actions
		if(indPoint==-30 || indPoint==-10){
			printf("(H %f dHby %f %f\n",H,dHby.Jr,dHby.Jz);
			printf("dXdJ %f %f %f %f %f %f %f %f\n",dXdJ.dbyJr.R,dXdJ.dbyJr.z,dXdJ.dbyJz.R,dXdJ.dbyJz.z,dXdJ.dbyJr.pR,dXdJ.dbyJr.pz,dXdJ.dbyJz.pR,dXdJ.dbyJz.pz);
			printf("[dHdX %f %f %f %f %f %f\n] ",dHdX.R,dHdX.z,dHdX.phi,dHdX.pR,dHdX.pz,dHdX.pphi);
		}
		if(dHdJ) *dHdJ = dHby;
		for(unsigned int p = 0; p<GFF->numParams(); p++) {
	    // derivs of toy actions by gen.fnc.params
			Actions dbyS = GFF->deriv(indPoint, p);
	    // derivs of Hamiltonian by gen.fnc.params
			double  dHdS = dHby.Jr * dbyS.Jr + dHby.Jz * dbyS.Jz + dHby.Jphi * dbyS.Jphi;
			derivGenFnc[p] = dHdS;
			if(std::isnan(dHdS)){
				printf("AtPoint  %f %f %f %f %f %f\n",dHby.Jr,dHby.Jz,dHby.Jphi,dbyS.Jr,dbyS.Jz,dbyS.Jphi);
				printf("%f %f %f %f %f %f\n",dXdJ.dbyJr.R,dXdJ.dbyJr.z,dXdJ.dbyJr.phi,
				      dXdJ.dbyJr.pR,dXdJ.dbyJr.pz,dXdJ.dbyJr.pphi);
				printf("%f %f %f %f %f %f\n",dXdJ.dbyJz.R,dXdJ.dbyJz.z,dXdJ.dbyJz.phi,
				       dXdJ.dbyJz.pR,dXdJ.dbyJz.pz,dXdJ.dbyJz.pphi);
				printf("%f %f %f %f %f %f\n",dXdJ.dbyJphi.R,dXdJ.dbyJphi.z,dXdJ.dbyJphi.phi,
				       dXdJ.dbyJphi.pR,dXdJ.dbyJphi.pz,dXdJ.dbyJphi.pphi);
				exit(0);
			}
		}
	}
//	if(derivGenFnc)printf("(%f %f)",H,derivGenFnc[10]);
	return H;
}
double TorusGenerator::Hamilton(const BaseTorus& T, const potential::BasePotential* ePot, const Angles& theta)
{
	coord::PosMomCyl Rz(T.from_true(theta));
	double H=.5*(pow_2(Rz.pR) + pow_2(Rz.pz) + pow_2(Rz.pphi/Rz.R)) + pot.value(Rz);
	coord::PosCyl pos(Rz.R,Rz.z,Rz.phi);
	return ePot? H+ePot->value(pos) : H;
}
PerturbingHamiltonian TorusGenerator::get_pH (const BaseTorus& T, int nf, bool ifp,
	const potential::BasePotential* ePot){//Fourier analyses H
	int nfr=nf, nfz=nf, nfp=nf/4;
	double* h = new double[nfr*nfz*nfp];
	Angles thetas;
	double dtr=2*M_PI/(double)nfr, dtz=2*M_PI/(double)nfz, dtp=2*M_PI/(double)nfp;
	for(int i=0; i<nfr; i++){
		thetas.thetar = i*dtr;
		for(int j=0; j<nfz/2; j++){
			thetas.thetaz = j*dtz;
			for(int k=0; k<nfp/2; k++){
				thetas.thetaphi = k*dtp;
				double a=Hamilton(T, ePot, thetas);
				//printf("%g ",a);
				h[nfz*nfp*i+nfp*j+k]=a;
				h[nfz*nfp*i+nfp*(nfz/2+j)+k]=a;//N-S symmetry
				h[nfz*nfp*i+nfp*j+nfp/2+k]=a; //bi-symmetry
				h[nfz*nfp*i+nfp*(nfz/2+j)+nfp/2+k]=a;//bi-symmetry
			}
		}
	}
	rlft3(h,nfr,nfz,nfp,1);
	std::vector<double> Hmods;
	GenFncIndices Hindices;
	int ntop=0;
	double hmax=0;
	for(int i=0; i<nfr; i++){//find largest perturbing terms
		for(int j=0; j<nfz; j++){
			for(int k=0; k<nfp/2; k++){
				double s=sqrt(pow_2(h[nfz*nfp*i+nfp*j+2*k])
					      +pow_2(h[nfz*nfp*i+nfp*j+2*k+1]));
				if(ntop<tmax || s>hmax){
					hmax=insertLine(ntop, tmax, s, i, j, k, Hmods, Hindices);
					//printf("%zd(%d %d %d)",Hindices.size(),Hindices.back().mr,Hindices.back().mz,Hindices.back().mphi);
				}
			}
		}
	}
	std::vector<std::complex<double> > Hvalues;
	double n=(nfr*nfz*nfp);
	for(int i=0; i<ntop; i++){
		Hvalues.push_back(std::complex<double>(h[nfz*nfp*Hindices[i].mr+nfp*Hindices[i].mz+2*Hindices[i].mphi]/n,
				 h[nfz*nfp*Hindices[i].mr+nfp*Hindices[i].mz+2*Hindices[i].mphi+1]/n));
	}
	if(ifp) printf("Terms in perturbing H (*100)\n");
	for(int i=0; i<ntop; i++){
		if(Hindices[i].mr>=nfr/2) Hindices[i].mr -= nfr;
		if(Hindices[i].mz>=nfz/2) Hindices[i].mz -= nfz;
		if(Hindices[i].mphi>=nfp/2) Hindices[i].mphi -= nfp;
		if(ifp) printf("(%3d %3d %3d) (%g %g)\n",
			       Hindices[i].mr, Hindices[i].mz, Hindices[i].mphi,
			       100*Hvalues[i].real(), 100*Hvalues[i].imag());
	}
	delete[] h;
	return PerturbingHamiltonian(Hindices,Hvalues);
}

void TorusGenerator::evalDeriv(const double params[],
			       double* deltaHvalues, double* dHdParams) const
{
	const unsigned int numPoints = GFF->numPoints();
	const unsigned int numParams = GFF->numParams();

    // we need to store the values of Hamiltonian at grid points even if this is not requested,
    // because they are used to correct the entries of the Jacobian matrix
    // to account for the fact that the mean <H> also depends on the parameters
	std::vector<double> Hvalues(numPoints);
	double Havg = 0;  // accumulator for the average Hamiltonian
    // loop over grid of toy angles
	for(unsigned int indPoint=0; indPoint < numPoints; indPoint++) {
		double H = computeHamiltonianAtPoint(params, indPoint, NULL,
			dHdParams ? dHdParams + indPoint * numParams : NULL);
		if(std::isnan(H)){
			Havg += 1;
			Hvalues[indPoint] = 1;
		} else {
	// accumulate the average value and store the output
			Havg += H;
			Hvalues[indPoint] = H;
		}
	}

    // convert from  H_k  to  deltaH_k = H_k - <H>
	Havg /= numPoints;
	if(deltaHvalues) {
		double disp = 0;
		for(unsigned int indPoint=0; indPoint < numPoints; indPoint++) {
			deltaHvalues[indPoint] = Hvalues[indPoint] - Havg;
			disp += pow_2(deltaHvalues[indPoint]);
		}
	}
    // convert derivatives:  d(deltaH_k) / dP_p = dH_k / dP_p - d<H> / dP_p
	if(dHdParams) {
		std::vector<double> dHavgdP(numPoints);
		for(unsigned int p=0; p<numParams; p++) dHavgdP[p]=0;
		for(unsigned int pp=0; pp < numPoints * numParams; pp++){
			dHavgdP[pp % numParams] += dHdParams[pp] / numPoints;
			if(std::isnan(dHavgdP[pp % numParams])){
				printf("nan@ %d %d %g\n",pp%numParams,pp/numParams,dHdParams[pp]);
				exit(0);
			}
		}
		for(unsigned int pp=0; pp < numPoints * numParams; pp++) {
			unsigned int indPoint = pp / numParams;
			unsigned int indParam = pp % numParams;
			dHdParams[pp] = dHdParams[pp] - dHavgdP[indParam];
		}
	}
}
double TorusGenerator::computeHamiltonianDisp(const std::vector<double> &params, double& Hbar)
{
	const unsigned int numParams = GFF->numParams();
	math::Averager Havg;
	int Nnan=0;
	for(unsigned int indPoint=0; indPoint < GFF->numPoints(); indPoint++){
		double H = computeHamiltonianAtPoint(&params[0], indPoint);
		Havg.add(H); if(std::isnan(H)) Nnan++;
	}
	Hbar = Havg.mean();
	NANfrac = (double)Nnan/(double)GFF->numPoints();
	return Havg.disp();    
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
double TorusGenerator::fitAngleMap(double params[], double& Hbar, Frequencies& freqs, GenFncDerivs& derivs) const
{
	unsigned int numPoints = GFF->numPoints();
	unsigned int numParams = GFF->numParams();
    // the matrix of coefficients shared between three linear systems
	math::Matrix<double> coefsdHdS(numPoints, numParams+1);
    // tmp storage for dH/dS
	std::vector<double> derivGenFnc(numParams+1);
    // derivs of Hamiltonian by toy actions (RHS vectors)
	std::vector<double> dHdJr(numPoints), dHdJz(numPoints), dHdJphi(numPoints);
    // accumulator for computing dispersion in H
	math::Averager Havg;

    // loop over grid of toy angles
	for(unsigned int indPoint=0; indPoint < numPoints; indPoint++) {
		Actions dHby;  // derivative of Hamiltonian by toy actions
		double H=computeHamiltonianAtPoint(&params[0], indPoint,
			&dHby, &derivGenFnc.front());
		if(!std::isnan(H)) Havg.add(H);
		else printf("nan @ %d ",indPoint);
	// fill the elements of each of three rhs vectors
		dHdJr  [indPoint] = dHby.Jr;
		dHdJz  [indPoint] = dHby.Jz;
		dHdJphi[indPoint] = dHby.Jphi;
	// fill the matrix row
		coefsdHdS(indPoint, 0) = 1;  // matrix coef for omega
		for(unsigned int p=0; p<numParams; p++)
			coefsdHdS(indPoint, p+1) = -derivGenFnc[p];  // matrix coef for dS_p/dJ
	}
	Hbar = Havg.mean();
    // solve the overdetermined linear system in the least-square sense:
    // step 1: prepare the SVD of coefs matrix
	math::SVDecomp SVD(coefsdHdS);

    // step 2: solve three linear systems with the same matrix but different rhs
	std::vector<double> dSdJr(SVD.solve(dHdJr)), dSdJz(SVD.solve(dHdJz)),
					dSdJphi(SVD.solve(dHdJphi));

    // store output
	freqs.Omegar   = dSdJr[0];
	freqs.Omegaz   = dSdJz[0];
	freqs.Omegaphi = dSdJphi[0];
	derivs.resize(numParams);
	for(unsigned int p=0; p<numParams; p++) {
		derivs[p].Jr   = dSdJr[p+1];
		derivs[p].Jz   = dSdJz[p+1];
		derivs[p].Jphi = dSdJphi[p+1];
	}
	return Havg.disp();
}
double TorusGenerator::getRsh(Actions& J){
	const double L = J.Jz + fabs(J.Jphi);
	return interpR.value(L, J.Jz/L);
}
void TorusGenerator::setConsts(actions::Actions J, double& Jscale){
	const double relToler=1e-4;
	const double L = fabs(J.Jphi) + J.Jz, Xi = J.Jz/L;
	const double Jtot = L + J.Jr;
	Jscale = J.Jr + J.Jz;
	Iso Is(L,J.Jr);
	Delta = interpD.value(L,Xi);
	cs = coord::UVSph(Delta);
	Rsh   = interpR.value(L,Xi);
	freqScale = potential::v_circ(pot, Rsh)/Rsh; //frequency scale set
	//For any Js b=Rsh/Is.g(Js) f_apo_peri=Is.f(b.Js,e) & where
	//this matches F_apo_peri in pot we pick Js and b
	double Jsmax=1.1*Jtot, Jsmin=.01*Jsmax;
	JsFinder JF(pot,Is,Rsh);
	Js_iso = math::findRoot(JF,Jsmin,Jsmax,relToler);
	b_iso = Rsh/Is.g(Js_iso);
	TM = Isochrone(Js_iso,b_iso);
	printf("Js: %6.3f, b: %6.3f, Delta: %6.3f, Rsh: %6.3f, freqScale: %6.3f\n",
	       Js_iso,b_iso,Delta,Rsh,freqScale);
#ifdef PLT
	const int ns=40;
	double F[2], Js[ns], f_apo_peri[ns], F_apo_peri[ns];
	for(int i=0; i<ns; i++){
		Js[i]=Jsmin+i*(Jsmax-Jsmin)/(double)(ns-1);
		double b=Rsh/Is.g(Js[i]);
		double e;
		f_apo_peri[i]=Is.f(b,Js[i],e);
		double c=Is.cob(Js[i])*b;
		for(int k=-1; k<2; k+=2){
			double u=1+k*e;
			double Phi, r=c*sqrt(u*(u+2*b/c));
			if(k==-1) printf("(%f ",r); else printf("%f ",r);
			coord::PosCyl Rz(r,0,0); coord::GradCyl grad;
			pot.eval(Rz,&Phi,&grad);
			F[(k+1)/2]=grad.dR;
		}
		F_apo_peri[i]=F[1]/F[0];
		printf("%f) ",b);
	}
	mgo::plt pl;
	pl.new_plot(Jsmin,Jsmax,0,1,"J\\ds","F\\da/F\\dp");
	pl.connect(Js,f_apo_peri,ns);
	pl.setcolour("red");
	pl.connect(Js,F_apo_peri,ns);
	pl.relocate(Jsmin,.1);pl.label(" true \\gF");
	pl.grend();
#endif
}
BaseTorus TorusGenerator::fitBaseTorus(const Actions& J){
	double Jscale, tolerance=1e-8;//controls optimisation of the given Sn
	setConsts(J, Jscale);//choose Rsh, Delta, freqScale, Js_iso, b_iso
	std::vector<double> params;
	int nrmax = 3, nzmax = 4;// nzmax must be even
	double Hbar, Hdisp = 1e20;
	bool converged = false;
	int Loop = 0, MaxLoop = 3, maxNumIter=15;
	GenFncIndices indices = makeGridIndices(nrmax, nzmax);
	do{
		GenFncFit GFFc(indices, J); GFF = &GFFc;
		params.resize(indices.size());
		try{
			//printf("%d values, %d params\n",numValues(),numVars());
			int numIter = math::nonlinearMultiFit(*this, &params[0], tolerance, maxNumIter, &params[0]);
			std::cout << numIter << " iterations; " << indices.size() << " GF terms;\n";
			//GFF->print(params);
			Hdisp = computeHamiltonianDisp(params, Hbar);
			converged = (Hdisp < tol*freqScale*Jscale);
			printf("Hbar: %f, Hdisp: %g with target %g ",Hbar, Hdisp,tol*freqScale*Jscale);
			if(NANfrac == 0) printf("\n");
			else printf("Fraction %7.3f of H vals NANs\n",NANfrac);
		}
		catch(std::exception& e) {
			std::cout << "Exception in fitTorus: " << e.what() <<'\n';
		}
		if(converged) break;
		Loop++;
		if(Loop<MaxLoop) indices = GFFc.expand(params);
	} while (Loop<MaxLoop);
	if(!converged) printf("fitTorus failed to converge: Hdisp = %g\n",Hdisp);
	GenFncFit GFFc(indices, J); GFF = &GFFc;
#ifdef TESTIT
	test_it(J,params);
#endif
	Frequencies freqs;
	GenFncDerivs derivs;
	Hdisp = fitAngleMap(&params[0], Hbar, freqs, derivs);
	//printf("Hbar: %f, freqs: %f\n", Hbar, freqs.Omegaz);
	GenFnc G(indices, params, derivs);
	BaseTorus T(J,freqs,G,TM,cs,Hbar);
	return T;
}
Torus TorusGenerator::fitTorus(const Actions& J){
	BaseTorus BT(fitBaseTorus(J));
	return Torus(BT);
}
eTorus TorusGenerator::fiteTorus(const Actions& J, const potential::BasePotential* ePhi){
	BaseTorus T=fitBaseTorus(J);
	PerturbingHamiltonian pH(get_pH(T,64,true,ePhi));
	return eTorus(T,pH);
}
double TorusGenerator::getDelta(double& L,double& Xi){
	return interpD.value(L,Xi);
}
double TorusGenerator::getDelta(actions::Actions& J){
	double L = fabs(J.Jphi)+J.Jz, Xi=J.Jz/L;
	return interpD.value(L,Xi);
}
std::vector<Torus> TorusGenerator::constE(const double Jrmin, const Actions& Jstart, const int Nsteps){
	double fac=exp(-log(Jstart.Jr/Jrmin)/(Nsteps-1));
	Actions Jnext(Jstart);
	std::vector<Torus> Tgrd;
	Torus T(fitTorus(Jnext));
	Tgrd.push_back(T);
	for(int i=1; i<Nsteps; i++){
		double frat=T.freqs.Omegar/T.freqs.Omegaz;
		double dJr=(1-fac)*Jnext.Jr;
		Jnext.Jr-=dJr; Jnext.Jz+=dJr*frat;
		printf("Next Jr: %f\n",Jnext.Jr);
		Torus T1(fitTorus(Jnext));
		T1*=0.5; T1+=T*0.5;
		double frat1=T1.freqs.Omegar/T1.freqs.Omegaz;
		Jnext.Jz=Jnext.Jz+dJr*frat-dJr*frat1;
		T=fitTorus(Jnext);
		Tgrd.push_back(T);
	}
	return Tgrd;
}

}//namespace
