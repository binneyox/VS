#pragma once
#include "potential_base.h"
#include "math_core.h"
#include "math_specfunc.h"
#include "actions_base.h"
#include "coord.h"
#include <iostream>
#include <cassert>

namespace actions {

class EXP Isochrone {
	public:
		double Js, b;
		Isochrone(double _Js=1, double _b=1) :
		    Js(_Js), b(_b) {}
		coord::PosMomSph aa2pq(const ActionAngles& aa, Frequencies* freqs=NULL,
				       DerivAct<coord::Sph>* dJ=NULL, DerivAng<coord::Sph>* dA=NULL) const;
		coord::PosMomSph aa2pq(const Actions& J, const Angles& theta, Frequencies* freqs=NULL,
				       DerivAct<coord::Sph>* dJ=NULL, DerivAng<coord::Sph>* dA=NULL) const{
			return aa2pq(ActionAngles(J,theta), freqs, dJ, dA);
		}
		ActionAngles pq2aa(const coord::PosMomSph& rtheta, Frequencies* freqs=NULL) const;
		Isochrone& operator *= (const double a) {
			Js *= a; b*=a;
			return *this;
		}
		Isochrone& operator += (const Isochrone& I){
			Js += I.Js; b+= I.b;
			return *this;
		}
		const Isochrone operator * (const double a) {
			Isochrone I2(Js*a,b*a);
			return I2;
		}
		const Isochrone operator + (const Isochrone I){
			Isochrone I2(Js+I.Js,b+I.b);
			return I2;
		}
};

}//namespace