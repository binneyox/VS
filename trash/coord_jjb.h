#include "coord.h"
#include "units.h"
#include "math_core.h"
#include <math.h>

namespace coord {
/* Classes for sky coordinates and proper motions. Units are degrees
 * and mas/yr rather than internal units
 */
struct VelSky{
	double mul, mub;
	VelSky(double _mul,double _mub) : mul(_mul), mub(_mub){}
};
struct PosVelSky{
	PosSky pos;
	VelSky pm;
	PosVelSky(PosSky _pos,VelSky _pm) : pos(_pos), pm(_pm){}
};

/* Class to convert obs coords (l,b in deg, mu in mas/yr, s in kpc,
 * vlos in km/s) into 6d Galactocentric phase space coords in int units
 * and back. Sun's Galactocentric phase-space coords should be in Kpc &
 * kms
 * */
class EXP solarShifter{
	private:
		PosVelCar Sun;
		double from_mas_per_yr, torad, from_Kpc, from_kms;
	public:
		solarShifter(const units::InternalUnits intUnits, PosVelCar* _Sun=NULL);
		PosVelCar toCar(const PosSky pos,double s,const VelSky pm,double vlos);
		PosVelCyl toCyl(const PosSky pos,double s,const VelSky pm,double vlos);
		VelSky toPM(const PosVelCar pos,double& vlos);
};

/* Transformations between equatorial & Galactic coords (units degreea)
 */
EXP PosSky from_RAdec(double ra,double dec);
EXP PosSky to_RAdec(double l,double b);
EXP PosVelSky from_muRAdec(double ra,double dec,double mura,double mudec);
}//namespace