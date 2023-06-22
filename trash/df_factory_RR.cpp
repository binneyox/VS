#include "df_factory_RR.h"
#include "df_disk_RR.h"
#include "df_halo_RR.h"
#include "utils.h"
#include "utils_config.h"
#include <cassert>
#include <stdexcept>

namespace df {

DoublePowerLawParam parseDoublePowerLawParams(
					      const utils::KeyValueMap& kvmap,
					      const units::ExternalUnits& conv)
{
	DoublePowerLawParam par;
	par.norm      = kvmap.getDouble("norm")    * conv.massUnit;
	par.J0        = kvmap.getDouble("J0")      * conv.lengthUnit * conv.velocityUnit;
	par.Jcutoff   = kvmap.getDouble("Jcutoff") * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0     = kvmap.getDouble("Jphi0")   * conv.lengthUnit * conv.velocityUnit;
	par.slopeIn   = kvmap.getDouble("slopeIn",   par.slopeIn);
	par.slopeOut  = kvmap.getDouble("slopeOut",  par.slopeOut);
	par.steepness = kvmap.getDouble("steepness", par.steepness);
	par.coefJphiIn  = kvmap.getDouble("coefJphiIn",  par.coefJphiIn);
	par.coefJzIn  = kvmap.getDouble("coefJzIn",  par.coefJzIn);
	par.coefJphiOut = kvmap.getDouble("coefJphiOut", par.coefJphiOut);
	par.coefJzOut = kvmap.getDouble("coefJzOut", par.coefJzOut);
	par.rotFrac   = kvmap.getDouble("rotFrac",   par.rotFrac);
	par.cutoffStrength = kvmap.getDouble("cutoffStrength", par.cutoffStrength);
	return par;
}
SoftDoublePowerLawParam parseSoftDoublePowerLawParams(
					      const utils::KeyValueMap& kvmap,
					      const units::ExternalUnits& conv)
{
	SoftDoublePowerLawParam par;
	par.norm      = kvmap.getDouble("norm")    * conv.massUnit;
	par.J0        = kvmap.getDouble("J0")      * conv.lengthUnit * conv.velocityUnit;
	par.Jcutoff   = kvmap.getDouble("Jcutoff") * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0     = kvmap.getDouble("Jphi0")   * conv.lengthUnit * conv.velocityUnit;
	par.slopeIn   = kvmap.getDouble("slopeIn",   par.slopeIn);
	par.slopeOut  = kvmap.getDouble("slopeOut",  par.slopeOut);
	par.steepness = kvmap.getDouble("steepness", par.steepness);
	par.coefJphiIn  = kvmap.getDouble("coefJphiIn",  par.coefJphiIn);
	par.coefJzIn  = kvmap.getDouble("coefJzIn",  par.coefJzIn);
	par.coefJphiOut = kvmap.getDouble("coefJphiOut", par.coefJphiOut);
	par.coefJzOut = kvmap.getDouble("coefJzOut", par.coefJzOut);
	par.rotFrac   = kvmap.getDouble("rotFrac",   par.rotFrac);
	par.cutoffStrength = kvmap.getDouble("cutoffStrength", par.cutoffStrength);
	par.h0	      = kvmap.getDouble("h0", par.h0) * conv.lengthUnit * conv.velocityUnit;
	return par;
}

QuasiIsothermalParam parseQuasiIsothermalParams(
    const utils::KeyValueMap& kvmap,
    const units::ExternalUnits& conv)
{
    QuasiIsothermalParam par;
    par.Sigma0  = kvmap.getDouble("Sigma0")  * conv.massUnit / pow_2(conv.lengthUnit);
    par.Rdisk   = kvmap.getDouble("Rdisk")   * conv.lengthUnit;
    par.Hdisk   = kvmap.getDouble("Hdisk")   * conv.lengthUnit;
    par.sigmar0 = kvmap.getDouble("sigmar0") * conv.velocityUnit;
    par.sigmaz0 = kvmap.getDouble("sigmaz0") * conv.velocityUnit;
    par.sigmamin= kvmap.getDouble("sigmamin")* conv.velocityUnit;
    par.Rsigmar = kvmap.getDouble("Rsigmar") * conv.lengthUnit;
    par.Rsigmaz = kvmap.getDouble("Rsigmaz") * conv.lengthUnit;
    par.coefJr  = kvmap.getDouble("coefJr", par.coefJr);
    par.coefJz  = kvmap.getDouble("coefJz", par.coefJz);
    par.Jmin    = kvmap.getDouble("Jmin") * conv.lengthUnit * conv.velocityUnit;
    par.beta_r    = kvmap.getDouble("beta_r", par.beta_r);
    par.beta_z    = kvmap.getDouble("beta_z", par.beta_z);
    par.Tsfr    = kvmap.getDouble("Tsfr", par.Tsfr);  // dimensionless! in units of Hubble time (galaxy age)
    par.sigmabirth = kvmap.getDouble("sigmabirth", par.sigmabirth);  // dimensionless ratio
    return par;
}

ExponentialParam parseExponentialParams(
					const utils::KeyValueMap& kvmap,
					const units::ExternalUnits& conv)
{
	ExponentialParam par;
	par.norm   = kvmap.getDouble("norm")   * conv.massUnit;
	par.Jr0    = kvmap.getDouble("Jr0")    * conv.lengthUnit * conv.velocityUnit;
	par.Jz0    = kvmap.getDouble("Jz0")    * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0  = kvmap.getDouble("Jphi0")  * conv.lengthUnit * conv.velocityUnit;
	par.pr     = kvmap.getDouble("pr");
	par.pz     = kvmap.getDouble("pz");
	par.addJden= kvmap.getDouble("addJden")* conv.lengthUnit * conv.velocityUnit;
	par.addJvel= kvmap.getDouble("addJvel")* conv.lengthUnit * conv.velocityUnit;
	return par;
}

newExpParam parsenewExpParams(
			      const utils::KeyValueMap& kvmap,
			      const units::ExternalUnits& conv)
{
	newExpParam par;
	par.norm   = kvmap.getDouble("norm")   * conv.massUnit;
	par.Jr0    = kvmap.getDouble("Jr0")    * conv.lengthUnit * conv.velocityUnit;
	par.Jz0    = kvmap.getDouble("Jz0")    * conv.lengthUnit * conv.velocityUnit;
	par.Jcut   = kvmap.getDouble("Jcut")   * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0  = kvmap.getDouble("Jphi0")  * conv.lengthUnit * conv.velocityUnit;
	par.Js     = kvmap.getDouble("Js",2000)* conv.lengthUnit * conv.velocityUnit;
	par.power  = kvmap.getDouble("power");
	par.addJden= kvmap.getDouble("addJden")* conv.lengthUnit * conv.velocityUnit;
	par.addJvel= kvmap.getDouble("addJvel")* conv.lengthUnit * conv.velocityUnit;
	par.coefJr = kvmap.getDouble("coefJr", par.coefJr);
	par.coefJz = kvmap.getDouble("coefJz", par.coefJz);
	par.beta_r = kvmap.getDouble("beta_r", par.beta_r);
	par.beta_z = kvmap.getDouble("beta_z", par.beta_z);
	par.Tsfr   = kvmap.getDouble("Tsfr", par.Tsfr);  // dimensionless! in units of Hubble time (galaxy age)
	par.sigmabirth = kvmap.getDouble("sigmabirth", par.sigmabirth);  // dimensionless ratio
	return par;
}

basicExpParam parsebasicExpParams(
			      const utils::KeyValueMap& kvmap,
			      const units::ExternalUnits& conv)
{
	basicExpParam par;
	par.norm   = kvmap.getDouble("norm")   * conv.massUnit;
	par.Jr0    = kvmap.getDouble("Jr0")    * conv.lengthUnit * conv.velocityUnit;
	par.Jz0    = kvmap.getDouble("Jz0")    * conv.lengthUnit * conv.velocityUnit;
	par.Jcut   = kvmap.getDouble("Jcut")   * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0  = kvmap.getDouble("Jphi0")  * conv.lengthUnit * conv.velocityUnit;
	par.pr	   = kvmap.getDouble("pr");
	par.pz	   = kvmap.getDouble("pz");
	par.addJden= kvmap.getDouble("addJden")* conv.lengthUnit * conv.velocityUnit;
	par.addJvel= kvmap.getDouble("addJvel")* conv.lengthUnit * conv.velocityUnit;
	return par;
}

inline void checkNonzero(const potential::BasePotential* potential, const std::string& type)
{
    if(potential == NULL)
        throw std::invalid_argument("Need an instance of potential to initialize "+type+" DF");
}

PtrDistributionFunction createDistributionFunction(
    const utils::KeyValueMap& kvmap,
    const potential::BasePotential* potential,
    const units::ExternalUnits& converter)
{
    std::string type = kvmap.getString("type");
    if(utils::stringsEqual(type, "DoublePowerLaw")) {
	    return PtrDistributionFunction(new DoublePowerLaw(parseDoublePowerLawParams(kvmap, converter)));
    }
    if(utils::stringsEqual(type, "SoftDoublePowerLaw")) {
	    return PtrDistributionFunction(new SoftDoublePowerLaw(parseSoftDoublePowerLawParams(kvmap, converter)));
    }
    if(utils::stringsEqual(type, "Exponential")) {
	    return PtrDistributionFunction(new Exponential(parseExponentialParams(kvmap, converter)));
    }
    if(utils::stringsEqual(type, "newExp")) {
	    return PtrDistributionFunction(new newExp(parsenewExpParams(kvmap, converter)));
    }
    if(utils::stringsEqual(type, "basicExp")) {
	    return PtrDistributionFunction(new basicExp(parsebasicExpParams(kvmap, converter)));
    }
    else if(utils::stringsEqual(type, "QuasiIsothermal")) {
        checkNonzero(potential, type);
        return PtrDistributionFunction(new QuasiIsothermal(parseQuasiIsothermalParams(kvmap, converter),
            potential::Interpolator(*potential)));
    }
    else
        throw std::invalid_argument("Unknown type of distribution function");
}

}; // namespace
