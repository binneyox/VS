/* \file    example_splitExpModel.cpp
    \author  Eugene Vasiliev adapted by James Binney July 2019
    \date    2015-2018

    This example demonstrates the machinery for constructing multicomponent self-consistent models
    specified by distribution functions in terms of actions.
    We create a six-component galaxy with disk, bulge and halo components defined by their DFs,
    and a static density profile of gas disk. The thin disk is split by age
    Then we perform several iterations of recomputing the density profiles of components from their DFs
    and recomputing the total potential.
    Finally, we create N-body representations of all mass components: dark matter halo,
    stars (bulge, thin and thick disks and stellar halo combined), and gas disk.
	
	edited 2021 for making queiroz plots
	needs to output a catalogue of stars for each component having e and zmax and j for each star in a line of sight.
	also needs to output integrted DFintegrandLOS to give values for each component for each point on sky used to 
	find proportions of each component for given (l,b) so that chemistry may be assigned.
*/
#include "galaxymodel_base.h"
#include "galaxymodel_selfconsistent.h"
#include "galaxymodel_velocitysampler.h"
#include "df_factory.h"
#include "potential_composite.h"
#include "potential_factory.h"
#include "potential_multipole.h"
#include "potential_utils.h"
#include "particles_io.h"
#include "math_core.h"
#include "math_spline.h"
#include "units.h"
#include "utils.h"
#include "utils_config.h"
#include "actions_staeckel.h"
//#include "press.h"
#include "actions_torus.h"
#include "orbit.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sys/stat.h>
#include <direct.h>
#include <algorithm>
#include <random>
using potential::PtrDensity;
using potential::PtrPotential;

inline bool exists(const std::string& filename){
	struct stat buffer;
	return (stat (filename.c_str(), &buffer) == 0);
}
// define internal unit system - arbitrary numbers here! the result should not depend on their choice
const units::InternalUnits intUnits(2.7183*units::Kpc, 3.1416*units::Myr);

// define external unit system describing the data (including the parameters in INI file)
const units::ExternalUnits extUnits(intUnits, 1.*units::Kpc, 1.*units::kms, 1.*units::Msun);

// used for outputting the velocity distribution (the value is read from the ini file)
double solarRadius = NAN;

// various auxiliary functions for printing out information are non-essential
// for the modelling itself; the essential workflow is contained in main()

//function for producing data for a scatter plot of Vphi against Fe/H for the thin disk as in Largarde et.al 2021
void writeLagarde(const std::string& filename,const galaxymodel::GalaxyModel& model, const std::vector<galaxymodel::GalaxyModel> modelComponents,
				 const  galaxymodel::SelfConsistentModel& SCM)
{
	std::ofstream strm(filename.c_str());
	std::cout << "Writing Lagarde data\n";
	int nc = model.distrFunc.numValues();//number of components
	std::vector<double> dens(nc),Vphi(nc);
	std::vector<coord::Vel2Cyl> sigmas(nc);

	int N=4265;//the number of stars at each location
	obs::solarShifter shift(intUnits);
	coord::PosCyl pos=coord::toPosCyl(shift.xyz());				//probably want int units
	galaxymodel::computeMoments(model, pos, &dens[0], &Vphi[0],		//finds the density at each point of each component
					&sigmas[0], intUnits,NULL, NULL, NULL, NULL, true, 1e-3,1e5,NULL,NULL,&shift); //added for sKpc
	//now sample velocities for each component
	double sum=0;
	for(int k=0;k<3;k++){		//only the thin disk
		sum+=dens[k];
	}
	for(int i=0;i<3;i++){//for each component in the thin disk
		double num=N*dens[i]/sum;		//not an int! may impact sample velocity? does it round up or down?
		std::vector<coord::VelCyl> vels = galaxymodel::sampleVelocity(modelComponents[i], pos,num);
		//for each sampled star
		for(int k=0; k<vels.size();k++){
			coord::PosVelCyl xv(pos, vels[k]);
			//calculate actions for the star
			actions::Actions J = SCM.actionFinder->actions(xv);
			//print results for each star 
			strm<< 
				utils::toString(pos.R * intUnits.to_Kpc)+"\t"+			
				utils::toString(pos.z* intUnits.to_Kpc)+"\t"+
				utils::toString(J.Jr * intUnits.to_Kpc_kms)+"\t"+
				utils::toString(J.Jz * intUnits.to_Kpc_kms)+"\t"+				
				utils::toString(J.Jphi * intUnits.to_Kpc_kms)+"\t"+	
				utils::toString(xv.vphi * intUnits.to_kms)+"\t"+	
				utils::toString(i)+"\t"+
				"\n";
		}
	}	
}

//funciton for producing hayden plots
void writeHayden(const std::string& filename,const galaxymodel::GalaxyModel& model, const std::vector<galaxymodel::GalaxyModel> modelComponents,
				 const  galaxymodel::SelfConsistentModel& SCM,const double bright = NULL, const double faint = NULL, obs::solarShifter* sun = NULL)
{
	std::ofstream strm(filename.c_str());
	std::cout << "Writing Hayden data\n";
	std::vector<double> Rs,zs;
	double z=0.2;
	for(int i=0;i<3;i++){
		if(i==1) z=0.7; else if(i==2) z=1.4;
		for(double R=4; R<=14; R+=2){
			Rs.push_back(R * intUnits.from_Kpc);
			zs.push_back(z * intUnits.from_Kpc);
		}
	}
	int nh = Rs.size();//number of points
	int nc = model.distrFunc.numValues();//number of components
	int n = nh*nc;
	std::vector<double> dens(n),Vphi(n);
	std::vector<coord::Vel2Cyl> sigmas(n);
/*#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif*/
	int N=5000;//the number of stars at each location
	for(int ih=0; ih<nh; ih++) {//for each location
		int ic=ih*nc;
		coord::PosCyl pos(Rs[ih],zs[ih],0);
		galaxymodel::computeMoments(model, pos, &dens[ic], &Vphi[ic],		//finds the density at each point of each component
			       &sigmas[ic],intUnits, NULL, NULL, NULL, NULL, true, 1e-3,1e5,NULL,NULL,sun); //added for sKpc
		//now sample velocities for each component
		double sum=0;
		for(int k=0;k<nc;k++){
			sum+=dens[ic+k];
		}
		for(int i=0;i<nc;i++){//for each component
			double num=N*dens[ic+i]/sum;		//not an int! may impact sample velocity? does it round up or down?
			std::vector<coord::VelCyl> vels = galaxymodel::sampleVelocity(modelComponents[i], pos,num);
			//for each sampled star
			for(int k=0; k<vels.size();k++){
				coord::PosVelCyl xv(pos, vels[k]);
				//calculate actions for the star
				actions::Actions J = SCM.actionFinder->actions(xv);
				//print results for each star 
				strm<< 
					utils::toString(pos.R * intUnits.to_Kpc)+"\t"+			
					utils::toString(pos.z* intUnits.to_Kpc)+"\t"+
					utils::toString(J.Jr * intUnits.to_Kpc_kms)+"\t"+
					utils::toString(J.Jz * intUnits.to_Kpc_kms)+"\t"+				
					utils::toString(J.Jphi * intUnits.to_Kpc_kms)+"\t"+	
					utils::toString(i)+"\t"+
					"\n";
			}
		}	
	}
}

///write PosVelCyl for each star in sample
void writePosVel(const std::string& directory, const std::string& l, const std::string& b, const std::vector<std::pair<std::vector<coord::PosVelCyl>, int>> results,
				const double bright = NULL, const double faint = NULL, obs::solarShifter* sun = NULL)
{
	//write xv to a file
	std::cout<<"printing positions and velocities\n";
	std::string filename = directory + "xvCyl";
	std::ofstream strm(filename.c_str());
	strm<<bright<<"\t"<<faint<<"\n";		//print bright and faint apparent mag values
	std::string fname = directory + "xvCar";
	std::ofstream strmCar(fname.c_str());
	strmCar<<bright<<"\t"<<faint<<"\n";		//print bright and faint apparent mag values
	//for each component
	for(int j=0; j<results.size();j++){
		//for each sampled star
		int N = results[j].first.size();		//do one component first
		for(int i=0; i<N;i++){
			coord::PosVelCyl xv = results[j].first[i];
			coord::PosCyl x(xv.R, xv.z, xv.phi);
			coord::PosVelCar xvCar(coord::toPosVelCar(xv));
			strm<< 
				utils::toString(sun->s(x))+"\t"+
				utils::toString(xv.R * intUnits.to_Kpc)+"\t"+
				utils::toString(xv.z * intUnits.to_Kpc)+"\t"+
				utils::toString(xv.phi)+"\t"+						
				utils::toString(xv.vR * intUnits.to_kms)+"\t"+
				utils::toString(xv.vz * intUnits.to_kms)+"\t"+
				utils::toString(xv.vphi * intUnits.to_kms)+"\t"+
				utils::toString(results[j].second)+"\t"+
				"\n";
			strmCar << 
				utils::toString(sun->s(x))+"\t"+
				utils::toString(xvCar.x * intUnits.to_Kpc)+"\t"+
				utils::toString(xvCar.y * intUnits.to_Kpc)+"\t"+
				utils::toString(xvCar.z * intUnits.to_Kpc)+"\t"+						
				utils::toString(xvCar.vx * intUnits.to_kms)+"\t"+
				utils::toString(xvCar.vy * intUnits.to_kms)+"\t"+
				utils::toString(xvCar.vz * intUnits.to_kms)+"\t"+
				utils::toString(results[j].second)+"\t"+
				"\n";
		}
	}
}

///write JEZ for each star in sample
void writeJEZ(const std::string& directory, const std::string& l, const std::string& b, const std::vector<std::pair<std::vector<coord::PosVelCyl>, int>> results,
				const galaxymodel::SelfConsistentModel& model, const double bright = NULL, const double faint = NULL)
{	//write actions, eccentricity and Zmax to a file
	std::cout<<"printing actions, eccentricity and Zmax...";
	std::string filename = directory + "JEZ("+l+","+b+")";
	std::ofstream strm(filename.c_str());
	//print position in degrees
	strm<<l+"\t"+
		b+"\t"+"\n";
	strm<<bright<<"\t"<<faint<<"\n";		//print bright and faint apparent mag values
	//for each component
	for(int k=0;k<results.size();k++){
		int N = results[k].first.size();
		//for each sampled star
		for(int i=0; i<N;i++){
			coord::PosVelCyl xv = results[k].first[i];
			//calculate actions for the star
			actions::Actions J = model.actionFinder->actions(xv);		
			double Zmax=0, Rmax=0, Rmin=0;
			//integrate orbits if they have problematic actions that torus mapper cannot handle
			if(J.Jr<0.05*J.Jz)					//this may not be the only problematic J value
			{ 	
				//find total_time
				double E = potential::totalEnergy(*model.totalPotential, xv);
				double total_time = 25*potential::T_circ(*model.totalPotential, E);		//increasing this will improve estimates of Z/Rmax
				//integrate the orbit
				double timestep=0.00005*total_time;
				std::vector<std::pair<coord::PosVelCyl, double> > traj =
					orbit::integrateTraj(xv, total_time, timestep, *model.totalPotential);
				//find Zmax
				std::vector<double> Zp;
				for(int j=0; j<traj.size()-1;  j++){
					if(traj[j].first.vz*traj[j+1].first.vz<0){				//check what happens with zeros
						double z0 = traj[j].first.z, z1 = traj[j+1].first.z;
						double vz0 = traj[j].first.vz, vz1 = traj[j+1].first.vz;
						//using cubic fitting for z=z0+at+bt^2+ct^3
						double a = vz0, b = 3*(z1-z0)-vz1-2*vz0, c = 2*(z0-z1)+vz1+vz0;
						double tp = (sqrt(pow(b,2)-3*a*c)-b)/(3*c);
						double zp = z0 + a*tp + b*pow(tp,2) + c*pow(tp,3);		//value of z for this maximum
						Zp.push_back(abs(zp));
					}
				}
				//find Zmax as the maximum value in Zp
				Zmax = *std::max_element(Zp.begin(), Zp.end());					
				
				//find Rmax and Rmin			
				std::vector<double> Rp;
				for(int j=0; j<traj.size()-1;  j++){
					if(traj[j].first.vR*traj[j+1].first.vR<0){				
						double R0 = traj[j].first.R, R1 = traj[j+1].first.R;
						double vR0 = traj[j].first.vR, vR1 = traj[j+1].first.vR;
						//using cubic fitting for R=R0+at+bt^2+ct^3
						double a = vR0, b = 3*(R1-R0)-vR1-2*vR0, c = 2*(R0-R1)+vR1+vR0;
						double tp = (sqrt(pow(b,2)-3*a*c)-b)/(3*c);
						double rp = R0 + a*tp + b*pow(tp,2) + c*pow(tp,3);		//value of R for this maximum/minimum
						Rp.push_back(abs(rp));
					}
				}
				//find Rmax and Rmin as the maximum and minimum values in Rp
				Rmax = *std::max_element(Rp.begin(), Rp.end());					
				Rmin = *std::min_element(Rp.begin(), Rp.end());
			}
			//if actions are ok then use torus to find Rmax/min and Zmax
			else{
				actions::ActionMapperTorus mapper(*model.totalPotential,J);	//<--error occurs in here if J values are problematic. !!note torus.h has been edited!!
				actions::Angles theta(M_PI,0.5*M_PI,0);
				actions::ActionAngles actAng(J,theta);
				Zmax = mapper.map(actAng).z;			
				actAng.thetaz=0;
				Rmax = mapper.map(actAng).R;					
				actAng.thetar=0; actAng.thetaz=.5*M_PI;
				Rmin = mapper.map(actAng).R;
			}
			double e = (Rmax-Rmin)/(Rmax+Rmin);
			//print info for each star
			strm<< 
				utils::toString(J.Jr * intUnits.to_Kpc_kms)+"\t"+
				utils::toString(J.Jz * intUnits.to_Kpc_kms)+"\t"+				
				utils::toString(J.Jphi * intUnits.to_Kpc_kms)+"\t"+	
				utils::toString(e)+"\t"+			
				utils::toString(Zmax * intUnits.to_Kpc)+"\t"+
				utils::toString(results[k].second)+"\t"+
				"\n";
		}
	}
	std::cout<<"done\n";
}

///write integrals for each component
std::vector<double> writeIntegral(const std::string& directory, const std::string& l, const std::string& b, 
				const obs::los* los, const galaxymodel::GalaxyModel& model,const double bright=NULL,
				const double faint=NULL, obs::solarShifter* sun = NULL)
{	
	std::cout<<"printing integrals along line of sight\n";
	std::string filename = directory + "Integrals("+l+","+b+")";
	std::ofstream strm(filename.c_str());
	//print position in degrees
	strm<<l+"\t"+
		b+"\t"+"\n";
	strm<<bright<<"\t"<<faint<<"\n";		//print bright and faint apparent mag values
	int nc = model.distrFunc.numValues();
	std::cout<< "we have "<<nc<<" components\n";
	std::vector<double> dens(nc),Vphi(nc);
	std::vector<coord::Vel2Cyl> sigmas(nc);
	int att=1;
	galaxymodel::computeMomentsLOS(model, los, &dens[0], &Vphi[0],
			       &sigmas[0], NULL, NULL, NULL, true, 1e-3, 1e5, bright, faint);
	std::cout<< "computed moments\n";
	double sum=0;
	for(int comp=0;comp<nc;comp++) sum+=dens[comp];
	for(int i = 0; i<nc; i++){						
		strm<<
			utils::toString(dens[i] * intUnits.to_Msun_per_Kpc2) + "\t";			
	}	
	return dens;
}
/// report progress after an iteration
void printoutInfo(const galaxymodel::SelfConsistentModel& model, const std::string& dir)
{
	const potential::BaseDensity& compDisk = *model.components[0]->getDensity();
	const potential::BaseDensity& compBulge= *model.components[1]->getDensity();
	const potential::BaseDensity& compHalo = *model.components[2]->getDensity();
	coord::PosCyl pt0(solarRadius * intUnits.from_Kpc, 0, 0);
	coord::PosCyl pt1(solarRadius * intUnits.from_Kpc, 1 * intUnits.from_Kpc, 0);
	std::cout <<
			"Disk total mass="      << (compDisk.totalMass()  * intUnits.to_Msun) << " Msun"
			", rho(Rsolar,z=0)="    << (compDisk.density(pt0) * intUnits.to_Msun_per_pc3) <<
			", rho(Rsolar,z=1kpc)=" << (compDisk.density(pt1) * intUnits.to_Msun_per_pc3) << " Msun/pc^3\n"
			"Halo total mass="      << (compHalo.totalMass()  * intUnits.to_Msun) << " Msun"
			", rho(Rsolar,z=0)="    << (compHalo.density(pt0) * intUnits.to_Msun_per_pc3) <<
			", rho(Rsolar,z=1kpc)=" << (compHalo.density(pt1) * intUnits.to_Msun_per_pc3) << " Msun/pc^3\n"
			"Potential at origin=-("<<
			(sqrt(-model.totalPotential->value(coord::PosCyl(0,0,0))) * intUnits.to_kms) << " km/s)^2"
			", total mass=" << (model.totalPotential->totalMass() * intUnits.to_Msun) << " Msun\n";
	writeDensity(dir + "dens_disk", compDisk, extUnits);
	writeDensity(dir + "dens_bulge", compBulge, extUnits);
	writeDensity(dir + "dens_halo", compHalo, extUnits);
	writePotential(dir + "potential", *model.totalPotential, extUnits);
	std::vector<PtrPotential> potentials(3);
	potentials[0] = dynamic_cast<const potential::CompositeCyl&>(*model.totalPotential).component(1);
	potentials[1] = potential::Multipole::create(compBulge, /*lmax*/6, /*mmax*/0, /*gridsize*/25);
	potentials[2] = potential::Multipole::create(compHalo,  /*lmax*/6, /*mmax*/0, /*gridsize*/25);
}

/// perform one iteration of the model
void doIteration(galaxymodel::SelfConsistentModel& model,const int iterationIndex, const std::string &dir)
{
    std::cout << "Starting iteration #" << iterationIndex << "\n";
    bool error=false;
    try {
        doIteration(model);
    }
    catch(std::exception& ex) {
        error=true;  // report the error and allow to save the results of the last iteration
        std::cout << "==== Exception occurred: \n" << ex.what();
    }
    printoutInfo(model, dir);
    if(error)
        exit(1);  // abort in case of problems
}


int prog(int narg,char **args)
{
	if(narg!=5){
		printf("You must specify:\n restart modelling? true/false\n read mode? true/false\n if read mode then l and b \n else then bright and faint\n"); return 0;
	}
	std::time_t duration = std::time(NULL); 
	bool newModel, lowZ=false;
	int tmp;
	sscanf(args[1],"%d",&tmp); newModel=tmp;
	const std::string dir = "C:/MinGW/msys/1.0/home/Thompson/agama/2022/";
    // read parameters from the INI file
	const std::string iniFileName = dir + "m.ini";
	utils::ConfigFile ini(iniFileName);
	utils::KeyValueMap
		    iniPotenThinDisk = ini.findSection("Potential thin disk"),
    iniPotenThickDisk= ini.findSection("Potential thick disk"),
    iniPotenGasDisk  = ini.findSection("Potential gas disk"),
    iniPotenBulge    = ini.findSection("Potential bulge"),
    iniPotenDarkHalo = ini.findSection("Potential dark halo"),
	iniDFyoungDisk    = ini.findSection("DF young disk"),
	iniDFmiddleDisk    = ini.findSection("DF middle disk"),
	iniDFoldDisk    = ini.findSection("DF old disk"),
	iniDFhighADisk   = ini.findSection("DF highA disk"),
	iniDFhighAlowZ   = ini.findSection("DF highAlowZ disk"),
        iniDFStellarHalo = ini.findSection("DF stellar halo"),
        iniDFBulge       = ini.findSection("DF bulge"),
	iniDFDarkHalo    = ini.findSection("DF dark halo"),
        iniSCMDisk       = ini.findSection("SelfConsistentModel disk"),
        iniSCMBulge      = ini.findSection("SelfConsistentModel bulge"),
        iniSCMHalo       = ini.findSection("SelfConsistentModel halo"),
        iniSCM           = ini.findSection("SelfConsistentModel");
    if(!iniSCM.contains("rminSph")) {  // most likely file doesn't exist
        std::cout << "Invalid INI file " << iniFileName << "\n";
        return -1;
    }
	std::cout<< "ini file read\n";
    if(!iniDFhighAlowZ.contains("mass")) lowZ = false;
    else lowZ = true;
    solarRadius = ini.findSection("Data").getDouble("SolarRadius", solarRadius);
	std::cout<<solarRadius<<"\n";

    // set up parameters of the entire Self-Consistent Model
    galaxymodel::SelfConsistentModel model;
    model.rminSph         = iniSCM.getDouble("rminSph") * extUnits.lengthUnit;
    model.rmaxSph         = iniSCM.getDouble("rmaxSph") * extUnits.lengthUnit;
    model.sizeRadialSph   = iniSCM.getInt("sizeRadialSph");
    model.lmaxAngularSph  = iniSCM.getInt("lmaxAngularSph");
    model.RminCyl         = iniSCM.getDouble("RminCyl") * extUnits.lengthUnit;
    model.RmaxCyl         = iniSCM.getDouble("RmaxCyl") * extUnits.lengthUnit;
    model.zminCyl         = iniSCM.getDouble("zminCyl") * extUnits.lengthUnit;
    model.zmaxCyl         = iniSCM.getDouble("zmaxCyl") * extUnits.lengthUnit;
    model.sizeRadialCyl   = iniSCM.getInt("sizeRadialCyl");
    model.sizeVerticalCyl = iniSCM.getInt("sizeVerticalCyl");
    model.useActionInterpolation = iniSCM.getBool("useActionInterpolation");
	

    // initialize density profiles of various components
    std::vector<PtrDensity> densityStellarDisk(2);
    PtrDensity densityBulge    = potential::createDensity(iniPotenBulge,    extUnits);
    PtrDensity densityDarkHalo = potential::createDensity(iniPotenDarkHalo, extUnits);
    densityStellarDisk[0]      = potential::createDensity(iniPotenThinDisk, extUnits);
    densityStellarDisk[1]      = potential::createDensity(iniPotenThickDisk,extUnits);
    PtrDensity densityGasDisk  = potential::createDensity(iniPotenGasDisk,  extUnits);

    // add components to SCM - at first, all of them are static density profiles
    model.components.push_back(galaxymodel::PtrComponent(
	    new galaxymodel::ComponentStatic(PtrDensity(
	    new potential::CompositeDensity(densityStellarDisk)), true)));
    model.components.push_back(galaxymodel::PtrComponent(
	    new galaxymodel::ComponentStatic(densityBulge, false)));
    model.components.push_back(galaxymodel::PtrComponent(
	    new galaxymodel::ComponentStatic(densityDarkHalo, false)));
    model.components.push_back(galaxymodel::PtrComponent(
	    new galaxymodel::ComponentStatic(densityGasDisk, true)));

    // initialize total potential of the model (first guess)
    updateTotalPotential(model);
    printoutInfo(model, "init");

    std::cout << "**** STARTING MODELLING ****\nInitial masses of density components: "
		    "Mdisk="  << (model.components[0]->getDensity()->totalMass() * intUnits.to_Msun) << " Msun, "
		    "Mbulge=" << (densityBulge->totalMass() * intUnits.to_Msun) << " Msun, "
		    "Mhalo="  << (densityDarkHalo->totalMass() * intUnits.to_Msun) << " Msun, "
		    "Mgas="   << (densityGasDisk->totalMass() * intUnits.to_Msun) << " Msun\n";

    // create the dark halo DF
    df::PtrDistributionFunction dfHalo = df::createDistributionFunction(
	    iniDFDarkHalo, model.totalPotential.get(), NULL, extUnits);
    printf("dark halo DF created\n");
    // same for the bulge
    df::PtrDistributionFunction dfBulge = df::createDistributionFunction(
	    iniDFBulge, model.totalPotential.get(), NULL, extUnits);
    printf("bulge DF created..");
    // same for the stellar components (thin/thick disks and stellar halo)
    std::vector<df::PtrDistributionFunction> dfStellarArray;
    dfStellarArray.push_back(df::createDistributionFunction(
	    iniDFyoungDisk, model.totalPotential.get(), NULL, extUnits));
    printf("young disc DF created..");
    dfStellarArray.push_back(df::createDistributionFunction(
	    iniDFmiddleDisk, model.totalPotential.get(), NULL, extUnits));
    printf("middle disc DF created..");
    dfStellarArray.push_back(df::createDistributionFunction(
	    iniDFoldDisk, model.totalPotential.get(), NULL, extUnits));
    printf("old disc DF created..");
    dfStellarArray.push_back(df::createDistributionFunction(
	    iniDFhighADisk, model.totalPotential.get(), NULL, extUnits));
    printf("highA disc DF created\n");
    if(lowZ){
	    dfStellarArray.push_back(df::createDistributionFunction(
		    iniDFhighAlowZ, model.totalPotential.get(), NULL, extUnits));
	    printf("HighAlowZ DF created\n");
    }
    dfStellarArray.push_back(df::createDistributionFunction(
	    iniDFStellarHalo, model.totalPotential.get(), NULL, extUnits));
    printf("Stellar halo DF created\n");
// composite DF of all stellar components except the bulge
    df::PtrDistributionFunction dfStellar(new df::CompositeDF(dfStellarArray));
//now create an all-star DF
    std::vector<df::PtrDistributionFunction> dfAllStarArray;
    for(int i=0; i<(int)dfStellarArray.size(); i++) dfAllStarArray.push_back(dfStellarArray[i]);
    dfAllStarArray.push_back(dfBulge);
    df::PtrDistributionFunction dfAllStar(new df::CompositeDF(dfAllStarArray));
//assemble total DF
    std::vector<df::PtrDistributionFunction> dfAllArray;
    dfAllArray.push_back(dfHalo);
    for(int i=0; i<(int)dfAllStarArray.size(); i++) dfAllArray.push_back(dfAllStarArray[i]);
    df::PtrDistributionFunction dfAll(new df::CompositeDF(dfAllArray));

    int ncmp=0;
    // replace the static disk density component of SCM with a DF-based disk component
    model.components[ncmp] = galaxymodel::PtrComponent(
        new galaxymodel::ComponentWithDisklikeDF(dfStellar, PtrDensity(),
        iniSCMDisk.getInt("mmaxAngularCyl"),
        iniSCMDisk.getInt("sizeRadialCyl"),
        iniSCMDisk.getDouble("RminCyl") * extUnits.lengthUnit,
        iniSCMDisk.getDouble("RmaxCyl") * extUnits.lengthUnit,
        iniSCMDisk.getInt("sizeVerticalCyl"),
        iniSCMDisk.getDouble("zminCyl") * extUnits.lengthUnit,
	    iniSCMDisk.getDouble("zmaxCyl") * extUnits.lengthUnit));
    // same for the bulge
    ncmp++;
    model.components[ncmp] = galaxymodel::PtrComponent(
        new galaxymodel::ComponentWithSpheroidalDF(dfBulge, potential::PtrDensity(),
        iniSCMBulge.getInt("lmaxAngularSph"),
        iniSCMBulge.getInt("mmaxAngularSph"),
        iniSCMBulge.getInt("sizeRadialSph"),
        iniSCMBulge.getDouble("rminSph") * extUnits.lengthUnit,
	    iniSCMBulge.getDouble("rmaxSph") * extUnits.lengthUnit));
    // same for the dark halo
    ncmp++;
    model.components[ncmp] = galaxymodel::PtrComponent(
        new galaxymodel::ComponentWithSpheroidalDF(dfHalo, potential::PtrDensity(),
        iniSCMHalo.getInt("lmaxAngularSph"),
        iniSCMHalo.getInt("mmaxAngularSph"),
        iniSCMHalo.getInt("sizeRadialSph"),
        iniSCMHalo.getDouble("rminSph") * extUnits.lengthUnit,
	    iniSCMHalo.getDouble("rmaxSph") * extUnits.lengthUnit));


    // we can compute the masses even though we don't know the density profile yet
    std::cout <<
        "Masses of DF components:"
		    "\n Mdisk + Mstel.halo=" << (dfStellar->totalMass() * intUnits.to_Msun) <<
		    " Msun\n (Myoung=" << (dfStellarArray[0]->totalMass() * intUnits.to_Msun) <<
		    ", Mmiddle="     << (dfStellarArray[1]->totalMass() * intUnits.to_Msun) <<
		    ", Mold="     << (dfStellarArray[2]->totalMass() * intUnits.to_Msun) <<
		    ", MHighA="     << (dfStellarArray[3]->totalMass() * intUnits.to_Msun);
    if(lowZ) std::cout <<
		    ", MHighAlowZ disk=" << (dfStellarArray[4]->totalMass() * intUnits.to_Msun) <<
		    ", Mstel.halo=" << (dfStellarArray[5]->totalMass() * intUnits.to_Msun) <<
		    ")\n Mbulge="    << (dfBulge->totalMass() * intUnits.to_Msun) << " Msun"
		    "; Mdark="      << (dfHalo ->totalMass() * intUnits.to_Msun) << " Msun\n";
    else  std::cout <<
		    ", Mstel.halo=" << (dfStellarArray[4]->totalMass() * intUnits.to_Msun) <<
		    ")\n Mbulge="    << (dfBulge->totalMass() * intUnits.to_Msun) << " Msun"
		    "; Mdark="      << (dfHalo ->totalMass() * intUnits.to_Msun) << " Msun\n";
  
    if(newModel){
	    if(utils::fileExists(dir + "potential"))
		    model.totalPotential = potential::readPotential(dir + "potential", extUnits);
	    // update the action finder
	    std::cout << "Updating action finder..."<<std::flush;
	    model.actionFinder.reset(new actions::ActionFinderAxisymFudge(model.totalPotential,
		    model.useActionInterpolation));
	    std::cout << "done"<<std::endl;
	    std::ofstream strm0(dir + "masses");
	    strm0 << "\n Mdisk + Mstel.halo=" << (dfStellar->totalMass() * intUnits.to_Msun) <<
			    " Msun\n (Myoung=" << (dfStellarArray[0]->totalMass() * intUnits.to_Msun) <<
			    ", Mmiddle="     << (dfStellarArray[1]->totalMass() * intUnits.to_Msun) <<
			    ", Mold="     << (dfStellarArray[2]->totalMass() * intUnits.to_Msun) <<
			    ", MHighA="     << (dfStellarArray[3]->totalMass() * intUnits.to_Msun) <<
			    ", Mstel.halo=" << (dfStellarArray[4]->totalMass() * intUnits.to_Msun) <<
			    ")\n Mbulge="    << (dfBulge->totalMass() * intUnits.to_Msun) << " Msun"
			    "; Mdark="      << (dfHalo ->totalMass() * intUnits.to_Msun) << " Msun\n";
	    strm0.close();
	    std::ofstream strm(dir + "DF.pars");


	    strm.close();// if(solarRadius>0) return 0;
    // do a few more iterations to obtain the self-consistent density profile for both disks
	    for(int iteration=1; iteration<=4; iteration++)
		    doIteration(model, iteration, dir);
	    duration -= std::time(NULL);
	    std::cout << -duration << " secs to build\n";
    } else {
	    model.totalPotential = potential::readPotential(dir + "potential", extUnits);
	    // update the action finder
	    std::cout << "Updating action finder..."<<std::flush;
	    model.actionFinder.reset(new actions::ActionFinderAxisymFudge(model.totalPotential,
		    model.useActionInterpolation));
	    std::cout << "done"<<std::endl;
	    }
    // output various profiles
    galaxymodel::GalaxyModel modelStars(*model.totalPotential, *model.actionFinder, *dfAllStar);
    galaxymodel::GalaxyModel modelAll(*model.totalPotential, *model.actionFinder, *dfAll);
	std::vector<galaxymodel::GalaxyModel> modelComponents; 
	for(int j=0;j<dfAllStarArray.size();j++){
		modelComponents.push_back(galaxymodel::GalaxyModel(*model.totalPotential, *model.actionFinder, *dfAllStarArray[j]));
	}
    std::cout << "Sampling the model\n";
	
	////POINT SAMPLING
	/*to produce hayden plot you need to sample at the specified points, the correct numbers from each component given by the density ratios
	at that point. */
	obs::solarShifter shifter(intUnits);
	coord::PosCyl sun = coord::toPosCyl(shifter.xyz());
	coord::GradCyl deriv;
	coord::HessCyl deriv2;
	double vc2;
	model.totalPotential->eval(sun, &vc2, &deriv, &deriv2);
	vc2 = sun.R*deriv.dR;
	std::cout<<"vc2="<<sqrt(vc2)*intUnits.to_kms <<"\n";
	//writeHayden("Hayden", modelStars, modelComponents, model,NULL, NULL, &shifter);		//account for bright/faint/dust??
	//writeLagarde("Lagarde", modelStars, modelComponents, model);
	
	////LINE OF SIGHT SAMPLING
	//if read mode then extract (l,b) and N for each line of sight from a file
	int mode, M;
	sscanf(args[2],"%i",&M); mode=M;
	if(mode==1){
		FILE* ifile;
		//read the first file to be mixed
		fopen_s(&ifile,"lines_of_sight/binnedCoords.dat","r");
		int i=0; 	//i will tell you the number of lines of sight for use in for loop
		double dat[2], glon[100], glat[100]; 
		int datN[1], Nstar[100];
		while(!feof(ifile)){
			if(3!=fscanf(ifile,"%lf %lf %i",dat,dat+1,datN)) break;
			glon[i]=dat[0];
			glat[i]=dat[1];
			Nstar[i]=datN[0];
			i++;
		}
		fclose(ifile);
		//read bright and faint as command line arguments
		double readbright, readfaint, bright, faint;
		sscanf(args[3],"%lf",&readbright); bright=readbright;
		sscanf(args[4],"%lf",&readfaint); faint=readfaint;
		std::cout<<"sampling with bright="<<bright<<" and faint="<<faint<<"\n";
		//for each line of sight
		for(int j = 0; j<i; j++){
			int N = Nstar[j];
			double l=glon[j], b=glat[j];
			std::string stringl = utils::toString(l), stringb = utils::toString(b);
			
			//convert l and b to radians
			obs::solarShifter shifter(intUnits);
			obs::PosSky pos(l*M_PI/180, b*M_PI/180);
			//	dust::dustModel dsty(7, .2, 1.6, shifter.xyz(), intUnits);//1.6 mag/kpc standard extinction in V
			dust::dustModel dsty(densityGasDisk, 1.6, &shifter, intUnits);//1.6 mag/kpc standard extinction in V
		//	obs::los los(pos, shifter.xyz()); //second argument should be a poscar of sun
			obs::los los(pos, shifter.xyz(), intUnits,&dsty); //second argument should be a poscar of sun //intUnits added for sKpc
			double s1=3.2;
			printf("A_H(%lf): %lf\n",s1,los.A_H(s1));
			printf("%f\n", los.A_V(0.001)/(0.001*intUnits.to_Kpc));
			
			std::cout<<"Sampling for line of sight l="+stringl+" " +"b=" + stringb + " with "<<N<<" samples\n"; //print to console
			//create a new directory here for each los
			std::string dir2(dir + "lines_of_sight/LOS("+stringl+","+stringb+")/");
			std::cout << "using " + dir2 + "\n";
			if(_mkdir(dir2.c_str())) printf("%s not created: already exists?\n",dir2.c_str());
			
			std::vector<double> dens = writeIntegral(dir2, stringl, stringb, &los, modelStars, bright, faint, &shifter);
			//sum for component fractions
			double sum = 0;
			for(int i=0;i<dens.size();i++) sum+=dens[i];
			//sample the model for each component
			std::vector<std::pair<std::vector<coord::PosVelCyl>, int>> results;
			std::vector<coord::PosVelCyl> posVels; 
			//edit to ensure sampling even for low star numbers
			std::random_device rd;
			std::default_random_engine eng(rd());
			std::uniform_real_distribution<float> distr(0, 1);
			for(int n =1; n<=N; n++){
				float r=distr(eng);
				if(r<=dens[0]/sum){
					posVels = galaxymodel::sampleLOS(modelComponents[0], &los, 1, bright, faint);
					results.push_back(std::make_pair(posVels,0));
				}
				else if(r<=(dens[0]+dens[1])/sum){
					posVels = galaxymodel::sampleLOS(modelComponents[1], &los, 1, bright, faint);
					results.push_back(std::make_pair(posVels,1));
				}
				else if(r<=(dens[0]+dens[1]+dens[2])/sum){
					posVels = galaxymodel::sampleLOS(modelComponents[2], &los, 1, bright, faint);
					results.push_back(std::make_pair(posVels,2));
				}
				else if(r<=(dens[0]+dens[1]+dens[2]+dens[3])/sum){
					posVels = galaxymodel::sampleLOS(modelComponents[3], &los, 1, bright, faint);
					results.push_back(std::make_pair(posVels,3));
				}
				else if(r<=(dens[0]+dens[1]+dens[2]+dens[3]+dens[4])/sum){
					posVels = galaxymodel::sampleLOS(modelComponents[4], &los, 1, bright, faint);
					results.push_back(std::make_pair(posVels,4));
				}
				else{
					posVels = galaxymodel::sampleLOS(modelComponents[5], &los, 1, bright, faint);
					results.push_back(std::make_pair(posVels,5));
				}
			}
/* 			for(int i=0;i<modelComponents.size();i++){
				posVels = galaxymodel::sampleLOS(modelComponents[i], &los, N*dens[i]/sum, bright, faint);
				results.push_back(std::make_pair(posVels,i));
			}	 */
			
			//write PosVelCyl to file
			writePosVel(dir2, stringl, stringb, results, bright, faint, &shifter);	
			//write actions, eccentricity and Zmax to a file
			writeJEZ(dir2, stringl, stringb, results, model, bright, faint);
			
			std::cout << "Results left in " + dir2 << '\n'; 
		}
	}
	else{
		//input l and b in degrees as command line arguments
		double L, B, l, b;
		sscanf(args[3],"%lf",&L); l=L;
		sscanf(args[4],"%lf",&B); b=B;
		int N = 500; //number of samples 	
		std::string stringl = utils::toString(l), stringb = utils::toString(b);
		//convert l and b to radians
		obs::solarShifter shifter(intUnits);
		obs::PosSky pos(l*M_PI/180, b*M_PI/180);
		//	dust::dustModel dsty(7, .2, 1.6, shifter.xyz(), intUnits);//1.6 mag/kpc standard extinction in V
		dust::dustModel dsty(densityGasDisk, 1.6, &shifter, intUnits);//1.6 mag/kpc standard extinction in V
	//	obs::los los(pos, shifter.xyz()); //second argument should be a poscar of sun
		obs::los los(pos, shifter.xyz(), intUnits,&dsty); //second argument should be a poscar of sun //intUnits added for sKpc
		double s1=3.2;
		printf("A_H(%lf): %lf\n",s1,los.A_H(s1));
		printf("%f\n", los.A_V(0.001)/(0.001*intUnits.to_Kpc));
		double bright=10, faint=15;
		std::cout<<"Sampling for line of sight l="+stringl+" " +"b=" + stringb + "\n"; //print to console
		
		//create a new directory here for each los
		std::string dir2(dir + "lines_of_sight/LOS("+args[3]+","+args[4]+")/");
		std::cout << "using " + dir2 + "\n";
		if(_mkdir(dir2.c_str())) printf("Did not mkdir %s\n",dir2.c_str());
		
		//write Integrals for each component to a file 
		std::vector<double> dens = writeIntegral(dir2, stringl, stringb, &los, modelStars, bright, faint, &shifter);
		std::cout<<"integral done\n";
		//sum for component fractions
		double sum = 0;
		for(int i=0;i<dens.size();i++) sum+=dens[i];
		//sample the model for each component
		std::vector<std::pair<std::vector<coord::PosVelCyl>, int>> results;
		std::vector<coord::PosVelCyl> posVels; 
		for(int i=0;i<modelComponents.size();i++){
			posVels = galaxymodel::sampleLOS(modelComponents[i], &los, N*dens[i]/sum, bright, faint);
			results.push_back(std::make_pair(posVels,i));
		}		
		
		//write PosVelCyl to file
		writePosVel(dir2, stringl, stringb, results, bright, faint, &shifter);	
		//write actions, eccentricity and Zmax to a file
		writeJEZ(dir2, stringl, stringb, results, model, bright, faint);		
		std::cout << "Results left in " + dir2 << '\n'; 
	}
	
///write potential to file so it does not have to be worked out from scratch each time
	if(newModel) writePotential(dir + "potential", *model.totalPotential, extUnits);
    duration -= std::time(NULL);
    std::cout << -duration << " secs for diagnostics\n";
	
	return 0;
}

int main(int narg, char **argv){
	try{
		prog(narg, argv);
	}
	catch(std::exception& e)
	{
		printf("%s", e.what());
		printf("hello");
	}
}
