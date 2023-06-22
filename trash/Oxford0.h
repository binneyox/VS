


/// Now the Pascale, Binney, Nipoti (2020 DF
double Oxford_rat(const double[],double,double,double);
double Bologna_rat(const double[],double,double);

class OB_syst: public math::IOdeSystem{
	private:
		int OB; double coefL, k, gamma0, gamma1;
	public:
		OB_syst(int _OB,double _coefL,double _k,double _gamma0,double _gamma1) :
		    OB(_OB), coefL(_coefL), k(_k), gamma0(_gamma0), gamma1(_gamma1){};
		void eval(const double,const double[],double[]) const;
		unsigned int size()  const {return 3;}
};

double get_h(int,double[],double,double,double,double);

struct OxfordParam{
	double
			mass,      ///< total mass in Msun
			J0,        ///< break action (defines the transition between inner and outer regions)
			Jcore,        ///< controls Cole&Binney core size
			Jcutoff,   ///< cutoff action (sets exponential suppression at J>Jcutoff, 0 to disable)
			slopeIn,   ///< power-law index for actions below the break action (Gamma)
			slopeOut,  ///< power-law index for actions above the break action (Beta)
			steepness, ///< steepness of the transition between two asymptotic regimes (eta)
			cutoffStrength, ///< steepness of exponential suppression at J>Jcutoff (zeta)
			coefLin,  ///< radial anisotropy
			kIn,      ///< flattening
			coefLout,  ///< radial anisotropy
			kOut,      ///< flattening
			L0,        ///< epicycle freqs vs r
			gamma0,     ///< epicycle freqs as r->0
			gamma1,     ///< epicycle freqs as r->0
			rotFrac,   ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
			Jphi0;     ///< controls the steepness of rotation and the size of non-rotating core

	OxfordParam() :  ///< set default values for all fields (NAN means that it must be set manually)
	    mass(NAN), J0(NAN), Jcore(0), Jcutoff(0), slopeIn(NAN), slopeOut(NAN), steepness(1), cutoffStrength(2),
	    coefLin(1), kIn(1), coefLout(1), kOut(1), L0(NAN), gamma0(0), gamma1(1), rotFrac(0), Jphi0(0) {}
};
class Oxford: public BaseDistributionFunction{
	const OxfordParam par;  ///< parameters of DF
	private:
		double norm,beta;
		double g(const double hJ) const;
		void set_beta(void);
		Oxford(const OxfordParam &params_,double beta_) : par(params_), beta(beta_){}
	public:
    /** Create an instance of ColeBinney double-power-law distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
		Oxford(const OxfordParam&);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
		virtual double value(const actions::Actions &J) const;
		virtual void set_norm(double);
		virtual void write_params(std::ofstream&,const units::InternalUnits&) const;
		virtual void tab_params(std::ofstream&,const units::InternalUnits&) const;
		double diff(const double) const;
		double intDiff(const double) const;
};

//We start with two helper classes used to determine Cole&Binney
//parameter beta

class OxfordDiff: public math::IFunctionNoDeriv{
	public:
		const Oxford& Sft;
		OxfordDiff(Oxford &Sft_): Sft(Sft_) {}
		virtual double value(const double hJ) const{
			return Sft.diff(hJ);
		}
};
class OxfordInt: public math::IFunctionNoDeriv{
	public:
		const Oxford& Sft;
		OxfordInt(Oxford &Sft_): Sft(Sft_) {}
		virtual double value(const double x) const{
			return Sft.intDiff(x);
		}
};

struct BolognaParam {
	double
			mass,      ///< total mass of a model
			J0,        ///< controls core size
			Jphi,      ///< controls the steepness of rotation and the size of non-rotating core
			L0,        ///< controls approx to Omega_phi/Omega_r
			gamma,     ///< epicycle ratio r->0
			L1,        ///< introduces harmonic core
			coefL,      ///< control radial anisotropy
			kIn, kOut,  ///< control flattening
			alpha,      ///< controls out density profile 
			rotFrac;    ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
	BolognaParam():     ///< set default values for all fields (NAN means that it must be set manually)
	    mass(NAN), J0(NAN), Jphi(0), L0(NAN), gamma(0), L1(0), coefL(0), kIn(1), kOut(1),
	    alpha(0.5), rotFrac(0) {}
};

class Bologna: public BaseDistributionFunction{
	const BolognaParam par;  ///< parameters of DF
	private:
		double norm;
	public:
    /** Create an instance of Pascale Exponential distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
		Bologna(const BolognaParam &params);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
		virtual double value(const actions::Actions &J) const;
};

