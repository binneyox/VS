
/*
   Based on program by D. Wood. modified by D. Spergel sept 82 & J Binney
   august 83 this defines a class to identify sinusoidal terms in the DFT of a time
   series. ts(NF/2) is the sine transform and tc(NF/2+1) is the cosine
   transform of a data stream NF long with data points separated by Dt so
   last point is (NF-1)Dt later than the first. Only terms with
   frequencies between frqmin (default 0) and frqmax (default 1e6) are added
   to the returned list of lines. The resid member of a line is
   the fraction of the original power left in the spectrum after that line's
   contributioon has been subtracted, so successful functioning of
   analyse() is very small resid for the last extracted (weakest) line.
   */
/*
 A class to hold information about a sinusoidal contribution to the
 spectrum
*/
class EXP line{
	private:
	public:
		line(void){};
		~line(void){};
		line& operator = (const line&);
		double A;
		double nu;
		double phi;
		double resid;
		unsigned int diag;//0 reg isol, 1 zero freq, 2 low freq, 3 pair
};
/*
   The class of spectrum analysers. The actual work is done by analyse()
*/
class EXP FrequencyFinder{
	private:
		bool delete_tc;
		int NF,NF5;
		double* ts;//input array of length NF/2
		double* tc;//input array of length NF/2+1 
		double* z;// work space seized and released 
		const double deltat;// interval between input times
		const double frqmin, frqmax;//lowest & highest frequencies to be saved
		double zs[3],zc[3];//2nd difference of spectrum near current peak
		std::vector<line> lines;//object returned by analyse()
		double wmin,/*freq resolution*/ discrm;//angle difference considered insignificant
		double pwr0,/*input power*/ sin0, cos0;//sin and cos pi/NF 
		double get_pwr();
		int difference(double&);//comput 2nd difference and return location and value of peak
		void subtract(double,double,double,int);//subtract a line's contribution to spectrum
		double isolated(int,std::vector<line>&);//analyse an isolated peak
		double lowFreq(int,std::vector<line>&);//analyse a peak near zero frequency
		double pair(int,std::vector<line>&);//analyse a peak possibly caused by two frequencies
		void order(std::vector<line>&);//amalgamate lines (not used)
	public:
		//Initialise with sine ts & cosine ts transforms
		FrequencyFinder(double* _ts,double* _tc,const int _NF,const double _dt,
				const double _frqmin=0,const double _frqmax=1e6);
		//Initialise with time series data
		FrequencyFinder(double* data,const int _NF,const double _dt,
				const double _frqmin=0,const double _frqmax=1e6);
		~FrequencyFinder(){
			delete[] z;
			if(delete_tc){
				delete[] tc; delete[] ts;
			}
		}
		std::vector<line> analyse();
};

