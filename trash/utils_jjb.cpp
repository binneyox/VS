#define SIGN(A,B) (std::signbit(B)?(-A):(A))
EXP std::string jjbcmp(double x){
	const int i1 = 90, i2 = 8100, i3 = 729000;
	const double at2t23 = 8388608., big=1.e28, small=1.e-28, aln2=log(2.);
	char m[6]={80,0,0,0,0,'\0'};
	int lnx,rem; double absx=fabs(x);
	if (x!=0){
		if (absx>big) x= SIGN(big, x);
		if (absx<small) x = SIGN(small, x);
		if (absx>=1.)
			lnx = (int)(log(absx)/aln2) + 129;
		else
			lnx = (int)(log(absx)/aln2) + 128;
		int j = lnx/i1;
		m[3] = lnx - j*i1;
		lnx = lnx - 129;
		if (lnx>0)
			rem = (int)((absx/pow(2.,(double)(lnx)) - 1.)*at2t23);
		else
			rem = (int)((absx*pow(2.,(double)(fabs(lnx))) - 1.)*at2t23);
		m[4] = rem/i3;
		rem = rem - m[4]*i3;
		m[2] = rem/i2;
		rem = rem - m[2]*i2;
		m[1] = rem/i1;
		m[0] = rem - m[1]*i1;
		m[4] = m[4] + j*12;
		if (x<0.) m[4] = m[4] + 40;
	}
	for (int i=0;i<5;i++){
		m[i] = m[i] + 33;
		if (m[i]>=94) m[i] = m[i] + 1;
		if (m[i]>=123) m[i] = m[i] + 1;
	}
	std::string st(m);
	return st;
}

EXP bool jjbget(FILE *tmpf,float *xg,int npt){
	float expon, at2t23= 8388608.;
	long int i1=90;
	unsigned char m[5];
	bool ok=true;
	for(int k=0; k<npt; k++){
		m[0]=getc(tmpf); if(m[0]=='\n') m[0]=getc(tmpf);
		if(feof(tmpf)){
			ok=false; break;
		}

		for(int j=1;j<5;j++) m[j]=getc(tmpf);
		for(int i = 0; i<5; i++){
			if (m[i]>=124) m[i] = m[i] - 1;
			if (m[i]>=95) m[i] = m[i] - 1;
			m[i] = m[i] - 33;
		}
		if (m[4]==80) xg[k] = 0;
		else{
			if (m[4]>=40){
				m[4] = m[4] - 40;
				xg[k] = -1.;
			}else xg[k] = 1;
			int j = m[4]/12;
			m[4] = m[4] - j*12;
			expon=(float)(j*i1+m[3]-129);
			xg[k]=xg[k]*(float)((((m[4]*i1+m[2])*i1+m[1])*i1
				      +m[0])/at2t23+ 1.);
			if (expon>0.)
				xg[k]=xg[k]*pow(2.,expon);
			else
				xg[k]=xg[k]/pow(2.,-expon);
		}
	}
	return ok;
}
