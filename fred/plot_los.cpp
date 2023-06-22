#include <stdio.h>
#include <math.h>
#include "/u/sm/mgo.h"

void plot_mags(char **argv, mgo::plt &pl){
	char fname[40];
	const int nb=40;
	float bins[nb],Ns[nb],faint=16,bright=8, db=(faint-bright)/(float)nb;
	bins[0]=bright; Ns[0]=0;
	for(int i=1;i<nb;i++){
		bins[i]=bins[i-1]+db; Ns[i]=0;
	}
	sprintf(fname,"lines_of_sight/los_%s_%s/mags",argv[1],argv[2]);
	printf("opening %s|\n",fname);
	FILE *ifile;
	if(fopen_s(&ifile,fname,"r")){
		printf("can't open %s\n",fname); return;
	}
	float mag;
	int ns=0;
	while(!feof(ifile)){
		if(1!=fscanf(ifile,"%f",&mag)) break;
		int i=0;
		while(mag>bins[i] && i<=nb) i++;
		if(i<nb){
			Ns[i]+=1; ns++;
		}
	}
	float ds=.5*(bins[1]-bins[0]),Nmax=ns/10;
	for(int i=0;i<nb;i++)bins[i]-=ds;
	pl.new_plot(bins[0],bins[nb-1],0,Nmax,"mag","N");
	pl.histogram(bins,Ns,nb);
	pl.grend();
}
		

int main(int narg,char** argv){
	char fname[40];
	sprintf(fname,"lines_of_sight/los_%s_%s/xvCar",argv[1],argv[2]);
	printf("opening %s|\n",fname);
	FILE *ifile;
	if(fopen_s(&ifile,fname,"r")){
		printf("can't open %s\n",fname); return 0;
	}
	float ell,b,torad=acos(-1)/180;
	sscanf(argv[1],"%f",&ell); sscanf(argv[2],"%f",&b);
	ell*=torad; b*=torad;
	const int nb=50;
	float bins[nb],db=20/(float)nb,Ns[nb];
	for(int i=0;i<nb;i++){
		bins[i]=(1+i)*db;Ns[i]=0;
	}
	float bright,faint;
	fscanf(ifile,"%f %f",&bright,&faint);
	double Xs=-8.27,Zs=0.025;
	float s0=-(Xs*cos(ell)*cos(b)+Zs*sin(b));
	int ns=0;
	while(!feof(ifile)){
		double x,y,z,vx,vy,vz;
		if(6!=fscanf(ifile,"%lf %lf %lf %lf %lf %lf",&x,&y,&z,&vx,&vy,&vz)) break;
		double s=sqrt(pow(x-Xs,2)+pow(y,2)+pow(z-Zs,2));
		int i=0;
		while(s>bins[i] && i<=nb) i++;
		if(i<nb){
			Ns[i]+=1; ns++;
		}
	}
	printf("%d stars binned %f\n",ns,Ns[0]);
	float ds=.5*(bins[1]-bins[0]),Nmax=ns/10;
	for(int i=0;i<nb;i++)bins[i]-=ds;
	mgo::plt pl;
	pl.new_plot(0,bins[nb-1],0,Nmax,"s/kpc","N");
	pl.histogram(bins,Ns,nb);
	sprintf(fname,"b:%2.0f f:%2.0f",bright,faint);
	pl.relocate(.8*bins[nb-1],.9*Nmax); pl.label(fname);
	sprintf(fname,"(%s,%s)",argv[1],argv[2]);
	pl.relocate(.8*bins[nb-1],.8*Nmax);
	pl.label(fname);
	pl.setltype(2); pl.relocate(s0,0); pl.draw(s0,Nmax);
	pl.grend();
	plot_mags(argv,pl);
}
		
	
	
	