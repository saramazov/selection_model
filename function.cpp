#include <cmath>
#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include "function.h"
#include <random>

using namespace std;




double rand2() {

random_device rd;unsigned long int seed=rd();
unsigned int m_a=1664525;
unsigned int m_c=1013904223;
unsigned int m_m=pow(2,31);

double d;
int n[2];

n[0]=seed;
n[1]=(m_a*n[0]+m_c)% m_m;
seed=n[1];
d=static_cast<double>(n[1])/static_cast<double>(m_m-1);


return d;
}


int uniform(int a,int b) {

double x=rand2();
double n=a+(b-a)*x;

n=n+0.5;
int m = (int)n;

return m;
}

int bernoulli(){

int pu=uniform(1,10);

if (pu<6){
	return 1;
}
else{
	return 0;
}

}


double* probability(double** x,int j,int n,double* sel,int n_sel,double* sigma_sel, double* prob, int choice2, double a, double s, int conf,int bins, double b, int power) {

//double *prob;
double tot=0;
int appo;
double appo2;
//double* histo;
//int* cont,where;

double* histo=new double[bins+1];
double* accu= new double[bins+1];
int* cont=new int[bins+1];
int* where=new int[n];
double* bias=new double[bins+1];

//double casual, appo;

//prob = new double[n+1];

for(int p=0; p<=bins; p++){

	histo[p]=0;
	cont[p]=0;
	bias[p]=0;
	accu[p]=0;

}

for(int p=0; p<n; p++){

	where[p]=0;

}

for(int p=0; p<=n; p++){

	prob[p]=0;

}

histo[bins]=1;

if((choice2==1)&&(conf==0)){

	for(int p=1; p<=n; p++){
		
		for(int u=0; u<n_sel; u++){

			appo2=exp((-pow((x[p-1][j-1]-sel[u]),2))/(2*pow(sigma_sel[u],2)));
			appo2=appo2/sigma_sel[u];
			prob[p]+=appo2;
				
		}

	tot+=prob[p];


	}
}


if((choice2==1)&&(conf==1)){


	for(int s=1; s<bins; s++){

		histo[s]=(static_cast<double>(s))/(static_cast<double>(bins));
	
	}

	for(int m=0; m<n; m++){

		for(int s=0; s<bins; s++){

			if((x[m][j-1]>=histo[s])&&(x[m][j-1]<histo[s+1])) {

				cont[s+1]++;
				where[m]=s+1;	
			
			}
			
			else if(x[m][j-1]==1){

				cont[bins]++;
				where[m]=bins;
			
			}

		}
	
	}

	for(int m=1; m<=bins; m++){

		bias[m]=pow((1+cont[m]*b),power);

	}


	for(int p=1; p<=n; p++){
			
		for(int u=0; u<n_sel; u++){

			appo=where[p-1];
			appo2=(bias[appo])*exp((-pow((x[p-1][j-1]-sel[u]),2))/(2*pow(sigma_sel[u],2)));
			appo2=appo2/sigma_sel[u];
			prob[p]+=appo2;			
				
		}

	tot+=prob[p];

	}
}


if((choice2==0)&&(conf==0)){


	for(int p=1; p<=n; p++){

		prob[p]=a+s*x[p-1][j-1];
		tot+=prob[p];

	}


}


if((choice2==0)&&(conf==1)){


	for(int r=1; r<bins; r++){

		histo[r]=(static_cast<double>(r))/(static_cast<double>(bins));
	
	}

	for(int m=0; m<n; m++){

		for(int r=0; r<bins; r++){

			if((x[m][j-1]>=histo[r])&&(x[m][j-1]<histo[r+1])) {

				cont[r+1]++;
				where[m]=r+1;	
			
			}
			
			else if(x[m][j-1]==1){

				cont[bins]++;
				where[m]=bins;
			
			}

		}
	
	}

	for(int m=1; m<=bins; m++){

		bias[m]=pow((1+cont[m]*b),power);

	}


	for(int p=1; p<=n; p++){
			
		for(int u=0; u<n_sel; u++){
			appo=where[p-1];
			prob[p]+=(bias[appo])*(a+s*x[p-1][j-1]);
				
		}

	tot+=prob[p];

	}
}



for(int p=1; p<=n; p++){

	prob[p]=prob[p]/tot;

}


if(conf==1){

	for(int q=1; q<=bins; q++){

		for(int p=1; p<=n; p++){

			if(where[p-1]==q) {

				accu[q]=accu[q]+prob[p];	

			}
		}

		//cout<<"Accu= "<<accu[q]<<endl;		
		
	}

//cout<<endl;

}


for(int p=1; p<=n; p++){

	prob[p]=prob[p]+prob[p-1];

}



delete[] histo;
delete[] cont;
delete[] where;
delete[] bias;
delete[] accu;

return prob;
}



double select(double* prob,int n){

double casual;
int appo=0;

casual=rand2();

//cout<<"casual= "<<casual<<endl;

for(int p=0; p<n; p++){

	if((casual>=prob[p]) && (casual<=(prob[p+1]))){
		appo=p;
		//cout<<"appo"<<endl;
		//cout<<"prob= "<<prob[p]<<endl;
	}	
}

//p è l'indice del modello da copiare

//cout<<"Appo= "<<appo<<endl;
return appo;
}



double blending(double* prob,int n,int models,double** x,double* w,int j){

double* casual;
double variant=0;
double accu=0;
int appo2=0;
int repeat=0;

casual= new double[models];
int *appo= new int[models];

for(int h=0; h<models; h++){

	appo[h]=0;

}

//cout<<"casual= "<<casual<<endl;

for(int h=0; h<models; h++){

	do{

		repeat=0;
		casual[h]=rand2();	

		for(int p=0; p<n; p++){

			if((casual[h]>=prob[p]) && (casual[h]<=(prob[p+1]))){
				appo[h]=p;
				//cout<<"prob= "<<prob[p]<<endl;
			}	
		}

		for(int k=0; k<h; k++){

			if(appo[k]==appo[h]){

				repeat=1;
			
			}	
		}	

	} while(repeat==1);

}
//p è l'indice del modello da copiare


for(int h=0; h<models; h++){

	appo2=appo[h];
	variant+=(x[appo2][j-1]*w[appo2]);
	accu+=w[appo2];
	//cout<<appo[h]<<" ";

}

//cout<<endl;

variant=variant/accu;
//cout<<"Variant= "<<variant<<endl;
delete[] casual;
delete[] appo;

return variant;
}



double attract(double** x,int i,int j,double* att,int n_att,double* sigma_att,int choice,int det,double* threshold,double* charge,double* elastic,double limit,double* range,int* strength,double* coefficient) {

double mean;
double coeff;
double a=0;
double sum=0;
int exponent;
double r;

/*if((det==1) && (choice==1)){

	for(int k=0; k<n_att; k++){
		

		if(x[i][j]>att[k]){

			if((x[i][j]-att[k])<limit){

				sum=-(x[i][j]-att[k]);
				break;
			
			}
			
			a = -charge[k]/(pow((x[i][j]-att[k]),2));

			//cout<<"a0= "<<a<<endl;

			if((x[i][j]-att[k]+a) >= 0){

				sum+=a;
			//	cout<<"sum0= "<<sum<<endl;
				
			}			
			
			else{

				sum+= -(x[i][j]-att[k]);
			//	cout<<"sum0= "<<sum<<endl;

			}

		}

		else {

			if((att[k]-x[i][j])<limit){

				sum=att[k]-x[i][j];
				break;
			
			}
			
			a = charge[k]/(pow((x[i][j]-att[k]),2));
		
			//cout<<"a1= "<<a<<endl;
			
			if((att[k]-x[i][j]-a) >= 0){

				sum+=a;
			//	cout<<"sum1= "<<sum<<endl;
				
			}			
			
			else{

				sum+= att[k]-x[i][j];
			//	cout<<"sum1= "<<sum<<endl;

			}

		}	
		
	
	}
}*/


if((det==1) && (choice==1)){

	for(int k=0; k<n_att; k++){
		

		if(x[i][j]>att[k]){

			if((x[i][j]-att[k])<limit){

				sum=-(x[i][j]-att[k]);
				break;
			
			}
			
			exponent=strength[k];				
			r=range[k];
			a=-(pow(r,exponent))/(pow((x[i][j]-att[k]),exponent));
			a=charge[k]*a;

			//cout<<"a0= "<<a<<endl;

			if((x[i][j]-att[k]+a) >= 0){

				sum+=a;
			//	cout<<"sum0= "<<sum<<endl;
				
			}			
			
			else{

				sum+= -(x[i][j]-att[k]);
			//	cout<<"sum0= "<<sum<<endl;

			}

		}

		else {

			if((att[k]-x[i][j])<limit){

				sum=att[k]-x[i][j];
				break;
			
			}
			
			exponent=strength[k];				
			r=range[k];
			a=(pow(r,exponent))/(pow((att[k]-x[i][j]),exponent));
			a=charge[k]*a;
		
			//cout<<"a1= "<<a<<endl;
			
			if((att[k]-x[i][j]-a) >= 0){

				sum+=a;
			//	cout<<"sum1= "<<sum<<endl;
				
			}			
			
			else{

				sum+= att[k]-x[i][j];
			//	cout<<"sum1= "<<sum<<endl;

			}

		}	
		
	
	}
}



/*if((det==0) && (choice==1)){

	for(int k=0; k<n_att; k++){
		
		if(x[i][j]>att[k]){
			
			if((x[i][j]-att[k])<limit){

				sum=-(x[i][j]-att[k]);
				break;
			
			}
						
			mean = -charge[k]/(pow((x[i][j]-att[k]),2));
			a = gauss(mean,sigma_att[k]);

			if((x[i][j]-att[k]+a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= -(x[i][j]-att[k]);

			}

		}

		else {

			if((att[k]-x[i][j])<limit){

				sum=att[k]-x[i][j];
				break;
			
			}
			
			mean = charge[k]/(pow((x[i][j]-att[k]),2));
			a = gauss(mean,sigma_att[k]);

			if((att[k]-x[i][j]-a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= att[k]-x[i][j];

			}

		}	
		
	
	}
}*/



if((det==0) && (choice==1)){

	for(int k=0; k<n_att; k++){
		
		if(x[i][j]>att[k]){
			
			if((x[i][j]-att[k])<limit){

				sum=-(x[i][j]-att[k]);
				break;
			
			}
			
			exponent=strength[k];
			r=range[k];			
			mean = -(pow(r,exponent))/(pow((x[i][j]-att[k]),exponent));
			//cout<<"mean= "<<mean<<endl;
			a = gauss(mean,sigma_att[k]);
			//cout<<"a1= "<<a<<endl;
			a=a*charge[k];
			//cout<<"a2= "<<a<<endl;

			if((x[i][j]-att[k]+a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= -(x[i][j]-att[k]);

			}

		}

		else {

			if((att[k]-x[i][j])<limit){

				sum=att[k]-x[i][j];
				break;
			
			}
			
			exponent=strength[k];
			r=range[k];			
			mean = (pow(r,exponent))/(pow((att[k]-x[i][j]),exponent));
			//cout<<"Mean= "<<mean<<endl;
			a = gauss(mean,sigma_att[k]);
			//cout<<"a1= "<<a<<endl;
			a=a*charge[k];
			//cout<<"a2= "<<a<<endl;

			if((att[k]-x[i][j]-a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= att[k]-x[i][j];

			}

		}	
		
	
	}
}



if((det==1) && (choice==0)){

	for(int k=0; k<n_att; k++){
		
		if((x[i][j]>att[k]) && (x[i][j]<=threshold[k+1])){

			if((x[i][j]-att[k])<limit){

				sum=-(x[i][j]-att[k]);
				break;
			
			}
			
			a = -elastic[k]*(x[i][j]-att[k]);
		
			if((x[i][j]-att[k]+a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= -(x[i][j]-att[k]);

			}

		}

		else if((x[i][j]<=att[k]) && (x[i][j]>=threshold[k])) {

			if((att[k]-x[i][j])<limit){

				sum=att[k]-x[i][j];
				break;
			
			}

			a = elastic[k]*(att[k]-x[i][j]);
			
			if((att[k]-x[i][j]-a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= att[k]-x[i][j];

			}

		}	
		
	
	}
}





/*if((det==0) && (choice==0)){

	for(int k=0; k<n_att; k++){
		
		if((x[i][j]>=att[k]) && (x[i][j]<=threshold[k+1])){

			if((x[i][j]-att[k])<limit){

				sum=-(x[i][j]-att[k]);
				break;
			
			}
			
			mean = -elastic[k]*(x[i][j]-att[k]);
			a = gauss(mean,sigma_att[k]);

			if((x[i][j]-att[k]+a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= -(x[i][j]-att[k]);

			}

		}

		else if ((x[i][j]<att[k]) && (x[i][j]>=threshold[k])) {

			if((att[k]-x[i][j])<limit){

				sum=att[k]-x[i][j];
				break;
			
			}

			mean = elastic[k]*(att[k]-x[i][j]);
			a = gauss(mean,sigma_att[k]);

			if((att[k]-x[i][j]-a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= att[k]-x[i][j];

			}

		}	
		
	
	}
}*/


if((det==0) && (choice==0)){

	for(int k=0; k<n_att; k++){
		
		if((x[i][j]>=att[k]) && (x[i][j]<=threshold[k+1])){

			if((x[i][j]-att[k])<limit){

				sum=-(x[i][j]-att[k]);
				break;
			
			}
			
			mean = -(x[i][j]-att[k]);
			a = gauss(mean,sigma_att[k]);
			a=a*elastic[k];

			if((x[i][j]-att[k]+a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= -(x[i][j]-att[k]);

			}

		}

		else if ((x[i][j]<att[k]) && (x[i][j]>=threshold[k])) {

			if((att[k]-x[i][j])<limit){

				sum=att[k]-x[i][j];
				break;
			
			}

			mean = (att[k]-x[i][j]);
			a = gauss(mean,sigma_att[k]);
			a=a*elastic[k];

			if((att[k]-x[i][j]-a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= att[k]-x[i][j];

			}

		}	
		
	
	}
}


if((det==1) && (choice==2)){

	for(int k=0; k<n_att; k++){
		

		if(x[i][j]>att[k]){

			if((x[i][j]-att[k])<limit){

				sum=-(x[i][j]-att[k]);
				break;
			
			}
			
			coeff=coefficient[k];				
			r=range[k];
			a=-coeff*exp(-(x[i][j]-att[k])/r);
			a=charge[k]*a;

			//cout<<"a0= "<<a<<endl;

			if((x[i][j]-att[k]+a) >= 0){

				sum+=a;
			//	cout<<"sum0= "<<sum<<endl;
				
			}			
			
			else{

				sum+= -(x[i][j]-att[k]);
			//	cout<<"sum0= "<<sum<<endl;

			}

		}

		else {

			if((att[k]-x[i][j])<limit){

				sum=att[k]-x[i][j];
				break;
			
			}
			
			coeff=coefficient[k];				
			r=range[k];
			a=coeff*exp(-(att[k]-x[i][j])/r);
			a=charge[k]*a;

		
			//cout<<"a1= "<<a<<endl;
			
			if((att[k]-x[i][j]-a) >= 0){

				sum+=a;
			//	cout<<"sum1= "<<sum<<endl;
				
			}			
			
			else{

				sum+= att[k]-x[i][j];
			//	cout<<"sum1= "<<sum<<endl;

			}

		}	
		
	
	}
}


if((det==0) && (choice==2)){

	for(int k=0; k<n_att; k++){
		
		if(x[i][j]>att[k]){
			
			if((x[i][j]-att[k])<limit){

				sum=-(x[i][j]-att[k]);
				break;
			
			}
			
			coeff=coefficient[k];				
			r=range[k];
			mean=-coeff*exp(-(x[i][j]-att[k])/r);
			//cout<<"mean= "<<mean<<endl;
			a = gauss(mean,sigma_att[k]);
			//cout<<"a1= "<<a<<endl;
			a=a*charge[k];
			//cout<<"a2= "<<a<<endl;

			if((x[i][j]-att[k]+a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= -(x[i][j]-att[k]);

			}

		}

		else {

			if((att[k]-x[i][j])<limit){

				sum=att[k]-x[i][j];
				break;
			
			}
			
			coeff=coefficient[k];				
			r=range[k];
			mean=coeff*exp(-(att[k]-x[i][j])/r);
			//cout<<"mean= "<<mean<<endl;
			a = gauss(mean,sigma_att[k]);
			//cout<<"a1= "<<a<<endl;
			a=a*charge[k];
			//cout<<"a2= "<<a<<endl;

			if((att[k]-x[i][j]-a) >= 0){

				sum+=a;
				
			}			
			
			else{

				sum+= att[k]-x[i][j];

			}

		}	
		
	
	}
}




//cout<<"sum= "<<sum<<endl;

return sum;

}



double gauss(double mean,double sigma){

double a=mean-5.0*sigma;
double b=mean+5.0*sigma;

double x,y,f;

//MonteCarlo method

	do {
		
		x=a+(b-a)*rand2();
		
		y=rand2()/sigma;
		f=exp((-pow((x-mean),2))/(2*pow(sigma,2)));
		f=f/sigma;
		
	}while (f < y); 

return x;
}
