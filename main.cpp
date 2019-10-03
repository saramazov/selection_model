#include <cmath>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <fstream>
#include "function.h"


//Conformism. Random selection predicts that, if an individual is closer to a selection peak, it will be more likely to be chosen as model. Adding conformism means that this probability is incremented, according to how mnay individuals possess that variant. Models with a common variant will therefore be disproportionaly likely to be copied, i.e. more likely than the selective gaussian function would predict. Example: if 5 models have 0.9 and 3 models 0.3, of course 0.9 is favoured because there are more models there. But this ratio is furtherly incremented by a bias coefficient. 
//when i use elastic attractor-selection with threshold 0.8 and conformism and selection=0.05, and elastic constrnt 0.5, then 0 prevails. when I use 0.6 threshold but conformism stronger than selecton (0.5 vs 0.05) then sometimes 1 prevails, sometimes 0. 
//Blending imitation: each model is chosen following the gaussian function (so different models have different prob of being selected). then, however, each model is chosen indipdendelty -i.e-, no assortative mating. Finally, the intrinsic weight of the model is indiepndent of its probability of being chosen as model. 
//with blending and everyhing linear, threshold 0.6 and paramaeters as in BR: it happens soething similar to conformism. sometimes 1 prevails(faster), sometimes 0 prevails(very fast). I used 5 models. Same with 2. actually i am not sure it even reaches 1 sometimes! maybe always 0. should check. looks like it always goes to 0| if I move the threshold to 0.5, it goes to 1.
//stenght of attractor has different meanings: it might be the magnitude of delta x, which is represented by charge. Or it might be the differential probability of being attracted with respect to being repelled; more specifically, it is the area of the gaussian prob function for dx positive with respect to the are with dx negative. How can this are be affected? Fixed x, if two attractor have the same mean, then the bigger is sigma, the bigger is left area, i.e. the weaker the attractor. This is true for every x, irrespective of how far you are from the attractor. Another way is to change the mean. this can be achieved in two ways: either the range is bigger (0.1/x) instead of (0.01/x), or the exponent is different (1\x)^2 instead of 1/x. However, changing the range or sigma means that one attractor is always stronger than another one. What if you want an attractor to be stronger when you are close, but weaker when you are far, i.e. waht if the strength of an attractor is inversely proportional to its range? Then you must change the exponent, and the range tells you at what point one attractor becomes stronger/weaker than the other. 
//BR claim that selection predicts the outcome, since, even in sprber model, it is the position of selection peak with respect to threshold that determines whether one attractor is favoured over another one. However, this is a idyosincratic property of their model, because they used attractor with threshold. If you use gravitational attractor, then something different happens. If 0 and 1 are attractors, selection peak is at 0.8, and 0 is much stronger than 1 (0.5 vs 0.05), and selection not too strong (range=0.1), then the final outcome is very close to 0 (albeit not exactly 0, but 0.05).Important:use deterministic model. with charge as specified, range=1(inessential), and exponent 2.Increasing the range of 1 attractor (2 instead of 1), with charge 0.1 and 0.01, makes 1 prevail. problem with all this models: mean is incredibly big!!!
//revised with exponential attractor and probabilistic (peak selection at 0.7). if charges are small (both 0.001), range is 0.1, and selection has range 0.1, then the final outcome is 1, but it takes almost 500 time-steps, which means that selection has a role (albeit marginal) for final outcome. If 0 increases its charge to 0.01 and its range to 0.5, then 0 prevails! Now I try still with 0 stronger than 1, but weak selection (range=1): many oscillations, but in the end it reaches 1 (after quite a long time) . Now I try with peak at 0.7 (and strong sel, i.e. range=0.1), and same chrge (0.001), but range of 0 bigger than 1 (0.5 vs 0.1); why? because in this way the two delta-x should be equal exactly at x=0.7, where selection peaks. wat happens? it takes a lot of time (500time steps), but 1 is reached. so it looks like selection peak is never stable. Now with CONFORMISM:I try again with small charges and strong selection that peaks at 0.7, conformism=0.1: 1 is still reached. now I increase the range of 0, and conforist=0.1:0 prevails as before. Now 0 is stronger but selection weak, conformist still 0.1: it loks variable; I mean, sometimes it reaches 1,sometimes 0, so it's weird. finally, I try strong selection and same charge but different ranges for 0 and 1: a lot of oscillation, selection lokks quite stable but in the end 1 prevails, even if after 500 time-steps still the average is 0.9. Now I try conformism, attractors with different range (0.5 vs 0.1), same charge=0.001, conformism=0.1, selecton strong=0.1,  and I increase variance of 0 to sigma=1 (the other one is 0.1): incredible result! Selection peak look stable! outcome is 0.67 even at 500 time-steps. So maybe increasing the variance (that is, making it probabilistic), increases the chances that selection is reached. Well, not surprising: when it is deterministic, only attractors are reached, and which attractors is eached depends on relative strenth, i.e. where the two mean delta-x are equal with respect to where the selection peak is (so, very similar to threshold of BR). But probabilistic might change it! Again: even trying without conformism in the last simulation, selection peak is reached!Even if the range is reduced to 0.4, still selection peak wins (with no conformism!).if you reduce range to 0.2, then 1 is almst reacvhed (0.9 after 500 timesteps); but if increase the variance then 0.76 is reached (variance =1 for both); now i try to increase both variances to 2: 
//note: TIME to reach equilibrium is important, because attractor might change through time (perhaps psychological preferences are culturall-depdendent), therefore selection is imortant in the short-run. 
//note. Assumptions of sperber. they changed the model of br, but: they didn't apply it to the actual model; they didn't consider conformism or blending; they didn't try other configurations or different initial conditions. About the cigarette-model; they still assume a threshold; they fixed the mean and (sort of) variance of their distribution, and didn't try to change those parameters; the attractors peak not exactly at 0 and 30, but in some point closer-unrealistic assumption; the didn't try with conformism and blending; results are very confused and unconvincing. 

using namespace std;


int main() {


int n,t,choice,choice2,conf,n_att,n_sel,det,bins,power,blend,models;
int index=0;
double* att;
double* sel;
double* sigma_sel;
double* prob;
double* charge;
double* threshold;
double* elastic;
double* sigma_att;
double* average;
double* weight;
double* range;
double* coefficient;
double** x;
double delta=0;
double a,s,bias;
double limit;

cout<<"Number of individuals?"<<endl;
cin>>n;

cout<<"timesteps?"<<endl;
cin>>t;

cout<<"How many attractors?"<<endl;
cin>>n_att;

att=new double[n_att];
sigma_att=new double[n_att];
charge=new double[n_att];
threshold=new double[n_att+1];
elastic=new double[n_att];
range=new double[n_att];
coefficient=new double[n_att];
int* strength=new int[n_att];
average=new double[t];

for(int i=0; i<t; i++){

	average[i]=0;

}

for(int i=0;i<n_att;i++) {

	cout<<"Position of attractor?"<<endl;
	cin>>att[i];

	
}

cout<<"How many selection peaks?"<<endl;
cin>>n_sel;

sel=new double[n_sel];
sigma_sel=new double[n_sel];


for(int i=0;i<n_sel;i++) {

	cout<<"Position of selection peak?"<<endl;
	cin>>sel[i];


}


cout<<"Elastic or gravitational attractor? 0 for elastic, 1 for gravitational, 2 for exponential"<<endl;
cin>>choice;

if (choice==0){

	for(int i=0;i<n_att;i++) {

	cout<<"Elastic constant?"<<endl;
	cin>>elastic[i];

	cout<<"Variance of attractor?"<<endl;
	cin>>sigma_att[i];

	}

	threshold[0]=0;
	threshold[n_att]=1;

	for(int i=1;i<n_att;i++) {

	cout<<"Threshold?"<<endl;
	cin>>threshold[i];
	
	}

cout<<"Limit?"<<endl;
cin>>limit;
}

if (choice==1){
	
	for(int i=0;i<n_att;i++) {

	cout<<"Charge?"<<endl;
	cin>>charge[i];

	cout<<"Variance of attractor?"<<endl;
	cin>>sigma_att[i];

	cout<<"Range of attractor?"<<endl;
	cin>>range[i];

	cout<<"Strength of attractor? (exponent of function)"<<endl;
	cin>>strength[i];

	}

cout<<"Limit?"<<endl;
cin>>limit;	

}

if (choice==2){
	
	for(int i=0;i<n_att;i++) {

	cout<<"Charge?"<<endl;
	cin>>charge[i];

	cout<<"Variance of attractor?"<<endl;
	cin>>sigma_att[i];

	cout<<"Range of attractor?"<<endl;
	cin>>range[i];

	cout<<"Coefficient of attractor? (Suggestion is to use 1/range"<<endl;
	cin>>coefficient[i];

	}

cout<<"Limit?"<<endl;
cin>>limit;	

}



cout<<"Linear or gaussian selection? 0 for linear, 1 for gaussian"<<endl;
cin>>choice2;

cout<<"Conformism? 1 for yes, 0 for no"<<endl;
cin>>conf;

if(conf==1){

	cout<<"How many bins?"<<endl;
	cin>>bins;
	cout<<"Conformist bias?"<<endl;
	cin>>bias;
	cout<<"Exponent of bias function?"<<endl;
	cin>>power;

}

if (choice2==0){

	cout<<"Baseline fitness?"<<endl;
	cin>>a;

	cout<<"Selective advantage?"<<endl;
	cin>>s;
	
	}



if (choice2==1){
	
	for(int i=0;i<n_sel;i++) {

	cout<<"Range of selection peak?"<<endl;
	cin>>sigma_sel[i];

	}
	

}

cout<<"Blending transmission? 1 for yes, 0 for no"<<endl;
cin>>blend;

if(blend==1){

	cout<<"Number of models?"<<endl;
	cin>>models;

}

cout<<"Deterministic or not? 1 for yes, 0 for no"<<endl;
cin>>det;

x = new double*[n];
	for(int i = 0; i < n; i++){
		x[i] = new double[t]; 
		for(int j = 0; j < t; j++){		
			x[i][j] = 0;
		}
	}

for(int i = 0; i < n; i++){

	//x[i][0]=rand2();
	x[i][0]=(static_cast<double>(i))/(static_cast<double>(n));
	//x[i][0]=0.8;	
	//cout<<x[i][0]<<endl;	
	average[0]=average[0]+x[i][0];

}

weight=new double[n];

for(int i=0; i<n; i++){

	weight[i]=rand2();
	//cout<<weight[i]<<endl;

}

average[0]=average[0]/(static_cast<double>(n));

prob = new double[n+1];

for(int i = 0; i <= n; i++){

	prob[i]=0;

}




for (int j=1; j<t; j++){

	prob = probability(x,j,n,sel,n_sel,sigma_sel,prob,choice2,a,s,conf,bins,bias,power);

	/*for (int m=1;m<=n;m++) {
		
		cout<<x[m-1][j-1]<<" ";
		cout<<prob[m]<<" ";
		cout<<endl;

	}*/
	

	for (int i=0; i<n; i++){

		if(blend==0){

			index=select(prob,n);
			x[i][j]=x[index][j-1];
			average[j]=average[j]+x[i][j];
			//cout<<x[i][j]<<endl;
		}

		else if(blend==1){

			x[i][j]=blending(prob,n,models,x,weight,j);
			average[j]=average[j]+x[i][j];
			//cout<<x[i][j]<<endl;
		
		}		

		delta=attract(x,i,j,att,n_att,sigma_att,choice,det,threshold,charge,elastic,limit,range,strength,coefficient);
		
		if((x[i][j]+delta)>=0 && (x[i][j]+delta)<=1){

		x[i][j]+=delta;

		}
		
		else if((x[i][j]+delta)<0) {
	
			x[i][j]=0;

		}

		else if((x[i][j]+delta)>1) {
	
			x[i][j]=1;

		}

		//cout<<x[i][j]<<endl;

	}

	/*for (int m=0;m<n;m++) {

		cout<<x[m][j]<<endl;

	}*/
average[j]=average[j]/(static_cast<double>(n));	

}

ofstream giuseppe;

giuseppe.open("output.txt");

/*for(int i=0;i<n;i++){
	for(int j=0;j<t;j++){
		giuseppe<<x[i][j]<<" ";
		
	}
	
	giuseppe<<endl;
}*/

for(int i=0; i<t; i++){

	giuseppe<<average[i]<<endl;
}

giuseppe<<endl;

for(int i=0; i<n; i++){

	giuseppe<<x[i][t-1]<<endl;

}

giuseppe.close();

/*for(int i=0; i<n; i++){

	cout<<x[i][t-1]<<endl;

}*/



delete[] average;
delete[] range;
delete[] att; 
delete[] sel;
delete[] sigma_att;
delete[] sigma_sel;
delete[] prob;
delete[] charge;
delete[] threshold;
delete[] elastic;
delete[] weight;
delete[] coefficient;


for(int i = 0; i < n; i++){

	delete[] x[i]; 
}

delete[] x;


return 0;

}

