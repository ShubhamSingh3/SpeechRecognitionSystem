
// 204101054_HMM.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include "vector"
#include "fstream"
#include "string"
#include "sstream"
#include "list"
#include "cmath"

using namespace std;

#define frameSize 320
#define thresholdEnergy 10000
#define PI 3.14159265
#define N 10 //number of states in system
#define M 32 //number of dustinct observation symbols
#define ModelIterate 25 //number of re-estimations preferred


//Global variables
long double pi[N];
long double A[N][N];
long double B[N][M];
long double pi_bar[N];
long double a_bar[N][N];
long double b_bar[N][M];
long double P; //to store the probability in apha, beta, gamma computations
long double **alpha;
long double **beta;
long double **gamma;
long double **delta;
int **psi;
long double ***zeta;
long double Pstar; //maximum probability in viterbi algo.
int *Qstar;//
int n = 0;
int totalFrame = 0;
int startMarker = 0;
int endMarker = 0;
vector<long double> sample;
vector<long double> a;
vector<int> o;
long double codebook[32][12];
long double r[13];
long double w[13];
vector<long double> energy(30000,0);
int T; //to store number of time frames

struct model //to represent each model
{
    long double pi[N];
	long double A[N][N];
	long double B[N][M];
};

struct model models[10];


// Function to calculate Tokhura Distance
long double tokhuraDistance(vector<long double> testCoeff, long double refCoeff[]){

	long double sum = 0;
	long double weight[13] = {0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

	for(int i = 1; i <= 12; i++){

		long double z = testCoeff[i] - refCoeff[i-1];
		sum +=  weight[i] * ( z * z );

	}

	return sum;
}


/*
	Function for pre-processing on the speech signal. This will apply DC Shift on signal and normalize the signal with maximum amplitude 10,000. 
	After normalization it apply framing on signal and find enegry of each frame.
	Frame Size = 320
*/
void preProcessing(){

	long double dcShift = 0.0;
	long double maxAmplitude = 0.0;
	long double factor = 0.0;
	long double sum = 0.0;
	startMarker = 0;
	endMarker = 0;
	//Sample size calculated
	n = sample.size();
	
	//DC Shfit Calculation
	for(int i = 0; i < 640; i++){

			sum += sample[i];
	}

	dcShift = sum / 640; //DC Shift Calculated
	//cout << "DC Shift: " << dcShift << endl; 
	//For maximum amplitude value 
	maxAmplitude = 0;

	for(int j = 0; j < n; j++){

		//DC Shift applied 
		sample[j] -= dcShift;

		//if current apmlitude is higher than maxAmplitude
		if(maxAmplitude < abs(sample[j])){

			maxAmplitude = abs(sample[j]); //maxAmplitude updated with higher one
		}
	}
	
	//Normalization
	factor = (long double)10000 / maxAmplitude;
	//cout <<"factor: " << factor << endl;
	//Normalizing whole signal
	for(int j = 0; j < n; j++){

		sample[j] *= factor;
		
	}

	
	
	totalFrame  = n / frameSize; //Total number of frame calculated

	//Framing and energy calculated for each frame
	for(int i = 0; i < totalFrame; i++){

		energy[i] = 0;

		for(int j = 0; j < 320; j++){

			int k = i*frameSize+j;
			energy[i] += sample[k]*sample[k];
			
			
		}

		energy[i] = energy[i]/frameSize;
	
	}
	
	int k = 0;

	for(k = 0; k < totalFrame-4; k++){

		long double ste = energy[k] + energy[k+1] + energy[k+2] + energy[k+3];
		ste = ste/4.0;

		if(ste > 50000.0){
			startMarker = k*frameSize;
			k++;
			break;
		}

	}
	
	for(; k < totalFrame - 4; k++){

		long double ste = energy[k] + energy[k+1] + energy[k+2] + energy[k+3];
		ste = ste/4.0;

		if(ste < 50000.0){
			endMarker = k * frameSize;
			break;
		}
	}
	
	/*

	for(k = 0; k < n-4; k++){

		long double ste = (sample[k]*sample[k]) + (sample[k+1]*sample[k+1]) + (sample[k+2]*sample[k+2]) + (sample[k+3]*sample[k+3]);
		ste = ste/4.0;

		if(ste > 1000){
			startMarker = k;
			break;
		}
	}

	for(; k < n-4; k++){

		long double ste = (sample[k]*sample[k]) + (sample[k+1]*sample[k+1]) + (sample[k+2]*sample[k+2]) + (sample[k+3]*sample[k+3]);
		ste = ste/4.0;

		if(ste < 1000){
			endMarker = k;
			break;
		}
	}
	*/
	//cout << "startMarker: " << startMarker <<endl;
	//cout <<  "endMarker: " << endMarker << endl; 
	
}

/*
	Function to compute Levinson Durbin coefficient from autocorrelation vector
*/
vector<long double> levinsonDurbin(long double r[], int n){

	//int n = r.size();
	int p = n-1;
	double sum = 0;
	vector<vector<long double> > alpha (n, vector<long double> (n,0));
	vector<long double> k(n,0);
	vector<long double> e(n,0);
	vector<long double> a(n,0);

	e[0] = r[0];
	
	for(int i = 1; i <= p; i++){

		sum = 0;

		for(int j = 1; j < i; j++){

			sum += alpha[i-1][j] * r[i-j];
		}

		k[i] = (r[i] - sum) / e[i-1];
		alpha[i][i] = k[i];

		for(int j = 1; j < i; j++){

			alpha[i][j] = alpha[i-1][j] - k[i] * alpha[i-1][i-j];
		}

		e[i] = (1 - (k[i]*k[i])) * e[i-1];
	}

	for(int i = 1; i <=p ; i++){

		a[i] = alpha[p][i];
	}

	return a;
}


// Function to compute  Cepstral Coefficient
vector<long double> cepstralCoefficient(){

	long double sum = 0.0;
	vector<long double> c(13,0);
	for(int i = 0; i < 13; i++){
		c[i] = 0.0;
	}
	c[0] = logl(r[0]);
	for(int m = 1; m <= 12; m++){

				sum = 0.0;

				for(int k = 1; k < m; k++){

					sum += (((long double)k/m))*c[k]*a[m-k]; 

				}
				
				c[m] = a[m] + sum;
				
		}

	return c;

}


void getObservationSequence(){

	
	int max = 0;
	long double sum = 0.0;
	vector<long double> coeff;
	long double min_dist = 0;
	int codebook_entry = 0;

	/*
	//Find frame with maximum energy
	for(int i = 0; i < n; i += 80){

			long double ste = 0;
			for(int j = i; j < i + 320; j++){

				ste += sample[j]*sample[j];
			}
			//cout << ste << "\t";
			if(ste > 500){

				max = i;
				break;

			}

	}
	*/
	//calcuating reference of 5 frame around the frame with highest energy
	//taking 2 frame before that highest energy frame and 2 frame after that frame for reference
	
	
	for(int i = startMarker; i < endMarker; i += 80){

		
		//int index = i*frameSize;
	
		for(int k = 0; k <= 12; k++){
			
			sum = 0;
	
			//calculating autocorrelation
			for(int m = 0; m <= 319 - k; m++){

				sum += sample[i+m]*sample[i+m+k];

			}
			
			
			r[k] = sum/320.0;

		}
		
		//levinson durbin coefficient for the autocorrelation
		a = levinsonDurbin(r, 13);

		//cepstral coefficient calculated 
		coeff = cepstralCoefficient();
		

		//applied raised sine window on cepstral coefficient
		for(int m = 1; m <= 12; m++){

			coeff[m] = coeff[m] * w[m];

		}	
		
		min_dist = tokhuraDistance(coeff, codebook[0]);
		codebook_entry = 0;
		for(int i = 1; i < 32; i++){
			
			long double dist = tokhuraDistance(coeff, codebook[i]);
			if(dist < min_dist){
				min_dist = dist;
				codebook_entry = i;
			}
		}
		//observationSeq[frameCount-1] = codebook_entry + 1;
		o.push_back(codebook_entry);
		//cout << o.back() << "\t";

	}

}


void getCodeBook(){

	string line = "";
	string temp;
	ifstream myfile;
	string fname = "codebook.txt";
	myfile.open (fname); 
		
	if(myfile.is_open()){
		
		int k = 0;
		while(!myfile.eof()){

			int j = 0;
			getline(myfile,line);
			stringstream ss(line);
			while(getline(ss,temp,'\t')){

				codebook[k][j++] = stod(temp.c_str());

			}
			k++;

		}
	}
	/*
	for(int i = 0; i < 32; i++){
		for(int j = 0; j < 12; j++){
			cout << codebook[i][j] << "\t";
		}
		cout << endl;
	}
	*/

}

void feed_forward_model() //this is biased model which works in lazy way
{
	//cout << "In forward procedure" << endl;
	//out << "In forward procedure" << endl;
	for(int i=0;i<N;i++) //assign the given values
		if(i==0) //make the first state as starting state
			pi[i]=1.0;
		else
			pi[i]=0;
    for(int i=0;i<N;i++) //assign the given values for transition probability distribution
        for(int j=0;j<N;j++)
			if(i==j&&i!=N-1)
				A[i][j]=0.8; //probability of being in current state
			else if(i==j&&i==N-1)
				A[i][j]=1; //forcing to end the transition in final state
			else if(j==i+1)
				A[i][j]=0.2; //probability to transit to next immediate state
			else
				A[i][j]=0; //probability to move to non immediate states
    for(int i=0;i<N;i++) //assign the given values for observation symbol probability distribution
        for(int j=0;j<M;j++)
            B[i][j]=1.0/M;
}


long double forward_proc() //alpha computation
{
    for(int i=0;i<N;i++) //initialization
        alpha[0][i]=pi[i]*B[i][o[0]];
	
    for(int t=0;t<T-1;t++) //induction
    {
        for(int j=0;j<N;j++)
        {
            double sum=0;
            for(int i=0;i<N;i++)
            {
                sum+=alpha[t][i]*A[i][j];
            }
            alpha[t+1][j]=sum*B[j][o[t+1]];
        }
    }
    P=0;
	
    for(int i=0;i<N;i++) //estimate what is the probability that the observation sequence is from the current model
	{
        P+=alpha[T-1][i];
		//cout <<"alpha: " << alpha[T-1][i] << endl;
	}
	
    return P;
}

 long double forward_proc2(long double AA[N][N], long double BB[N][M]) //alpha computation
{
    for(int i=0;i<N;i++) //initialization
        alpha[0][i]=pi[i]*BB[i][o[0]];
    for(int t=0;t<T-1;t++) //induction
    {
        for(int j=0;j<N;j++)
        {
            double sum=0;
            for(int i=0;i<N;i++)
            {
                sum+=alpha[t][i]*AA[i][j];
            }
            alpha[t+1][j]=sum*BB[j][o[t+1]];
        }
    }
    P=0;
    for(int i=0;i<N;i++) //estimate what is the probability that the observation sequence is from the current model
	{
        P+=alpha[T-1][i];
		//cout <<"alpha: " << alpha[T-1][i] << endl;
	}
    return P;
}


double backward_proc() //beta computation
{
    for(int i=0;i<N;i++) //initialization
        beta[T-1][i]=1;
    for(int t=T-2;t>=0;t--) //induction
    {
        for(int i=0;i<N;i++)
        {
            long double sum=0;
            for(int j=0;j<N;j++)
            {
                sum+=A[i][j]*B
					
					[j][o[t+1]]*beta[t+1][j];
            }
            beta[t][i]=sum;
        }
    }
    P=0;
    for(int i=0;i<N;i++) //this is for self assessment and not required as there is no further use
    {
        P+=beta[0][i];
		//cout << beta[0][i] << endl;
    }
   return P;
}

double gamma_proc() //gamma computation
{
	int *q,argmax=0; //q--> store the state which has maximum probability of occurence at time t.
	q=new int[T];
	long double devider=0; //used as devider in baye's theorem for computation of gamma
	for(int t=0;t<T;t++)
	{
		for(int i=0;i<N;i++) //compute it once for t
		{
			devider+=alpha[t][i]*beta[t][i];
		}
		argmax=0;
		for(int i=0;i<N;i++)
		{
			gamma[t][i]=alpha[t][i]*beta[t][i]/devider;
			if(gamma[t][argmax]<gamma[t][i])
				argmax=i;
		}
		q[t]=argmax;
		devider=0;
	}
	P=1;
	for(int t=0;t<T;t++)
	{
		P*=gamma[t][q[t]];
		//cout << gamma[t][q[t]] << endl;
	}
	return P;
}

long double viterbi_proc() //delta and psi computation
{
	Qstar=new int[T];
	int argmax=0; // to store the argument which gives maximum probability
	for(int i=0;i<N;i++) //initialization
	{
		delta[0][i]=pi[i]*B[i][o[0]];
		//cout  << delta[0][i] << endl;
		psi[0][i]=-1; //indicates, no state is yet assigned
	}
	for(int t=1;t<T;t++) //recursion over time sequence
	{
		for(int j=0;j<N;j++) // to store maximum probabilities
		{
			argmax=0; //assume the first state gives maximum probability
			for(int i=1;i<N;i++)
			{
				if(delta[t-1][i]*A[i][j] > delta[t-1][argmax]*A[argmax][j]) //checking for maximum score
					argmax=i; //argument which gives maximum score
			}
			delta[t][j]=delta[t-1][argmax]*A[argmax][j]*B[j][o[t]];
			psi[t][j]=argmax; //largest argument index which gives maximum probability
		}
	}
	argmax=0;
	for(int i=1;i<N;i++) //to find the argmax for last time frame
		if(delta[T-1][i] > delta[T-1][argmax])
			argmax=i;
	Pstar=delta[T-1][argmax];
	//cout << "state sequence with back tracking" << endl;
	Qstar[T-1]=argmax;
	//cout << T-1+1 << " ---> " << Qstar[T-1]+1 << endl;
	for(int t=T-2;t>=0;t--) //back tracking the path
	{
		Qstar[t]=psi[t+1][Qstar[t+1]];
		//cout << t+1 << " --> " << Qstar[t]+1 << endl;
	}
	return Pstar;
}



long double viterbi_proc_alt() //alternate viterbi method for delta and psi computation
{
	long double pii[N];
	long double AA[N][N];
	long double BB[N][M];
	Qstar=new int[T];
	int argmax=0; // to store the argument which gives maximum probability

	//preprocessing
	for(int i=0; i<N;i++){
		pii[i] = logl(pi[i]);
		for(int j=0; j<N; j++){

			AA[i][j] = logl(A[i][j]);
		}
	}
	for(int i=0; i < N;i++){
		for(int j=0; j< T; j++){
			BB[i][o[j]] = logl(B[i][o[j]]);
		}
	}


	for(int i=0;i<N;i++) //initialization
	{
		delta[0][i]=pii[i]+BB[i][o[0]];
		//cout  << delta[0][i] << endl;
		psi[0][i]=-1; //indicates, no state is yet assigned
	}
	for(int t=1;t<T;t++) //recursion over time sequence
	{
		for(int j=0;j<N;j++) // to store maximum probabilities
		{
			argmax=0; //assume the first state gives maximum probability
			for(int i=1;i<N;i++)
			{
				if(delta[t-1][i]+AA[i][j] > delta[t-1][argmax]+AA[argmax][j]) //checking for maximum score
					argmax=i; //argument which gives maximum score
			}
			delta[t][j]=delta[t-1][argmax]+AA[argmax][j]+BB[j][o[t]];
			psi[t][j]=argmax; //largest argument index which gives maximum probability
		}
	}
	argmax=0;
	for(int i=1;i<N;i++) //to find the argmax for last time frame
		if(delta[T-1][i] > delta[T-1][argmax])
			argmax=i;
	Pstar=delta[T-1][argmax];
	//cout << "state sequence with back tracking" << endl;
	Qstar[T-1]=argmax;
	//cout << T-1+1 << " ---> " << Qstar[T-1]+1 << endl;
	for(int t=T-2;t>=0;t--) //back tracking the path
	{
		Qstar[t]=psi[t+1][Qstar[t+1]];
		//cout << t+1 << " --> " << Qstar[t]+1 << endl;
	}
	return Pstar;
}






void disp_state_seq() //display the state sequence followed
{
	std :: cout << "\nobserved state sequence: " << endl;
	for(int t=0;t<T;t++)
	{
		cout << Qstar[t] << '\t';
	}
	cout << endl;
}


void baum_welch_proc() //zeta computaion
{
	long double devider=0; //used as common devider just like in baye's theorem
	for(int t=0;t<T-1;t++) //repeat this for all T-1 state transitions
	{
		devider=0;
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
				devider+=alpha[t][i]*A[i][j]*B[j][o[t+1]]*beta[t+1][j];
		}
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
				zeta[t][i][j]=alpha[t][i]*A[i][j]*B[j][o[t+1]]*beta[t+1][j]/devider;
		}
	}
}


void re_estimation() //re-estimate transition probabilities and observation symbol probabilities
{
	long double numerator, denominator; //for re-estimation of transition probabilities
	for(int i=0;i<N;i++) //re-estimation of pi as pi_bar
		pi_bar[i]=gamma[0][i];
	for(int i=0;i<N;i++) //re-estimation of a as a_bar
	{
		for(int j=0;j<N;j++)
		{
			numerator=0;
			denominator=0;
			for(int t=0;t<T-2;t++)
			{
				numerator+=zeta[t][i][j];
				denominator+=gamma[t][i];
			}
			a_bar[i][j]=numerator/denominator;
		}
	}
	for(int j=0;j<N;j++) //re-estimation of b as b_bar
	{
		for(int k=0;k<M;k++)
		{
			numerator=0;
			denominator=0;
			for(int t=0;t<T;t++)
			{
				if(o[t]==k)
					numerator+=gamma[t][j];
			}
			for(int t=0;t<T-1;t++)
			{
				denominator+=gamma[t][j];
			}
			b_bar[j][k]=numerator/denominator;
		}
	}
}


void diplay_model(long double pi[N], long double a[N][N], long double b[N][M]) //display the whole model
{
	cout << "Initial State probabilities: " << endl;
	//out << "Initial State probabilities: " << endl;
	for(int i=0;i<N;i++)  //display initial state state distributions
	{
        cout << pi[i] << "\t";
		//out << pi[i] << "\t";
	}
	cout << "\ntransition probabilities: " << endl;
	//out << "\ntransition probabilities: " << endl;
    for(int i=0;i<N;i++) //display state transition probabilities
	{
        for(int j=0;j<N;j++)
		{
            cout << A[i][j] << "\t";
			//out << a[i][j] << "\t";
		}
		cout << endl;
		//out << endl;
	}
	cout << "observation symbol probability: " << endl;
	//out << "observation symbol probability: " << endl;
    for(int i=0;i<N;i++) //display observation symbol probability distribution
	{
        for(int j=0;j<M;j++)
		{
            cout << B[i][j] << "\t";
			//out << b[i][j] << "\t";
		}
		cout << endl;
		//out << endl;
	}
}

void replace_model() //replace old model (lambda) with new model (lambda-dash)
{
	cout << "replacing the old model with new model" << endl;
	//out << "replacing the old model with new model" << endl;
	for(int i=0;i < N;i++) //assign the given values
        pi[i]=pi_bar[i];
    for(int i=0;i < N;i++) //assign the given values for transition probability distribution
        for(int j=0;j < N;j++)
            A[i][j]=a_bar[i][j];
    for(int i=0;i < N;i++) //assign the given values for observation symbol probability distribution
        for(int j=0;j < M;j++)
            B[i][j]=b_bar[i][j];
}



void createCodeBook(){

	for(int i = 0; i <= 9; i++){


	}

}

int _tmain(int argc, _TCHAR* argv[])
{

	for(int m = 1; m <= 12; m++){

			w[m] = 1 + 6*sin ( PI * m / 2 );

	}

	//get code book
	getCodeBook();

	long double AA[20][N][N];
	long double BB[20][N][M];
	long double pii[20][N];

	ofstream myfile;
	string fname = "output.txt";
	myfile.open (fname);
	//get cepstral coefficient
	myfile <<"-----------------------------------------------TRAINING:----------------------------------------\n";
	myfile <<"File name with correponding P* value:\n\n";
	myfile <<"File Name\t\t\t   P*\n";
	for(long long l=0; l<=9; l++){
		for(long long k=1; k<=20; k++){

			sample.clear();
			string filename = "./digits/204101054_" + to_string(l) + "_" + to_string(k)+".txt";
			cout <<"FILE: " << filename << endl;
			ifstream inFile;
			inFile.open(filename); 
			long double x;
		
			if(inFile.is_open()){

				while ( inFile >> x ) {

					sample.push_back(x);

				}

				n = sample.size();           //size of signal
				//cout << "n = " << n << endl;
				//cout << "PreProcessing\n";
				preProcessing();             //pre processing of signal
				//cout << "PreProcessing Done\n";
				o.clear();
				getObservationSequence(); //compute coefficients for testing
				//cout << "getObservationSequence Done\n";
				T = o.size();
				//cout << "T :" << T << endl;

				//////////////////All declarations///////////////////////////////////////
				alpha=new long double*[T];
				for(int i=0;i < T;i++){

					alpha[i]=new long double[N];

				}

				beta=new long double*[T];
				for(int i=0;i < T;i++)
					beta[i]=new long double[N];

				gamma=new long double*[T];
				for(int i=0;i < T;i++)
					gamma[i]=new long double[N];

				delta=new long double*[T];
				for(int i=0;i < T;i++)
					delta[i]=new long double[N];

				psi=new int*[T];
				for(int i=0;i < T;i++)
					psi[i]=new int[N];

				zeta=new long double**[T];
				for(int i=0;i < T;i++)
				{
					zeta[i]=new long double*[N];
					for(int j=0;j < N;j++)
						zeta[i][j]=new long double[N];
				}

				feed_forward_model(); //biased transitions (matrices A,B,pi)
				cout << "feed_forward_model \n";
				///////////////////////////////////////////////////////////alpha, beta, gamma, zeta, psi computation//////////////////////////////////////////////////////////

				int counter=0;
				long double thrshold = pow(10.0,-60);
				long double old_Pstar = 0;
				Pstar = 1;
				//diplay_model(pi, A, B);
				do
				{
					//cout << "Iteration = " << counter << endl;
					cout  <<  "Probability of an utterance being from the model: " << forward_proc() << endl;
					//cout  <<  "Probability in backward procedure: " << backward_proc() << endl;
					//out  <<  "Probability in backward procedure: " << backward_proc() << endl;
					backward_proc();
					//cout  <<  "Probability in gamma procedure: " << gamma_proc() << endl;
					//out  <<  "Probability in gamma procedure: " << gamma_proc() << endl;
					gamma_proc();
					old_Pstar = Pstar;
					cout  <<  "Probability in viterbi procedure: " << viterbi_proc() << endl;
					//out  <<  "Probability in viterbi procedure: " << viterbi_proc() << endl;
					disp_state_seq();

					/////////////////////////////////////////////////////re-estimation solution///////////////////////////////////////////////////////////
					baum_welch_proc();
					re_estimation();

					//cout << "\n=====================MODEL before re-estimation================================" << endl;
					//out << "\n=====================MODEL before re-estimation================================" << endl;
					//diplay_model(pi, A, B);
					//cout << "\n=====================MODEL after re-estimation================================" << endl;
					//out << "\n=====================MODEL after re-estimation================================" << endl;
					//diplay_model(pi_bar, a_bar, b_bar);
					if(Pstar < old_Pstar)
						replace_model();
				}while(Pstar > thrshold && Pstar < old_Pstar);

				myfile << filename << "\t" << Pstar << "\n\n";
				//cout << "\n=====================FINAL MODEL after re-estimation================================" << endl;
				//diplay_model(pi_bar, a_bar, b_bar);
				//replace_model();

				
				//cout  <<  "Probability of an utterance being from the model: " << forward_proc() << endl;
				backward_proc();
				gamma_proc();
				//disp_state_seq();
				

				for(int i = 0; i < N; i++){

					pii[k-1][i] = pi[N];
					for(int j = 0; j < N; j++){

						AA[k-1][i][j] = A[i][j];
					}
				}

				for(int i = 0; i < N; i++){

					for(int j = 0; j < M; j++){

						BB[k-1][i][j] = B[i][j];
					}
				}

			}
			else{
				cout << "UNABLE TO OPEN FILE" << endl;
			}
			inFile.close();  //close file
			//Average all model
			for(int i = 0; i < N; i++){

					for(int j = 0; j < N; j++){

						long double sum = 0;
						for(int x = 0; x < 20; x++){

							sum += AA[x][i][j]; 
						}
						sum = sum/20.0;
						models[l].A[i][j] = sum;
					}
				}

				for(int i = 0; i < N; i++){

					for(int j = 0; j < M; j++){

						long double sum = 0;
						for(int x = 0; x < 20; x++){

							sum += BB[x][i][j]; 
						}
						sum = sum/20.0;
						models[l].B[i][j] = sum;
					}
				}

				for(int i = 0; i < N; i++){

					long double sum = 0;
					for(int x = 0; x < 20; x++){

							sum += pii[x][i]; 
					}
					sum = sum/20.0;
					models[l].pi[i] = sum;
				}


		}
	}

	/////////////////////////////////////////////////****TESTING********//////////////////////////////////////////////////////////////////////////
	int sum = 0, count= 0;
	myfile << "-------------------------------------------TESTING--------------------------------------------\n";
	myfile <<"File name with all 10 probailities and winning score\n\n";
	for(long long k = 0; k <=9; k++){
		count = 0;
		for(long long j = 21; j <= 30; j++){
			sample.clear();
			string filename = "./digits/204101054_" + to_string(k) + "_" + to_string(j) +".txt";
			cout <<"FILE: " << filename << endl;
			myfile << filename <<"\n";
			ifstream inFile;
			inFile.open(filename); 
			long double x;
		
			if(inFile.is_open()){

				while ( inFile >> x ) {

					sample.push_back(x);

				}

				n = sample.size();           //size of signal
				//cout << "n = " << n << endl;
				//cout << "PreProcessing\n";
				preProcessing();             //pre processing of signal
				//cout << "PreProcessing Done\n";
				o.clear();
				getObservationSequence(); //compute coefficients for testing
				T = o.size();
				alpha=new long double*[T];
					for(int i=0;i < T;i++){

						alpha[i]=new long double[N];

					}
				long double max_pstar = forward_proc2(models[0].A, models[0].B);
				int digit = 0;
				myfile << max_pstar << "\n";
				for(int i = 1; i <= 9; i++){

					long double temp = forward_proc2(models[i].A, models[i].B);
					myfile << temp << "\n";
					if(temp > max_pstar){
						max_pstar = temp;
						digit = i;
					}

				}
				myfile << "WINNING SCORE: " << max_pstar << "\n\n\n";
				cout << "digit : " << digit << endl;
				if(digit == k){
					count++;
				}
					
			}
			else{
				cout << "Unable to open file "<< endl;
			}
		}
		sum+=count;
		cout << "No. of correct detection: " << count << endl;
	}
	cout << "Total correct detection: " << sum << endl;
	myfile <<"\n\n ACCURACY: " << sum << " / 100\n"; 
	system("pause");
	return 0;
}