#include<iostream>
#include<armadillo>
#include<stdlib.h>
#include<omp.h>
#include<iomanip>

#define DONT_USE_QICLIB_NLOPT
//Include the path of the QIClib file
#include".../QIClib"

using namespace std;
using namespace arma;
using namespace qic;


//The spin quantum number 's'
double spin = 1.5;
    
//Corresponding dimension 'd' of the Hilbert space
int dim = (int)(2.0*spin + 1.0);

//Entry for the off-diagonal terms of spin-x and spin-y matrices.
double paulientry(double j, double s){
	double b;
	b = sqrt((s+j)*(s+1-j));
	return b/2;
}

//Spin-x matrix
cx_mat paulix(double s){
	int dim = (int)(2.0*s + 1.0);

	cx_mat sigma(dim, dim, fill::zeros);
    
    //Setting the off-diagonal terms of the spin-x matrix.	
	for (int p=0; p<dim-1; p++){
		sigma(p, p+1) += paulientry(s-double(p), s);
		sigma(p+1, p) += paulientry(s-double(p), s);
	}
	return sigma;
}//A "d x d" matrix.

//Spin-y matrix
cx_mat pauliy(double s){
	int dim = (int)(2.0*s + 1.0);
	complex<double> i(0,1);

	cx_mat sigma(dim, dim, fill::zeros);
	
    //Setting the off-diagonal terms of the spin-y matrix.	
	for (int p=0; p<dim-1; p++){
		sigma(p, p+1) += -i*paulientry(s-double(p), s);
		sigma(p+1, p) += i*paulientry(s-double(p), s);
	}
	return sigma;
}//A "d x d" matrix.

//Spin-z matrix
cx_mat pauliz(double s){
	int dim = (int)(2.0*s + 1.0);

	cx_mat sigma(dim, dim, fill::zeros);
	
    //Setting the diagonal terms of the spin-z matrix,
    //as per convention.	
	for (int p=0; p<dim; p++){
		sigma(p,p) += s-double(p);
	}
	return sigma;
}//A "d x d" matrix.




//Function to set-up a "\sum S_a \tensor S_a"  nearest neighbor chain
//where 'a' denotes 'x', 'y', or 'z'.
cx_mat chain(int L, cx_mat sigma){
	cx_mat spinchain(pow(dim,L), pow(dim,L), fill::zeros);
	cx_mat ten, iden = eye<cx_mat>(dim, dim);
	vector<cx_mat> sig;

	for (int p=0; p<L; p++) sig.push_back(iden);

	for (int p=0; p<L-1; p++){
		sig[p] = sigma; sig[p+1] = sigma;

		ten = kron(sig[0], sig[1]);
		for (int k=2; k<L; k++) ten = kron(ten, sig[k]); 

		spinchain += ten;

		sig[p] = iden, sig[p+1] = iden;
	}

	return spinchain;
}//Returns a matrix of size "L^d x L^d".


//Overloaded function for chain "\sum S_a \tensor S_a" for
//variable range.
cx_mat chain(int L, cx_mat sigma, int r){
	cx_mat spinchain(pow(dim,L), pow(dim,L), fill::zeros);
	cx_mat ten, iden = eye<cx_mat>(dim, dim);
	vector<cx_mat> sig;

	for (int p=0; p<L; p++) sig.push_back(iden);

	for (int p=0; p<L-r; p++){
		sig[p] = sigma; sig[p+r] = sigma;

		ten = kron(sig[0], sig[1]);
		for (int k=2; k<L; k++) ten = kron(ten, sig[k]); 

		spinchain += ten;

		sig[p] = iden, sig[p+r] = iden;
	}

	return spinchain;
}//Returns a matrix of size "L^d x L^d".


//Function to set-up a "\sum S_a" matrix.
cx_mat bfield(int L, cx_mat sigma){
	cx_mat sumsigma(pow(dim, L), pow(dim, L), fill::zeros);
	cx_mat ten, iden = eye<cx_mat>(dim, dim);
	vector<cx_mat> sig;

	for (int p=0; p<L; p++) sig.push_back(iden);

	for (int p=0; p<L; p++){
		sig[p] = sigma;

		ten = kron(sig[0], sig[1]);
		for (int k=2; k<L; k++) ten = kron(ten, sig[k]);

		sumsigma += ten;
		sig[p] = iden;
	}
	return sumsigma;
}//Returns a matrix of size "L^d x L^d".


//Function to calculate a thermal state of the
//provided Hamiltonian. Basically calculate "exp(-\beta H)".
cx_mat thermalstate(cx_mat H, double b){
	cx_mat state;
	int l = sqrt(H.size());

	vec eigval; cx_mat eigvec;
	eig_sym(eigval, eigvec, H);

	cx_mat Yexp(l, l, fill::zeros);

	for (int s=0; s<l; s++){
		double expo = 0.0;
		for (int e=0; e<l; e++){
			if (s!=e){
				expo += exp(-b*(eigval(e) - eigval(s)));
			}
		}
		Yexp(s,s) = 1.0/(1.0 + expo);
	}	
	state = eigvec*Yexp*eigvec.t();

	return state;
}//Returns a matrix of size "L^d x L^d".


//Setting-up projector of "(|0> + |(d-1)>)/\sqrt(2)" on the 
//first particle of the chain.
cx_mat projectx(int L){
	cx_mat z = pauliz(spin);
	cx_mat iden = eye<cx_mat> (pow(dim, L-1), pow(dim, L-1));
	vec eigval; cx_mat eigvec;
	eig_sym(eigval, eigvec, z);

	complex<double> ax(0.5,0.0),bx(0.5,0.0);
	cx_mat projector = ax*eigvec.col(0) + bx*eigvec.col(dim-1);

	double norm = sqrt(real(as_scalar(projector.t()*projector)));
	projector = (1/norm)*projector;
	projector = projector*projector.t();
	//cout << projector << endl;
	
	return kron(projector, iden);
}//Returns a matrix of size "L^d x L^d".

//Setting-up projector of "(|0> + i|(d-1)>)/\sqrt(2)" on the 
//first particle of the chain.
cx_mat projecty(int L){
	cx_mat z = pauliz(spin);
	cx_mat iden = eye<cx_mat> (pow(dim, L-1), pow(dim, L-1));
	vec eigval; cx_mat eigvec;
	eig_sym(eigval, eigvec, z);

	complex<double> ax(0.5,0.0),bx(0.0,0.5);
	cx_mat projector = ax*eigvec.col(0) + bx*eigvec.col(dim-1);

	double norm = sqrt(real(as_scalar(projector.t()*projector)));
	projector = (1/norm)*projector;
	projector = projector*projector.t();
	//cout << projector << endl;
	
	return kron(projector, iden);
}//Returns a matrix of size "L^d x L^d".



//Main Function
int main(int argc, char* argv[]){
    
    //User-input for the size of the chain.	
	int L = strtol(argv[1], NULL, 10);	

	//User-input for the fall-off rate.
	double alpha = atof(argv[2]); 

	//User-input for the coordination number.	
	int Z = atof(argv[3]);

    //User-input for the magnetic field strength in x-direction.
	double hx = atof(argv[4]);

    //User-input for the beginning time for evaluation.
	double tinit = atof(argv[5]); 

	//Defining the spin-x and spin-z matrices using 
    //function 'paulix' and 'pauliz'.
	cx_mat sigmax = paulix(spin), sigmaz = pauliz(spin);    
    
    //A 'd x d' identity matrix.
	cx_mat iden = eye<cx_mat>(dim,dim);
    
    //Complex number 'iota'.
	complex<double> i(0,1);
    
    //Coupling strength for spin-spin interaction, magnetic field 
    //strength(parameter) to be encoded onto the state, inverse temperature, 
    //respectively.
	double J = -1.0, omega = 1e-06, beta=10;

    //Small-increment in parameter, for calculating 
    //the derivative w.r.t. parameter.
	const double domega = 1e-07;

    //one step and two step increment in parameter for 
    //four-point derivative of the probability. 	
	double omegap = omega + domega, omegam = omega - domega;
	double omegap2 = omega + 2*domega, omegam2 = omega - 2*domega;
    
    //The value of 'pi'.
    double pi = 4.0*atan(1.0);

    //Parameter encoding time.
	double Tint = 500*pi;
    
    //The Kac-factor.
	double A = 0;
	for (int p=1; p<=1; p++){
		A += 1.0/pow(double(p), alpha);
	}

	/*
	cx_mat zzchain(pow(dim,L),pow(dim,L),fill::zeros);
	for(int p = 1; p <= Z; p++){
		zzchain += chain(L,sigmaz,p)/pow(p, alpha);
	}*/
    
    //Ising interaction in 'z' direction.	
    //A matrix of 'L^d x L^d' size.
	cx_mat zzchain(pow(dim,L),pow(dim,L),fill::zeros);
	zzchain += chain(L, sigmaz, 1);
	//zzchain += chain(L, sigmaz);
	//zzchain+=chain(L, sigmaz, sigmax, 2)/pow(2, alpha);
    
    //The magnetic field operator in 'z' direction. 
    //A matrix of "L^d x L^d" size.
	cx_mat fieldz = bfield(L, sigmaz);
    
    //The magnetic field operator in 'x' direction.
    //A matrix of "L^d x L^d" size.
	cx_mat fieldx = bfield(L, sigmax);
    
    //Ising Hamiltonian with a Kac-factor and coupling strength.	
	cx_mat hising = (J/A)*zzchain;
    
    //The transverse field Ising (TFI) Hamiltonian.	
	cx_mat htfi = hising + hx*fieldx;
    
    //The magnetic field operator with the magnetic field strength.
	cx_mat homega = omega*fieldz;

    //The magnetic field operator with the magnetic field strength 
    //with one step increment in parameter. 
	cx_mat homegap = omegap*fieldz;

	//The magnetic field operator with the magnetic field strength 
    //with one step decrement in parameter. 
	cx_mat homegam = omegam*fieldz;
    
    //The magnetic field operator with the magnetic field strength 
    //with two step increment in parameter. 
	cx_mat homegap2 = omegap2*fieldz;
    
    //The magnetic field operator with the magnetic field strength 
    //with two step decrement in parameter. 	
	cx_mat homegam2 = omegam2*fieldz;
    
    //Corresponding unitaries for paramter encoding.	
	cx_mat unitaryomega = expm_sym(homega, -i*Tint);
	cx_mat unitaryomegap = expm_sym(homegap, -i*Tint);
	cx_mat unitaryomegam = expm_sym(homegam, -i*Tint);
	cx_mat unitaryomegap2 = expm_sym(homegap2, -i*Tint);
	cx_mat unitaryomegam2 = expm_sym(homegam2, -i*Tint);

    //Thermal state of the transverse Ising chain.
	cx_mat rhobeta = thermalstate(htfi, beta);

    //Projector for the first particle of the chain.
	cx_mat px = projectx(L), py = projecty(L);	
    
    //Initial projector on the first particle of thermal state.
	cx_mat rho0 = px*rhobeta*px.t();
	rho0 = rho0/trace(rho0);
   
     
    double delta;         	
    int num_data_points=10000;
	vec deltaomega(num_data_points);	
	
    
    //The name of the data file.
	//string filename = to_string(spin).substr(0,4);
	//filename += "_offising_";
	//filename += to_string(int(tinit))+"_"+to_string(L).substr(0,4);
	//filename += "_"+to_string(alpha).substr(0,4)+"_";
	//filename += to_string(Z)+".dat";

	//Enter a suitable name for the file.
    string filename = "test.dat";
    
    //Object for the file output.
	ofstream file(filename);
	

	cout << "Calculating standard deviation in the estimation...\n";
	#pragma omp parallel for
	for (int p=0; p<num_data_points; p++){
		cx_mat rhot, rhott, rhotomega, rhotomegap, rhottomegap;
		cx_mat rhotomegam, rhottomegam;
		cx_mat rhotomegam2, rhottomegam2;
		cx_mat rhotomegap2, rhottomegap2;
    
		double prob, probomegap, probomegam, dprob, delta;
        double probomegap2, probomegam2;

        //Time for which the state is evolved with TFI.
		double tstar = double(p)/100 + tinit;
		cx_mat unitary = expm_sym(htfi, -i*tstar);

        //Evolution of projected state with TFI.	
		rhot = unitary*rho0*unitary.t();
	    
        //Evolution with magnetic field and its 
        //one and two step increment and decrement. 	
		rhotomegap = unitaryomegap*rhot*unitaryomegap.t();
		rhotomegam = unitaryomegam*rhot*unitaryomegam.t();
		rhotomegap2 = unitaryomegap2*rhot*unitaryomegap2.t();
		rhotomegam2 = unitaryomegam2*rhot*unitaryomegam2.t();
		rhotomega = unitaryomega*rhot*unitaryomega.t();
    
        //Evolution with reversed time uintaries. 
		rhott = unitary.t()*rhotomega*unitary;
		rhottomegap = unitary.t()*rhotomegap*unitary;
		rhottomegam = unitary.t()*rhotomegam*unitary;
		rhottomegap2 = unitary.t()*rhotomegap2*unitary;
		rhottomegam2 = unitary.t()*rhotomegam2*unitary;
	    
        //Calculating the probability for the parameter
        //and its increment/decrement for the derivative.	
		prob = real(trace(rhott*py));
		probomegap = real(trace(rhottomegap*py));
		probomegam = real(trace(rhottomegam*py));	
		probomegap2 = real(trace(rhottomegap2*py));
		probomegam2 = real(trace(rhottomegam2*py));	
        
        //The four-point derivative at the defined value. 
        dprob = (probomegam2 - 8*probomegam + 8*probomegap - probomegap2);
        dprob = dprob/(6.0*(omegap - omegam));

        //Two-point derivative.
		//dprob = (probomegap - probomegam)/(omegap - omegam);
		//dprob = (probomegap - prob)/(omegap - omega);	

        //The uncertainty.
		delta = sqrt(prob - prob*prob)/abs(dprob);
		deltaomega(p) = delta*sqrt(Tint + 2*tstar);
	}	

	cout << "Caculation complete.\n";	
	
	cout << "Minimum std. dev.: " << deltaomega.min() << endl;
	cout << "At " << deltaomega.index_min() << endl;

    //Output on the data file.	
	double tstar = 0;
	for (int p=0; p<num_data_points; p++){
		tstar = double(p)/100 + tinit;
		file << tstar << "	" << deltaomega(p) << endl;
	}

	file.close();
	cout << filename << " completed.\n";

	return 0;
}


