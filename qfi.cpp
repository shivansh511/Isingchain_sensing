#include<iostream>
#include<armadillo>
#include<stdlib.h>

#define DONT_USE_QICLIB_NLOPT
//Include the path of the QIClib file
#include".../QIClib"

using namespace std;
using namespace arma;
using namespace qic;

//Spin quantum number "s"
double spin = 1.5;
//Corresponding dimension "d" of the Hilbert space.
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
}

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
}

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
}



//Function to set-up a "\sum S_a \tensor S_a" chain
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



//Calculate the Qauntum Fisher Information of a state "rho"
//w.r.t. a generator "j".
double qfi(cx_mat rho, cx_mat j){
	double fisherinfo = 0, factor, prod;

    //Requires all the eigenvectors and eigenvectors of 
    //the state.
	vec eigval; cx_mat eigvec;
	eig_sym(eigval, eigvec, rho);
	int l = sqrt(rho.size());
	
	cx_mat m,n;
    
	for (int p=0; p<l; p++){
		m = eigvec.col(p);
		for (int q=0; q<l; q++){
        
            //Skip the step if the sum of eigenvalue pair is zero. 
            //Check definition of QFI.       
			if (eigval(p) + eigval(q) < 1e-03) continue;

			n = eigvec.col(q);
			prod = real(as_scalar((m.t()*j*n)*(n.t()*j*m)));
			factor = 2.0*pow(eigval(p) - eigval(q), 2)/(eigval(p) + eigval(q));

			fisherinfo += prod*factor;
		}
	}

	return fisherinfo;
}


//Main Function
int main(int argc, char* argv[]){
	
	int L = strtol(argv[1], NULL, 10);	
	//spin = 1.5;
	cx_mat sigmax = paulix(spin), sigmaz = pauliz(spin);

	//Complex number "iota".
	complex<double> i(0,1);
	
	double hx = 0.1, J = -1.0, beta=10;
	
    //double pi = 4.0*atan(1.0);
    double tstar = 0;

	cx_mat zzchain = chain(L, sigmaz);
	cx_mat fieldx = bfield(L, sigmax);
	cx_mat fieldz = bfield(L, sigmaz);
    
    //Transverse field Ising (TFI) chain.
	cx_mat htfi = J*zzchain + hx*fieldx;

	//Thermal state of TFI Hamiltonian.
	cx_mat rhobeta = thermalstate(htfi, beta);
	
	//Projector measurement on the first spin of chain.
    cx_mat px = projectx(L);
	cx_mat rho0 = px*rhobeta*px.t();
	rho0 = rho0/trace(rho0);

	//Unitary evlution with TFI Hamiltonian.
	cx_mat unitary = expm_sym(htfi, -i*tstar);
	cx_mat dunitary = expm_sym(htfi, -i*0.01);

	cx_mat rhot;

	vec fisherinfo(10000);	
    
	cout << "Calculating QFI w.r.t S_z ...\n";

	for (int p=0; p<10000; p++){
		//Increment in unitary after with time step.
		unitary = unitary*dunitary;
		//tstar += 0.01;
		
		//Evolving the state with unitary.
		rhot = unitary*rho0*unitary.t();

		//Calculating QFI
		fisherinfo(p) = qfi(rhot, fieldz);
		
	}	

	cout << "Caculation complete.\n";	
	//cout << tstar << endl;
	
	cout << "Maximum QFI w.r.t. S_z: " << fisherinfo.max() << endl;

	//string filename = to_string(spin)+"_qfi_"+to_string(L)+".dat";
	
	//Enter a suitable name for the data file.
	string filename = to_string(spin).substr(0,4)+"_qfi_test.dat";
	
	ofstream file(filename);
	
	for (int p=0; p<10000; p++){
		tstar += 0.01;
		file << tstar << "	" << fisherinfo(p) << endl;
	}

	file.close();
	cout << filename << " completed.\n" << endl;

	return 0;
}

