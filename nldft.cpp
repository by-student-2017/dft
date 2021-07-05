#include <fstream>   // for file in and out
#include <iostream>  // for cout
#include <cmath>     // for log, exp
#include <sstream>   // for read parameters

//#include "maxwell_construction.h"

using namespace std;

//non-local smoothed density approximation：SDA
//non-local density functional theory（NLDFT)
//Reference: https://www.j-ad.org/adsorption_news/30_1.pdf

// Note
// This routine assumes that rho, etc is same value in x-y plane.
// Because of cut off (rc), it calculate circle in x-y plane.
// Units are fundamentaly [K] and [nm] in this routine.

// There are many imperfections, so I hope someone can make it better with a CC0 license. 
// It seems that this code is the first in the world at present (2021/7/5) to be released on CC0 even in NLDFT. 

// compiling: c++ nldft.cpp
// usage: ./a.out

// debag mode
// compiling: c++ nldft.cpp -g3 -ggdb
// run: gdb ./a.out
//      (gdb) r
//      (gdb) backtrace

// Adsorbent 
//double H = 1.00; //distace of slit [nm]
//double sigma_ss = 0.34; // [nm]
//int nstep = 100;
//double w_pw = (H-sigma_ss); // pore width [nm]
//double dr = w_pw/double(nstep);
double H;
double sigma_ss;
int nstep;
double w_pw;
double dr;

// assume rho is same value in x-y plane.
// cylinder and normalization, because of cut off (rc).
//int nrmesh = 20; //rho_si and xi function
int nrmesh;

// iteration of rho
//int cycle_max = 50;
int cycle_max;

//Carbon dioxide 253.9  [K](epsilon), 0.3454 [nm](sigma), 0.3495 [nm](d_hs)
//Argon          118.05 [K](epsilon), 0.3305 [nm](sigma), 0.3390 [nm](d_hs)
//Nitrogen        94.45 [K](epsilon), 0.3575 [nm](sigma), 0.3575 [nm](d_hs)
//extern double epsilon_ff = 94.45;
//double sigma_ff = 0.3575;
//extern double d_hs = 0.3575; // Maxwell_construction()
//double rc = 1.28; // [nm],cut off, (12.8 [A])
double epsilon_ff;
double sigma_ff;
double d_hs;
double rc;

//double rm = std::pow(2.0,1.0/6.0)*sigma_ff; //minimum position of LJ
//double rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
double rm;

// Carbon dioxide/Carbon slit 81.5  [K](epsilon), 0.3430 [nm](sigma)
// Nitrogen/Carbon slit       53.72 [K](epsilon), 0.3508 [nm](sigma)
//double epsilon_sf = 53.72; // [K] 
//double sigma_sf = 0.3508; // [nm]
double epsilon_sf;
double sigma_sf;

// slit pore (graphite)
//double delta = 0.335; // [nm]
//double rho_ss = 114.0; // [nm^-3], [molecules/nm3]?
double delta;
double rho_ss;

//double m = 14.0067*2.0/(6.02214076e23)/1000; // N2 = 4.65173e-26 [kg]
//double m = 4.65173e-26; //[kg] (N2) (e.g., Ar = 6.63e-26 [kg])
double m;
double k = 1.0;
double kb = 1.38e-23; //[J/K] (8.61733262e-5 [eV/K])
//extern double T = 77.347; //[K]
double T;
double h = 6.63e-34; //[Js] (4.135667696e-15 [eVs])
// thermal de Broglie wavelength
//extern double lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm], Maxwell_construction()
double lam;
// Ref: https://www1.doshisha.ac.jp/~bukka/lecture/statistic/pdftext/std-07.pdf

// alpha = integal phi_att * -1.0
//extern double alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
//	( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
double alpha;

void read_parameters(void){
	std::ifstream ifs("parameters.txt");
	std::string str;
	double num[20];
	int i,j;
	j = 0;
	while(getline(ifs,str)){
		std::string tmp;
		std::istringstream stream(str);
		i = 0;
		while(getline(stream,tmp,'=')){
			if (i == 1){
				num[j] = atof(tmp.c_str());
				//std::cout<< num[j] << std::endl;
			}
			i++;
		}
		j++;
	}
	//
	H = num[0]; //distace of slit [nm]
	sigma_ss = num[1]; // [nm]
	nstep = int(num[2]);
	nrmesh = int(num[3]);
	cycle_max = int(num[4]);
	epsilon_ff = num[5]; // [K]
	sigma_ff = num[6]; // [nm]
	d_hs = num[7]; // [nm]
	rc = num[8]; // [nm],cut off, (12.8 [A])
	epsilon_sf = num[9]; // [K]
	sigma_sf = num[10]; // [nm]
	delta = num[11]; // nm
	rho_ss = num[12]; // [nm^-3], [mulecules/nm3]
	m = num[13]; // [kg]
	T = num[14]; // [K]
	
	w_pw = (H-sigma_ss); // pore width [nm]
	dr = w_pw/double(nstep);
	rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
	
	// thermal de Broglie wavelength
	lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm], Maxwell_construction()
	
	// alpha = integal phi_att * -1.0
	alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
		( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
	
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "thermal de Broglie wavelength = " << lam << " [nm]" << std::endl;
	std::cout << "integal phi_att * -1.0 = alpha = " << alpha << std::endl;
}

double ingegral_simpson(double *f, int n, double dx){
	if( (n+1)%2 == 1 ){
		std::cout << "Error, plase change number of data to even" << std::endl;
	}
	double sum;
	sum = f[0] + f[n];
	int i;
	for(i=1; i<n; i+=2){
		sum += 4.0 * f[i];
	}
	for(i=2; i<n; i+=2){
		sum += 2.0 * f[i];
	}
	return (dx/3.0)*sum;
}

//Barker-Henderson model
double d_bh(void){
	double epsilon = 94.45; //dummy
	double sigma = 0.3575; //dummy
	//Lstoskie et al.,
	double xi1 = 0.3837;
	double xi2 = 1.035;
	double xi3 = 0.4249;
	double xi4 = 1.0;
	double d_bh_out;
	d_bh_out = (xi1*k*T/epsilon+xi2)/(xi3*k*T/epsilon+xi4)*sigma;
	return d_bh_out;
}

double phi_att(double r){
	double e;
	// WCA (Weeks-Chandler-Anderson) type
	if (r < rm){
		e = -1.0*epsilon_ff;
	}else if (rm <= r && r <= rc){
		// Lennard-Jones（LJ) potential
		e = 4.0*epsilon_ff*( std::pow((sigma_ff/r),12.0) - std::pow((sigma_ff/r),6.0) );
	}else if (rc < r){
		e = 0.0;
	}
	//std::cout << e << std::endl;
	return e;
}

// Tarazona theory, Percus-Yevick approximation
double wi(double r, int i){
	double wi_out;
	switch(i){
		case 0:
			if (r <= d_hs){
				wi_out = 3.0/(4.0*M_PI*std::pow(d_hs,3.0));
			} else {
				wi_out = 0.0;
			}
			break;
		case 1:
			if (r <= d_hs) {
				wi_out = 0.475-0.648*(r/d_hs)+0.113*std::pow((r/d_hs),2.0);
			} else if (d_hs < r && r <= 2.0*d_hs) {
				wi_out = 0.288*(d_hs/r)-0.924+0.764*(r/d_hs)-0.187*std::pow((r/d_hs),2.0);
			} else {
				wi_out = 0.0;
			}
			break;
		case 2:
			if (r <= d_hs) {
				wi_out = 5.0*M_PI*std::pow(d_hs,3.0)/144.0 * (6.0-12.0*(r/d_hs)+5.0*std::pow((r/d_hs),2.0));
			} else {
				wi_out = 0.0;
			}
			break;
		default:
			std::cout << "Error: " << i << std::endl;
			break;
	}
	return wi_out;
}

double rho_si(double *rho, double r1, double *r, int i){
	int j,k;
	double rho_si_out;
	double ra;
	rho_si_out = 0.0;
	double rho_si_int_j[nstep];
	double rho_si_int_k[nrmesh];
	for (j=0; j<nstep; j++) {
		for (k=0; k<=nrmesh; k++) {
			ra = std::pow((r1-r[j]),2.0) + std::pow((double(k)*rc/double(nrmesh)),2.0);
			ra = std::pow(ra,0.5);
			//std::cout << ra << std::endl;
			//
			//rho_si_out = rho_si_out + rho[j]*wi(std::abs(r1-r[j]),i)*(4.0*M_PI*r[j]*r[j])*dr;
			//rho_si_out = rho_si_out + rho[j]*wi(ra,i)*(2.0*M_PI*(double(k)*rc/double(nrmesh))*(rc/double(nrmesh)))*dr;
			rho_si_int_k[k] = rho[j]*wi(ra,i)*(2.0*M_PI*(double(k)*rc/double(nrmesh)));
		}
		//ingegral_simpson(double *f, int n, double dx)
		rho_si_int_j[j] = ingegral_simpson(rho_si_int_k, nrmesh, (rc/double(nrmesh)));
	}
	//ingegral_simpson(double *f, int n, double dx)
	rho_si_out = ingegral_simpson(rho_si_int_j, nstep, dr);
	//
	//rho_si_out = rho_si_out / (M_PI*std::pow((rc),2.0)) / (nstep*dr);
	return rho_si_out;
}

double rho_s(double *rho, double r1, double *r){
	double rho_den1, rho_den2, rho_s_out;
	rho_den1 = std::pow((1.0 - rho_si(rho,r1,r,1)),2.0);
	rho_den2 = std::pow((rho_den1 - 4.0*rho_si(rho,r1,r,0)*rho_si(rho,r1,r,2)),0.5);
	rho_s_out = 2.0*rho_si(rho,r1,r,0)/(1.0 - rho_si(rho,r1,r,1)+rho_den2);
	return rho_s_out;
}

double f_ex(double rho_s){
	double eta, f_ex_out;
	eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
	f_ex_out = k*T*eta*(4.0-3.0*eta)/std::pow((1.0-eta),2.0);
	return f_ex_out;
}

//Steele 10-4-3 potential
double phi_sf(double z){
	double phi_sf_out;
	phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*std::pow(sigma_sf,2.0)*delta*
				( (2.0/5.0)*std::pow((sigma_sf/z),10.0)-std::pow((sigma_sf/z),4.0)-std::pow(sigma_sf,4.0)/
				(3.0*delta*std::pow((0.61*delta+z),3.0)) );
	return phi_sf_out;
}

//e.g., Carbon slit
double phi_ext(double z){
	double phi_ext_out;
	phi_ext_out = phi_sf(z) + phi_sf(H-z);
	//std::cout << phi_ext_out << std::endl;
	return phi_ext_out;
}

double mu_ex(double rho_b){
	double y, mu_ex_out;
	y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	mu_ex_out = k*T*(8.0*y-9.0*y*y+3.0*y*y*y)/std::pow((1.0-y),3.0);
	return mu_ex_out;
}

double mu_b(double rho_b){
	double mu_id, mu_hs, mu_b_out;
	mu_id = k*T*std::log(std::pow(lam,3.0)*rho_b);
	mu_hs = mu_id + mu_ex(rho_b);
	mu_b_out = mu_hs - rho_b*alpha;
	return mu_b_out;
}

// d(f_ex)/d(rho_s)
double dfex_per_drhos(double rho_s){
	double dfex_per_drhos_out;
	double eta;
	eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
	dfex_per_drhos_out = k*T*(4.0-2.0*eta)/std::pow((1-eta),3.0)*M_PI*std::pow(d_hs,3.0)/6.0;
	return dfex_per_drhos_out;
}

// d(rho_s)/d(rho)
double drhos_per_drho(double *rho, double r1, double r2, double *r, double ra){
	double w, drhos_per_drho_out;
	// Percus-Yevick approximation, Tarazona theory
	//w = wi(std::abs(r1-r2),0) + wi(std::abs(r1-r2),1)*rho_s(rho,r1,r) + wi(std::abs(r1-r2),2)*std::pow(rho_s(rho,r1,r),2.0);
	w = wi(ra,0) + wi(ra,1)*rho_s(rho,r1,r) + wi(ra,2)*std::pow(rho_s(rho,r1,r),2.0);
	drhos_per_drho_out = w/(1.0-rho_si(rho,r2,r,1)-2.0*rho_si(rho,r2,r,2)*rho_s(rho,r2,r));
	return drhos_per_drho_out;
}

// xi include k*T*(std::log(rho_b)) type.
// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
double xi(double *rho, double r1, double rho_b, double *r){
	int j,k;
	double rho_dfex_int, rho_phi_int;
	double xi_out;
	double ra;
	rho_dfex_int = 0.0;
	rho_phi_int  = 0.0;
	double rho_dfex_int_j[nstep], rho_phi_int_j[nstep];
	double rho_dfex_int_k[nrmesh], rho_phi_int_k[nrmesh];
	for (j=0; j<nstep; j++) {
		for (k=0; k<nrmesh; k++) {
			ra = std::pow((r1-r[j]),2.0) + std::pow((double(k)*rc/double(nrmesh)),2.0);
			ra = std::pow(ra,0.5);
			//std::cout << ra << std::endl;
			//
			// d(f_ex)/d(rho) = d(f_ex)/d(rho_s) * d(rho_s)/d(rho)
			rho_dfex_int = rho_dfex_int + rho[j]*dfex_per_drhos(rho_s(rho,r[j],r))*drhos_per_drho(rho,r1,r[j],r,ra)*(2.0*M_PI*(double(k)*rc/double(nrmesh))*(rc/double(nrmesh)))*dr;
			rho_phi_int  = rho_phi_int  + rho[j]*phi_att(ra)*(2.0*M_PI*(double(k)*rc/double(nrmesh))*(rc/double(nrmesh)))*dr;
			//std::cout << rho_dfex_int << ", " << rho_phi_int << std::endl;
			//std::cout << dfex_per_drhos(rho_s(rho,r[j],r)) << ", " << drhos_per_drho(rho,r1,r[j],r) << std::endl;
			//
			rho_dfex_int_k[k] = rho[j]*dfex_per_drhos(rho_s(rho,r[j],r))*drhos_per_drho(rho,r1,r[j],r,ra)*(2.0*M_PI*(double(k)*rc/double(nrmesh)));
			rho_phi_int_k[k]  = rho[j]*phi_att(ra)*(2.0*M_PI*(double(k)*rc/double(nrmesh)));
		}
		//ingegral_simpson(double *f, int n, double dx)
		rho_dfex_int_j[j] = ingegral_simpson(rho_dfex_int_k, nrmesh, (rc/double(nrmesh)));
		rho_phi_int_j[j]  = ingegral_simpson(rho_phi_int_k, nrmesh, (rc/double(nrmesh)));
	}
	//ingegral_simpson(double *f, int n, double dx)
	rho_dfex_int = ingegral_simpson(rho_dfex_int_j, nstep, dr);
	rho_phi_int  = ingegral_simpson(rho_phi_int_j, nstep, dr);
	//
	//rho_dfex_int = rho_dfex_int / (M_PI*std::pow((rc),2.0)) / (nstep*dr);
	//rho_phi_int  = rho_phi_int  / (M_PI*std::pow((rc),2.0)) / (nstep*dr);
	//
	xi_out = k*T*std::log(rho_b) + mu_ex(rho_b) - rho_b*alpha - phi_ext(r1) - f_ex(rho_s(rho,r1,r)) - rho_dfex_int - rho_phi_int;
	//std::cout << "xi, (k*T)*log(rho_b), mu_ex(rho_b), -rho_b*alpha, -phi_ext(r1), -f_ex(rho_s(rho,r1,r)), -rho_dfex_int, -rho_phi_int" << std::endl;
	//std::cout << xi_out << ", " << k*T*std::log(rho_b) << ", " << mu_ex(rho_b) << ", " << -rho_b*alpha << ", " << -phi_ext(r1) << ", " << -f_ex(rho_s(rho,r1,r)) << ", " << -rho_dfex_int << ", " << -rho_phi_int << std::endl;
	return xi_out;
}

double press_hs(double rho_b){
	double y, press_hs_out;
	y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	press_hs_out = rho_b*k*T* (1.0 + y + y*y - y*y*y)/std::pow((1.0-y),3.0);
	return press_hs_out;
}

double Maxwell_construction(double *r){
	int i,j;
	int iter_max_drhob0 = 250000;
	int iter_max_dmue = 1500;
	double drhob0 = 0.0001;
	double dmue = 0.01;
	double threshold_diff = 0.3;
	double threshold_find = 0.3;
	//
	double mu_b_per_epsilon_ff[iter_max_drhob0];
	double mu_e_per_epsilon_ff;
	double diff,diffp;
	int flag;
	double rho_b0;
	//
	// rho_b vs. mu_b/epsilon_ff
	std::ofstream ofs("./Maxwell_construction_data.txt");
	ofs << "Chemical_potential(mu_b/epsilon_ff), Density(rho_b*d_hs^3)" << std::endl;
	for (i=0; i<iter_max_drhob0; i++){
		rho_b0 = drhob0*double(i+1.0);
		mu_b_per_epsilon_ff[i] = mu_b(rho_b0)/epsilon_ff;
		ofs << mu_b_per_epsilon_ff[i] << ", " << rho_b0*std::pow(d_hs,3.0) << std::endl;
		//std::cout << "rho_b0 = "<< rho_b0 << ", mu_b/epsilon_ff = " << mu_b_per_epsilon_ff[i] << std::endl;
	}
	// Maxwell equal area rule
	for (j=0; j<iter_max_dmue; j++){
		mu_e_per_epsilon_ff = dmue*double(j+1.0) - 12.0;
		diff = 0.0;
		flag = 0;
		for (i=0; i<iter_max_drhob0; i++){
			diffp = mu_b_per_epsilon_ff[i] - mu_e_per_epsilon_ff;
			if (diffp > 0.0 && flag != 2) {
				diff = diff + diffp*drhob0;
				flag = 1;
			} else if (diffp <= 0.0 && flag != 0) {
				diff = diff + diffp*drhob0;
				flag = 2;
			}
			//std::cout << diffp << std::endl;
		}
		//std::cout << "mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << ", diff = " << diff << std::endl;
		if (std::abs(diff) <= threshold_diff) {
			//std::cout << "mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << ", diff = " << diff << std::endl;
			break;
		}
	}
	// find rho_b0
	for (i=0; i<iter_max_drhob0; i++){
		rho_b0 = drhob0*double(i+1.0);
		//if ( std::abs(mu_b(rho_b0)/epsilon_ff - mu_e_per_epsilon_ff) <= threshold_find &&
		//	 0.05 <= rho_b0*std::pow(d_hs,3.0) &&  rho_b0*std::pow(d_hs,3.0) <= 0.75) {
		if ( std::abs(mu_b(rho_b0)/epsilon_ff - mu_e_per_epsilon_ff) <= threshold_find ) {
			//std::cout << "rho_b0 = " << rho_b0 << ", rho_b0*d_hs^3 = " << rho_b0*std::pow(d_hs,3.0) << std::endl;
			break;
		}
	}
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "Maxwell equal area rule" << std::endl;
	std::cout << "chemical potential, mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << std::endl;
	std::cout << "density, rho_b0*d_hs^3 = " << rho_b0*std::pow(d_hs,3.0) << std::endl;
	std::cout << "rho_b0 = " << rho_b0 << std::endl;
	return rho_b0;
}

int main(){
	int i,j,k;
	double w = 0.3;
	double diff;
	double v_gamma;
	double press_b, press_b0, pp0;
	double rho_b, rho_b0;
	//
	read_parameters();
	double r[nstep];
	double rho[nstep], rho_new[nstep];
	// set dr
	for (i=0; i<nstep; i++){
		//r[i] = sigma_ss/2.0 + (H-sigma_ss)/double(nstep)*double(i);
		r[i] = sigma_ss/2.0 + dr*double(i); // dr=(H-sigma_ss)/double(nstep)
		//std::cout << i << ", " << r[i] << std::endl;
	}
	// set rho_b0
	rho_b0 = Maxwell_construction(r);
	//std::cout << rho_b0 << std::endl;
	// initialization
	for (i=0; i<nstep; i++){
		rho[i] = rho_b0/(nstep*dr);
		rho_new[i] = 0.0;
	}
	// volume and pressure
	std::ofstream ofsppov("./PP0_vs_Vgamma_data.txt");
	ofsppov << "# w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov << "# P/P0, Vgamma" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	std::cout << "P/P0, Vgamma" << std::endl;
	for (k=0; k<100; k++){
		rho_b = rho_b0 * std::exp(-(20.0-2.0*double(k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		for (j=0; j<cycle_max; j++){
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dr*nstep/2. 
			for (i=0; i<=(nstep-1)/2; i++){
				//rho_new[i] = rho_b*std::exp(xi(rho,r[i],rho_b,r)/(k*T)); // this equation occure inf.
				rho_new[i] = std::exp(xi(rho,r[i],rho_b,r)/(k*T)); // xi include k*T*(std::log(rho_b)) type.
				//std::cout << "num of cycle i, r[i], rho_new[i], rho[i]" << std::endl;
				//std::cout << i << ", " << r[i] << ", "<< rho_new[i] << ", " << rho[i] << std::endl;
			}
			diff = 0.0;
			for (i=0; i<=(nstep-1)/2; i++){
				rho[i] = w*rho_new[i] + (1.0-w)*rho[i];
				rho[nstep-i] = rho[i]; // The rest is filled with mirror symmetry. 
				diff = diff + 2.0*std::abs((rho_new[i]-rho[i])/rho[i]);
			}
			if ( (diff/(nstep+1)*100.0) < 5.0 ) {
				break;
			}
			//std::cout << "--------------------------------------------------" << std::endl;
			//std::cout << "cycle=" << j << ", diff=" << diff << ", rho[nstep/2]=" << rho[nstep/2] << std::endl;
		}
		//
		v_gamma = 0.0;
		for (i=0; i<nstep/2; i++){
			//std::cout << r[i] << ", " << rho[i] << std::endl;
			v_gamma = v_gamma + 2.0*rho[i]*dr;
		}
		v_gamma = v_gamma/(H-sigma_ss) - rho_b;
		if (v_gamma < 0) { v_gamma = 0.0; }
		//std::cout << "V= " << v_gamma << std::endl;
		//
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		//std::cout << "P= " << press_b << std::endl;
		//std::cout << "P0= " << press_b0 << std::endl;
		pp0 = press_b/press_b0;
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov << pp0 << ", "<< v_gamma << std::endl;
		std::cout << pp0 << ", "<< v_gamma << std::endl;
	}
	return 0;
}
