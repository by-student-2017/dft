#include <fstream>   // for file in and out
#include <iostream>  // for cout
#include <cmath>     // for log, exp
#include <sstream>   // for read parameters
#include <omp.h>     // OpenMP (c++ nldft.cpp -fopenmp) (set OMP_NUM_THREADS=4)
#include <mpi.h>     // OpenMPI (mpic++ nldft.cpp) (set OMP_NUM_THREADS=1)

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

// compiling: mpic++ nldft_sda_slit_openmpi.cpp -O2
// usage: ./a.out

// debag mode
// compiling: c++ nldft_sda_slit_openmpi.cpp -g -Wall -O0
// run: gdb ./a.out
//      (gdb) run

// ---------- ----------- ------------ ------------
// Adsorbent 
//double H = 1.00; //distace of slit [nm]
//double sigma_ss = 0.34; // [nm]
//int nstep = 100;
//double w_pw = (H-sigma_ss); // pore width [nm]
//double dr = w_pw/double(nstep);
double H;
double sigma_ss;
//#define nstep=1001;
//constexpr int nstep = 1001;
int nstep;
double w_pw;
double dr;
// ---------- ----------- ------------ ------------
// assume rho is same value in x-y plane.
// cylinder and normalization, because of cut off (rc).
//int nrmesh = 20; //rho_si and xi function
int nrmesh;
//int ndmesh = d_hs*nrmesh/rc
int ndmesh;
//double drc = rc/double(nrmesh-1);
double drc;
//double dd = 2.0*d_hs/double(ndmesh-1);
double dd;
// ---------- ----------- ------------ ------------
// iteration of rho
//int cycle_max = 50;
int cycle_max;
//double wmixing = 0.005;
double wmixing;
// ---------- ----------- ------------ ------------
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
// ---------- ----------- ------------ ------------
//double rm = std::pow(2.0,1.0/6.0)*sigma_ff; //minimum position of LJ
//double rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
double rm;
// ---------- ----------- ------------ ------------
// Carbon dioxide/Carbon slit 81.5  [K](epsilon), 0.3430 [nm](sigma)
// Nitrogen/Carbon slit       53.72 [K](epsilon), 0.3508 [nm](sigma)
//double epsilon_sf = 53.72; // [K] 
//double sigma_sf = 0.3508; // [nm]
double epsilon_sf;
double sigma_sf;
// ---------- ----------- ------------ ------------
// slit pore (graphite)
//double delta = 0.335; // [nm]
//double rho_ss = 114.0; // [nm^-3], [molecules/nm3]?, 0.114 [A^-3]
double delta;
double rho_ss;
// ---------- ----------- ------------ ------------
//double m = 14.0067*2.0/(6.02214076e23)/1000; // N2 = 4.65173e-26 [kg]
//double m = 4.65173e-26; //[kg] (N2) (e.g., Ar = 6.63e-26 [kg])
double m;
double kb1 = 1.0;
double kb = 1.38e-23; //[J/K] (8.61733262e-5 [eV/K])
//extern double T = 77.347; //[K]
double T;
double h = 6.63e-34; //[Js] (4.135667696e-15 [eVs])
// thermal de Broglie wavelength
//extern double lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm], Maxwell_construction()
double lam;
// Ref: https://www1.doshisha.ac.jp/~bukka/lecture/statistic/pdftext/std-07.pdf
// ---------- ----------- ------------ ------------
// alpha = integal phi_att * -1.0
//extern double alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
//	( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
double alpha;
// ---------- ----------- ------------ ------------

//Barker-Henderson (BH) theory
double d_bh_calc(double epsilon, double sigma){
	//double epsilon = 94.45;
	//double sigma = 0.3575;
	//Lstoskie et al.,
	double xi1 = 0.3837;
	double xi2 = 1.035;
	double xi3 = 0.4249;
	double xi4 = 1.0;
	double d_bh_out;
	d_bh_out = (xi1*kb1*T/epsilon+xi2)/(xi3*kb1*T/epsilon+xi4)*sigma;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "d = d_hs = " << d_bh_out << " [nm] at " << T << " [K] from Barker-Henderson (BH) theory" << std::endl;
	return d_bh_out;
}

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
	// ---------- ----------- ------------ ------------
	H = num[0]; //distace of slit [nm]
	// ---------- ----------- ------------ ------------
	sigma_ss = num[1]; // [nm]
	// ---------- ----------- ------------ ------------
	nstep = int(num[2]);
	if ( nstep == 0 ) {
		nstep = int((H-sigma_ss)/0.02 + 0.5);
		if ( nstep%2 == 1 ){
			nstep = nstep + 1;
		}
		std::cout << "autoset nstep = " << nstep << std::endl;
	}
	// ---------- ----------- ------------ ------------
	cycle_max = int(num[3]);
	// ---------- ----------- ------------ ------------
	wmixing = num[4];
	// ---------- ----------- ------------ ------------
	epsilon_ff = num[5]; // [K]
	// ---------- ----------- ------------ ------------
	sigma_ff = num[6]; // [nm]
	// ---------- ----------- ------------ ------------
	d_hs = num[7]; // [nm]
	//if ( d_hs == 0.0 ) { d_hs = d_bh_calc(epsilon_ff, sigma_ff); }
	// move below (T)
	// ---------- ----------- ------------ ------------
	rc = num[8]; // [nm], cut off
	if ( rc == 0.0 ) { 
		rc = 5.0*sigma_ff;
		std::cout << "autoset (cut off) rc = " << rc << " [nm]" << std::endl;
	}
	// ---------- ----------- ------------ ------------
	nrmesh = int(num[9]);
	if ( nrmesh == 0 ) {
		nrmesh = int(rc/0.08 + 0.5);
		if ( nrmesh%2 == 1 ){
			nrmesh = nrmesh + 1;
		}
		std::cout << "autoset nrmesh = " << nrmesh << std::endl;
	}
	// ---------- ----------- ------------ ------------
	epsilon_sf = num[10]; // [K]
	// ---------- ----------- ------------ ------------
	sigma_sf = num[11]; // [nm]
	// ---------- ----------- ------------ ------------
	delta = num[12]; // nm
	// ---------- ----------- ------------ ------------
	rho_ss = num[13]; // [nm^-3], [mulecules/nm3]
	// ---------- ----------- ------------ ------------
	m = num[14]; // [kg]
	// ---------- ----------- ------------ ------------
	T = num[15]; // [K]
	if ( d_hs == 0.0 ) { d_hs = d_bh_calc(epsilon_ff, sigma_ff); }
	// ---------- ----------- ------------ ------------
	
	w_pw = (H-sigma_ss); // pore width [nm]
	dr = (H-sigma_ss)/double(nstep-1);
	rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
	
	// ---------- ----------- ------------ ------------
	
	
	//ndmesh = int(2*d_hs*nrmesh/rc); // why ? this setting occures nan.
	//if ( ndmesh < 9 ) { 
	//	ndmesh = 9;
	//	std::cout << "autoset ndmesh = " << ndmesh << std::endl;
	//}
	//if ( ndmesh%2 == 1 ) { ndmesh = ndmesh + 1; }
	//dd = 2.0*d_hs/double(ndmesh-1); // rho_si(), xi()
	ndmesh = nrmesh;
	drc = rc/double(nrmesh-1); // xi()
	dd = drc;
	
	// ---------- ----------- ------------ ------------
	
	// thermal de Broglie wavelength
	//lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm], Maxwell_construction()
	lam = h/std::sqrt(2.0*M_PI*m*kb*T)*1e9; //[nm], Maxwell_construction()
	
	// alpha = integal phi_att * -1.0
	alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
		( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
	// rm = rm when the potential is split according to the WCA schem and rm = simga_ff when the LJ potential is split according to the BH decomposition.
	
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "thermal de Broglie wavelength = " << lam << " [nm]" << std::endl;
	std::cout << "integal phi_att * -1.0 = alpha = " << alpha << std::endl;
}

double integral_simpson(double *f, int n, double dx){
	//if( (n+1)%2 == 1 ){
	//	std::cout << "Error, plase change number of data to even ( = array[odd] )" << std::endl;
	//}
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

double phi_att(double r){
	double e;
	// WCA (Weeks-Chandler-Anderson) type
	if (r < rm){
		e = - epsilon_ff;
	}else if (rm <= r && r <= rc){
		// Lennard-Jones（LJ) potential
		//e = 4.0*epsilon_ff*( std::pow((sigma_ff/r),12.0) - std::pow((sigma_ff/r),6.0) );
		e = std::pow((sigma_ff/r),6.0);
		e = 4.0*epsilon_ff*( e*e - e );
	}else {
		e = 0.0;
	}
	//std::cout << e << std::endl;
	return e;
}

// Percus-Yevick (PY) two-particle direct correlation function of the homogeneous hard-sphere fluid
double wi(double r, int i){
	double wi_out;
	double rpdhs;
	switch(i){
		case 0:
			if (r <= d_hs){
				//wi_out = 3.0/(4.0*M_PI*std::pow(d_hs,3.0));
				wi_out = 3.0/(4.0*M_PI*(d_hs*d_hs*d_hs));
			} else {
				wi_out = 0.0;
			}
			break;
		case 1:
			rpdhs = r/d_hs;
			if (r <= d_hs) {
				//wi_out = 0.475-0.648*(r/d_hs)+0.113*std::pow((r/d_hs),2.0);
				//wi_out = 0.475-0.648*(r/d_hs)+0.113*((r/d_hs)*(r/d_hs));
				wi_out = 0.475-0.648*(rpdhs)+0.113*(rpdhs*rpdhs);
			} else if (d_hs < r && r <= 2.0*d_hs) {
				//wi_out = 0.288*(d_hs/r)-0.924+0.764*(r/d_hs)-0.187*std::pow((r/d_hs),2.0);
				//wi_out = 0.288*(d_hs/r)-0.924+0.764*(r/d_hs)-0.187*((r/d_hs)*(r/d_hs));
				wi_out = 0.288/rpdhs-0.924+0.764*(rpdhs)-0.187*(rpdhs*rpdhs);
			} else {
				wi_out = 0.0;
			}
			break;
		case 2:
			rpdhs = r/d_hs;
			if (r <= d_hs) {
				//wi_out = 5.0*M_PI*std::pow(d_hs,3.0)/144.0 * (6.0-12.0*(r/d_hs)+5.0*std::pow((r/d_hs),2.0));
				//wi_out = 5.0*M_PI*(d_hs*d_hs*d_hs)/144.0 * (6.0-12.0*(r/d_hs)+5.0*((r/d_hs)*(r/d_hs)));
				wi_out = 5.0*M_PI*(d_hs*d_hs*d_hs)/144.0 * (6.0-12.0*(rpdhs)+5.0*(rpdhs*rpdhs));
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

// Tarazona theory
double rho_si(double *rho, double r1, double *r, int i){
	int j,k;
	double ra;
	double raj;
	double rak;
	//double ndmesh = 2*d_hs*nrmesh/rc;
	//double dd = 2.0*d_hs/double(ndmesh-1);
	//dd = drc;
	double tpidd = 2.0*M_PI*dd;
	double rho_si_out;
	double rho_si_int_j[nstep];
	double rho_si_int_k[nrmesh];
	rho_si_int_k[0] = 0.0;
// #pragma omp parallel for  private(k) // Slow down (I do not recommend)
	for (j=0; j<nstep; j++) {
		raj = (r1-r[j]);
		for (k=1; k<ndmesh; k++) {
			rak = dd*double(k);
			//ra = std::pow((r1-r[j]),2.0) + std::pow((double(k)*dd),2.0);
			//ra = (r1-r[j])*(r1-r[j]) + (double(k)*dd)*(double(k)*dd);
			ra = raj*raj + rak*rak;
			//ra = std::pow(ra,0.5);
			ra = std::sqrt(ra);
			//std::cout << ra << std::endl;
			//
			//rho_si_int_k[k] = rho[j]*wi(ra,i)*(2.0*M_PI*(double(k)*dd)); // old ver.1.1.0
			//rho_si_int_k[k] = wi(ra,i)*(2.0*M_PI*(double(k)*dd));
			rho_si_int_k[k] = wi(ra,i)*(tpidd*double(k));
		}
		//integral_simpson(double *f, int n, double dx)
		//rho_si_int_j[j] = integral_simpson(rho_si_int_k, ndmesh-1, dd); // old ver.1.1.0
		rho_si_int_j[j] = rho[j]*integral_simpson(rho_si_int_k, ndmesh-1, dd);
	}
	//integral_simpson(double *f, int n, double dx)
	rho_si_out = integral_simpson(rho_si_int_j, nstep-1, dr);
	//
	return rho_si_out;
}

// smoothed density approximation (SDA)
//double rho_s(double *rho, double r1, double *r){
//	double rho_den1, rho_den2, rho_s_out;
//	//rho_den1 = std::pow((1.0 - rho_si(rho,r1,r,1)),2.0);
//	rho_den1 = (1.0 - rho_si(rho,r1,r,1));
//	rho_den1 = rho_den1 * rho_den1;
//	//rho_den2 = std::pow((rho_den1 - 4.0*rho_si(rho,r1,r,0)*rho_si(rho,r1,r,2)),0.5);
//	rho_den2 = std::sqrt(rho_den1 - 4.0*rho_si(rho,r1,r,0)*rho_si(rho,r1,r,2));
//	rho_s_out = 2.0*rho_si(rho,r1,r,0)/(1.0 - rho_si(rho,r1,r,1)+rho_den2);
//	return rho_s_out;
//}

// smoothed density approximation (SDA), modified version
double rho_s(double *rho, double *r, double *rho_sj, double *rho_s0j, double *rho_s1j, double *rho_s2j){
	int j;
	double rho_den1j, rho_den2j;
//#pragma omp parallel for // An abnormal value occurred in OpenMP setting
	for (j=0; j<nstep; j++) {
		rho_s0j[j] = rho_si(rho, r[j], r, 0);
		rho_s1j[j] = rho_si(rho, r[j], r, 1);
		rho_s2j[j] = rho_si(rho, r[j], r, 2);
		//rho_den1j = std::pow((1.0 - rho_s1j[j]),2.0);
		rho_den1j = (1.0 - rho_s1j[j]);
		rho_den1j = rho_den1j * rho_den1j;
		//rho_den2j = std::pow((rho_den1j - 4.0*rho_s0j[j]*rho_s2j[j]),0.5);
		//rho_den2j = std::sqrt(rho_den1j - 4.0*rho_s0j[j]*rho_s2j[j]);
		rho_den2j = rho_den1j - 4.0*rho_s0j[j]*rho_s2j[j];
		// to avoide nan
		if ( rho_den2j > 0 ) {
			rho_den2j = std::sqrt(rho_den2j);
		} else {
			rho_den2j = 0.0;
		}
		rho_sj[j] = 2.0*rho_s0j[j]/(1.0 - rho_s1j[j]+rho_den2j);
		//std::cout << j << ", " << rho[j] << ", " << rho_sj[j] << ", " << rho_s0j[j] << ", " << rho_s1j[j] << ", " << rho_s2j[j] << std::endl;
		//std::cout << rho_den1j << ", " << rho_den2j << std::endl;
	}
	return 0;
}

// Steele 10-4-3 potential
double phi_sf(double z){
	double phi_sf_out;
	double sigma_sf2 = sigma_sf*sigma_sf;
	double sfpz = (sigma_sf/z);
	double sfpz2 = sfpz*sfpz;
	double dez = (0.61*delta+z);
	//phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*std::pow(sigma_sf,2.0)*delta*
	//			( (2.0/5.0)*std::pow((sigma_sf/z),10.0)-std::pow((sigma_sf/z),4.0)-std::pow(sigma_sf,4.0)/
	//			(3.0*delta*std::pow((0.61*delta+z),3.0)) );
	phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*(sigma_sf2)*delta*
				( (2.0/5.0)*std::pow(sfpz2,5.0)-(sfpz2*sfpz2)-(sigma_sf2*sigma_sf2)/
				(3.0*delta*(dez*dez*dez)) );
	return phi_sf_out;
}

// e.g., wall potential (Carbon slit)
double phi_ext(double z){
	double phi_ext_out;
	phi_ext_out = phi_sf(z) + phi_sf(H-z);
	//std::cout << phi_ext_out << std::endl;
	return phi_ext_out;
}

// from Carnahan-Starling (CS) equation of state
double mu_ex(double rho_b){
	double y, mu_ex_out;
	//y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
	double den1y = (1.0-y);
	//mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/std::pow((1.0-y),3.0);
	//mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/((1.0-y)*(1.0-y)*(1.0-y));
	mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/(den1y*den1y*den1y);
	return mu_ex_out;
}

double mu_b(double rho_b){
	double mu_id, mu_hs, mu_b_out;
	//mu_id = kb1*T*std::log(std::pow(lam,3.0)*rho_b);
	mu_id = kb1*T*std::log((lam*lam*lam)*rho_b);
	mu_hs = mu_id + mu_ex(rho_b);
	mu_b_out = mu_hs - rho_b*alpha;
	return mu_b_out;
}

double f_ex(double rho_s){
	double eta, f_ex_out;
	//eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
	eta = M_PI*rho_s*(d_hs*d_hs*d_hs)/6.0;
	double den1e = (1.0-eta);
	//f_ex_out = kb1*T*eta*(4.0-3.0*eta)/std::pow((1.0-eta),2.0);
	//f_ex_out = kb1*T*eta*(4.0-3.0*eta)/((1.0-eta)*(1.0-eta));
	f_ex_out = kb1*T*eta*(4.0-3.0*eta)/(den1e*den1e);
	return f_ex_out;
}

// d(f_ex)/d(rho_s)
double dfex_per_drhos(double rho_s){
	double dfex_per_drhos_out;
	double eta;
	//eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
	eta = M_PI*rho_s*(d_hs*d_hs*d_hs)/6.0;
	double den1e = (1.0-eta);
	//dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/std::pow((1.0-eta),3.0)*M_PI*std::pow(d_hs,3.0)/6.0;
	//dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/((1.0-eta)*(1.0-eta)*(1.0-eta))*M_PI*(d_hs*d_hs*d_hs)/6.0;
	dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/(den1e*den1e*den1e)*(M_PI*(d_hs*d_hs*d_hs)/6.0);
	return dfex_per_drhos_out;
}

// d(rho_s)/d(rho)
//double drhos_per_drho(double *rho, double r1, double r2, double *r, double ra){
//	double w, drhos_per_drho_out;
//	// Percus-Yevick approximation, Tarazona theory
//	w = wi(ra,0) + wi(ra,1)*rho_s(rho,r1,r) + wi(ra,2)*std::pow(rho_s(rho,r1,r),2.0);
//	drhos_per_drho_out = w/(1.0-rho_si(rho,r2,r,1)-2.0*rho_si(rho,r2,r,2)*rho_s(rho,r2,r));
//	return drhos_per_drho_out;
//}

// d(rho_s)/d(rho), modified version
double drhos_per_drho_j(double ra, double rho_sj, double rho_s1j, double rho_s2j){
	double w, drhos_per_drho_out;
	// Percus-Yevick approximation, Tarazona theory
	//w = wi(ra,0) + wi(ra,1)*rho_sj + wi(ra,2)*std::pow(rho_sj,2.0);
	w = wi(ra,0) + wi(ra,1)*rho_sj + wi(ra,2)*(rho_sj*rho_sj);
	drhos_per_drho_out = w/(1.0-rho_s1j-2.0*rho_s2j*rho_sj);
	return drhos_per_drho_out;
}

double calc_alpha(double *r){
	int i,j,k;
	double ra;
	double raj;
	double rak;
	//double drc = rc/double(nrmesh-1);
	double tpidrc = 2.0*M_PI*drc;
	double alpha_other_method;
	double alpha_int_j[nstep];
	double alpha_int_k[nrmesh];
	alpha_int_k[0] = 0.0;
//#pragma omp parallel for private(k)
	for (i=0; i<=(nstep-2)/2; i++){
		for (j=0; j<nstep; j++) {
			raj = (r[i]-r[j]);
			for (k=1; k<nrmesh; k++) {
				rak = drc*double(k);
				//ra = std::pow((r[i]-r[j]),2.0) + std::pow((double(k)*drc),2.0);
				//ra = (r[i]-r[j])*(r[i]-r[j]) + (double(k)*drc)*(double(k)*drc);
				ra = raj*raj + rak*rak;
				//ra = std::pow(ra,0.5);
				ra = std::sqrt(ra);
				//std::cout << ra << std::endl;
				//alpha_int_k[k]  = -phi_att(ra)*(2.0*M_PI*(double(k)*drc));
				alpha_int_k[k]  = -phi_att(ra)*(tpidrc*double(k));
			}
			//integral_simpson(double *f, int n, double dx)
			alpha_int_j[j]  = integral_simpson(alpha_int_k, nrmesh-1, drc);
		}
		//integral_simpson(double *f, int n, double dx)
		//alpha_other_method  = alpha_other_method + integral_simpson(alpha_int_j, nstep-1, dr)*2.0*dr;
		alpha_other_method  = alpha_other_method + integral_simpson(alpha_int_j, nstep-1, dr);
	}
	alpha_other_method  = alpha_other_method * 2.0 * dr / (H-sigma_ss);
	//std::cout << "--------------------------------------------------" << std::endl;
	//std::cout << "average alpha of other method = " << alpha_other_method << " in (carbon) slit" << std::endl;
	return alpha_other_method;
}

double phi_att_int(double *r, double *phi_att_int_ij){
	int i,j,k;
	double ra;
	double raj;
	double rak;
	//double drc = rc/double(nrmesh-1);
	//dd = drc;
	double phi_int_k[nrmesh];
	double tpidrc = 2.0*M_PI*drc;
	phi_int_k[0] = 0.0;
//#pragma omp parallel for private(k) // -nan is occurred
	for (i=0; i<nstep; i++) {
		for (j=0; j<nstep; j++) {
			raj = (r[i]-r[j]);
			for (k=1; k<nrmesh; k++) {
				rak = drc*double(k);
				//ra = std::pow((r[i]-r[j]),2.0) + std::pow((double(k)*drc),2.0);
				//ra = (r[i]-r[j])*(r[i]-r[j]) + (double(k)*drc)*(double(k)*drc);
				ra = raj*raj + rak*rak;
				//ra = std::pow(ra,0.5);
				ra = std::sqrt(ra);
				//std::cout << ra << std::endl;
				//rho_phi_int_k[k]  = rho[j]*phi_att(ra)*(2.0*M_PI*(double(k)*drc)); // old ver.1.1.0
				//phi_int_k[k]  = phi_att(ra)*(2.0*M_PI*(double(k)*drc));
				phi_int_k[k]  = phi_att(ra)*(tpidrc*double(k));
			}
			phi_att_int_ij[i*nstep+j] = integral_simpson(phi_int_k, nrmesh-1, drc);
		}
	}
	return 0;
}

// xi include kb1*T*(std::log(rho_b)) type.
// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
double xi(double *rho, double *r, int i, double rho_b, double *rho_sj, double *rho_s0j, double *rho_s1j, double *rho_s2j, double *phi_att_int_ij){
	int j,k;
	double ra;
	double raj;
	double rak;
	//double ndmesh = 2*d_hs*nrmesh/rc;
	//double dd = 2.0*d_hs/double(ndmesh-1);
	//double drc = rc/double(nrmesh-1);
	//dd = drc;
	double tpidd = 2.0*M_PI*dd;
	double rho_dfex_int_j[nstep];
	double rho_phi_int_j[nstep];
	double rho_dfex_int_k[nrmesh];
	double rho_phi_int_k[nrmesh]; // old ver.1.1.1
	rho_phi_int_k[0] = 0.0;
#pragma omp parallel for private(k)
	for (j=0; j<nstep; j++) {
		raj = (r[i]-r[j]);
		for (k=1; k<ndmesh; k++) {
			rak = dd*double(k);
			//ra = std::pow((r[i]-r[j]),2.0) + std::pow((double(k)*dd),2.0);
			//ra = (r[i]-r[j])*(r[i]-r[j]) + (double(k)*dd)*(double(k)*dd);
			ra = raj*raj + rak*rak;
			//ra = std::pow(ra,0.5);
			ra = std::sqrt(ra);
			//std::cout << ra << std::endl;
			//
			// d(f_ex)/d(rho) = d(f_ex)/d(rho_s) * d(rho_s)/d(rho)
			//rho_dfex_int_k[k] = rho[j]*dfex_per_drhos(rho_sj[j])*drhos_per_drho_j(ra, rho_sj[j], rho_s1j[j], rho_s2j[j])*(2.0*M_PI*(double(k)*dd)); // old ver.1.1.0
			//rho_dfex_int_k[k] = drhos_per_drho_j(ra, rho_sj[j], rho_s1j[j], rho_s2j[j])*(2.0*M_PI*(double(k)*dd));
			rho_dfex_int_k[k] = drhos_per_drho_j(ra, rho_sj[j], rho_s1j[j], rho_s2j[j])*(tpidd*double(k));
		}
		//integral_simpson(double *f, int n, double dx)
		//rho_dfex_int_j[j] = integral_simpson(rho_dfex_int_k, ndmesh-1, dd); // old ver.1.1.1
		rho_dfex_int_j[j] = rho[j]*dfex_per_drhos(rho_sj[j])*integral_simpson(rho_dfex_int_k, ndmesh-1, dd);
		//
		//for (k=1; k<nrmesh; k++) { // old ver.1.1.1
		//	//ra = std::pow((r[i]-r[j]),2.0) + std::pow((double(k)*drc),2.0);
		//	ra = (r[i]-r[j])*(r[i]-r[j]) + (double(k)*drc)*(double(k)*drc);
		//	//ra = std::pow(ra,0.5);
		//	ra = std::sqrt(ra);
		//	//std::cout << ra << std::endl;
		//	//rho_phi_int_k[k]  = rho[j]*phi_att(ra)*(2.0*M_PI*(double(k)*drc)); // old ver.1.1.0
		//	rho_phi_int_k[k]  = phi_att(ra)*(2.0*M_PI*(double(k)*drc));
		//} // old ver.1.1.1
		//integral_simpson(double *f, int n, double dx)
		//rho_phi_int_j[j]  = integral_simpson(rho_phi_int_k, nrmesh-1, drc); // old ver.1.1.0
		//rho_phi_int_j[j]  = rho[j]*integral_simpson(rho_phi_int_k, nrmesh-1, drc); // old ver.1.1.1
		rho_phi_int_j[j]  = rho[j]*phi_att_int_ij[i*nstep+j];
	}
	double rho_dfex_int;
	double rho_phi_int;
	//integral_simpson(double *f, int n, double dx)
	rho_dfex_int = integral_simpson(rho_dfex_int_j, nstep-1, dr);
	rho_phi_int  = integral_simpson(rho_phi_int_j, nstep-1, dr);
	//
	double xi_out;
	//xi_out = kb1*T*std::log(rho_b) + mu_ex(rho_b) - rho_b*alpha - phi_ext(r[i]) - f_ex(rho_sj[i]) - rho_dfex_int - rho_phi_int; // old ver.1.1.1
	xi_out = ( - rho_b*alpha - rho_dfex_int - f_ex(rho_sj[i]) ) + ( mu_ex(rho_b) - rho_phi_int ) + ( kb1*T*std::log(rho_b) - phi_ext(r[i]) );
	// debug
	//std::cout << "xi, (kb1*T)*log(rho_b), mu_ex(rho_b), -rho_b*alpha, -phi_ext(r[i]), -f_ex(rho_s(rho,r[i],r)), -rho_dfex_int, -rho_phi_int" << std::endl;
	//std::cout << xi_out << ", " << kb1*T*std::log(rho_b) << ", " << mu_ex(rho_b) << ", " << -rho_b*alpha << ", " << -phi_ext(r[i]) << ", " << -f_ex(rho_sj[i]) << ", " << -rho_dfex_int << ", " << -rho_phi_int << std::endl;
	//if ( std::isnan(rho_sj[i]) || std::isnan(f_ex(rho_sj[i])) || std::isnan(rho_dfex_int) || std::isnan(rho_phi_int) ){
	//	std::cout << i << ", " << rho[i] << ", " << rho_sj[i] << ", " << f_ex(rho_sj[i]) << ", " << rho_dfex_int << ", " << rho_phi_int << std::endl;
	//	std::exit(1);
	//}
	return xi_out;
}

double press_hs(double rho_b){
	double y, press_hs_out;
	//y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
	double den1y = (1.0-y);
	//press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/std::pow((1.0-y),3.0);
	//press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/((1.0-y)*(1.0-y)*(1.0-y));
	press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/(den1y*den1y*den1y);
	return press_hs_out;
}

double Maxwell_construction(void){
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
	double rho_b0_gas, rho_b0_metastable, rho_b0_liquid;
	double press_b0;
	//
	// rho_b vs. mu_b/epsilon_ff
	std::ofstream ofs("./Maxwell_construction_data.txt");
	ofs << "Chemical_potential(mu_b/epsilon_ff), Density(rho_b*d_hs^3)" << std::endl;
	for (i=0; i<iter_max_drhob0; i++){
		rho_b0 = drhob0*double(i+1.0);
		mu_b_per_epsilon_ff[i] = mu_b(rho_b0)/epsilon_ff;
		//ofs << mu_b_per_epsilon_ff[i] << ", " << rho_b0*std::pow(d_hs,3.0) << std::endl;
		ofs << mu_b_per_epsilon_ff[i] << ", " << rho_b0*(d_hs*d_hs*d_hs) << std::endl;
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
		rho_b0 = drhob0*double(j+1.0);
		if (std::abs(diff) <= threshold_diff) {
			//std::cout << "mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << ", diff = " << diff << std::endl;
			break;
		}
	}
	// find rho_b0
	flag = 0;
	for (i=0; i<iter_max_drhob0; i++){
		rho_b0 = drhob0*double(i+1.0);
		//if ( std::abs(mu_b(rho_b0)/epsilon_ff - mu_e_per_epsilon_ff) <= threshold_find &&
		//	 0.05 <= rho_b0*std::pow(d_hs,3.0) &&  rho_b0*std::pow(d_hs,3.0) <= 0.75) {
		if ( std::abs(mu_b(rho_b0)/epsilon_ff - mu_e_per_epsilon_ff) <= threshold_find ) {
			//std::cout << "rho_b0 = " << rho_b0 << ", rho_b0*d_hs^3 = " << rho_b0*std::pow(d_hs,3.0) << std::endl;
			if ( flag == 0 ){
				rho_b0_gas = rho_b0;
				flag = 1;
			} else if ( flag == 1 && 5.0*rho_b0_gas < rho_b0 ){
				rho_b0_metastable = rho_b0;
				flag = 2;
			} else if ( flag == 2 && 1.5*rho_b0_metastable < rho_b0 ){
				rho_b0_liquid = rho_b0;
				break;
			}
		}
	}	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "Maxwell construction (Maxwell equal area rule)" << std::endl;
	std::cout << "chemical potential, mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << std::endl;
	rho_b0 = rho_b0_gas;
	//std::cout << "density, rho_b0*d_hs^3 = " << rho_b0*std::pow(d_hs,3.0) << std::endl;
	std::cout << "density, rho_b0*d_hs^3 = " << rho_b0*(d_hs*d_hs*d_hs) << std::endl;
	//press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
	press_b0 = press_hs(rho_b0) - 0.5*(rho_b0*rho_b0)*alpha;
	std::cout << "Bulk pressure, P0 = " << press_b0 << " (rho_b0 = " << rho_b0 << ")" <<std::endl;
	std::cout << std::endl;
	std::cout << "gas phase   : rho_b0_gas        = " << rho_b0_gas        << ", rho_b0_gas*d_hs^3        = " << rho_b0_gas*std::pow(d_hs,3.0)        << std::endl;
	std::cout << "metastable  : rho_b0_metastable = " << rho_b0_metastable << ", rho_b0_metastable*d_hs^3 = " << rho_b0_metastable*std::pow(d_hs,3.0) << std::endl;
	std::cout << "liquid phase: rho_b0_liquid     = " << rho_b0_liquid     << ", rho_b0_liquid*d_hs^3     = " << rho_b0_liquid*std::pow(d_hs,3.0)     << std::endl;
	return rho_b0;
}

//int main(int argc, char **argv){
//MPI::Init(argc,argv);
int main(void){
MPI::Init();
	int i,j,k;
	double diff;
	double v_gamma;
	double press_b, press_b0, pp0;
	double rho_b, rho_b0;
	//
	read_parameters();
	double r[nstep];
	double rho[nstep], rho_new[nstep];
	//
#pragma omp parallel for
	for (i=0; i<nstep; i++){
		// r[i] = sigma_ss/2.0 + (H-sigma_ss)/double(nstep)*double(i);
		// 1.72 times is escape nan, etc from positive value of wall potential
		//r[i] = sigma_ss*1.74/2.0 + dr*double(i) + dr/2.0; // dr = (H-sigma_ss*1.74)/double(nstep+1);
		r[i] = sigma_ss/2.0 + dr*double(i); // dr = (H-sigma_ss)/double(nstep+1);
		//std::cout << i << ", " << r[i] << std::endl;
	}
	
	// show alpha
	//calc_alpha(r);
	// alpha = calc_alpha(r);
	
	// set rho_b0
	rho_b0 = Maxwell_construction();
	
	//std::cout << rho_b0 << std::endl;
	// initialization
#pragma omp parallel for
	for (i=0; i<nstep; i++){
		rho[i] = rho_b0/(nstep*dr);
		rho_new[i] = 0.0;
	}
	// volume and pressure
	std::ofstream ofsppov("./PP0_vs_Vgamma_data.txt");
	ofsppov << "# w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov << "# P/P0, V[molecules/nm^3]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	std::cout << "P/P0, V[molecules/nm^3]" << std::endl;
	double rho_sj[nstep];
	double rho_s0j[nstep];
	double rho_s1j[nstep];
	double rho_s2j[nstep];
	double phi_att_int_ij[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	phi_att_int(r, phi_att_int_ij); // calculate integral phi_att at r[i]
	for (k=0; k<100; k++){
		rho_b = rho_b0 * std::exp(-(20.0-2.0*double(k+1.0)/10.0));
		//rho_b = rho_b0 * std::exp(-(20.0-2.0*double(99.0-k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		for (j=0; j<cycle_max; j++){
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dr*nstep/2. 
			rho_s(rho, r, rho_sj, rho_s0j, rho_s1j, rho_s2j); 
//Since the function contains many variables specified by "private", the setting of "OpenMP" in this loop is not suitable.
//#pragma omp parallel for // An abnormal value occurred in OpenMP setting
			for (i=0; i<=(nstep-2)/2; i++){
				//rho_new[i] = rho_b*std::exp(xi(rho,r[i],rho_b,r)/(kb1*T)); // this equation occure inf.
				rho_new[i] = std::exp(xi(rho,r,i,rho_b, rho_sj, rho_s0j, rho_s1j, rho_s2j, phi_att_int_ij)/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				//std::cout << "num of cycle i, r[i], rho_new[i], rho[i]" << std::endl;
				//std::cout << i << ", " << r[i] << ", "<< rho_new[i] << ", " << rho[i] << std::endl;
				//std::cout << i << ", " << rho[i] << ", " << rho_sj[i] << ", " << rho_s0j[i] << ", " << rho_s1j[i] << ", " << rho_s2j[i] << std::endl;
			}
			diff = 0.0;
#pragma omp parallel for
			for (i=0; i<=(nstep-2)/2; i++){
				rho[i] = wmixing*rho_new[i] + (1.0-wmixing)*rho[i];
				rho[(nstep-1)-i] = rho[i]; // The rest is filled with mirror symmetry. 
				diff = diff + 2.0*std::abs((rho_new[i]-rho[i])/rho[i]);
			}
			if ( diff/nstep < 0.005) {
				break;
			}
			//std::cout << "--------------------------------------------------" << std::endl;
			//std::cout << "cycle=" << j << ", diff=" << diff << ", rho[nstep/2]=" << rho[nstep/2] << std::endl;
		}
		//for (i=0; i<nstep; i++){
		//	std::cout << "--------------------------------------------------" << std::endl;
		//	std::cout << "cycle=" << j << ", r[" << i << "]" << r[i] << " , -ext " << - phi_ext(r[i]) << ", rho[" << i << "]=" << rho[i] << std::endl;
		//}
		//
		//v_gamma = 0.0;
		//for (i=0; i<=(nstep-2)/2; i++){
			//std::cout << i << ", " << r[i] << ", " << rho[i] << std::endl;
			//v_gamma = v_gamma + 2.0*rho[i]*dr;
		//}
		v_gamma = integral_simpson(rho, nstep-1, dr);
		v_gamma = v_gamma/(H-sigma_ss) - rho_b;
		if (v_gamma < 0) { v_gamma = 0.0; }
		//v_gamma = v_gamma * (0.8064/28.0134/1e21*6.02214e23)/rho_b;
		// N2(77K): 0.8064 g/mL, 0.8064/28.0134 mol/mL, 0.8064/28.0134/1e21 mol/nm3, 0.8064/28.0134/1e21*6.02214e23 molecules/nm3
		//std::cout << "V= " << v_gamma << std::endl;
		// press_hs(rho_b) from Carnahan-Starling (CS) equation of state
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		//std::cout << "P= " << press_b << std::endl;
		//std::cout << "P0= " << press_b0 << std::endl;
		pp0 = press_b/press_b0;
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov << pp0 << ", "<< v_gamma << std::endl;
		std::cout << pp0 << ", "<< v_gamma << std::endl;
	}
	// reverse
	for (k=0; k<100; k++){
		//rho_b = rho_b0 * std::exp(-(20.0-2.0*double(k+1.0)/10.0));
		rho_b = rho_b0 * std::exp(-(20.0-2.0*double(99.0-k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		for (j=0; j<cycle_max; j++){
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dr*nstep/2. 
			rho_s(rho, r, rho_sj, rho_s0j, rho_s1j, rho_s2j);
//Since the function contains many variables specified by "private", the setting of "OpenMP" in this loop is not suitable.
//#pragma omp parallel for // An abnormal value occurred in OpenMP setting
			for (i=0; i<=(nstep-2)/2; i++){
				//rho_new[i] = rho_b*std::exp(xi(rho,r[i],rho_b,r)/(kb1*T)); // this equation occure inf.
				rho_new[i] = std::exp(xi(rho,r,i,rho_b, rho_sj, rho_s0j, rho_s1j, rho_s2j, phi_att_int_ij)/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				//std::cout << "num of cycle i, r[i], rho_new[i], rho[i]" << std::endl;
				//std::cout << i << ", " << r[i] << ", "<< rho_new[i] << ", " << rho[i] << std::endl;
			}
			diff = 0.0;
#pragma omp parallel for
			for (i=0; i<=(nstep-2)/2; i++){
				rho[i] = wmixing*rho_new[i] + (1.0-wmixing)*rho[i];
				rho[(nstep-1)-i] = rho[i]; // The rest is filled with mirror symmetry. 
				diff = diff + 2.0*std::abs((rho_new[i]-rho[i])/rho[i]);
			}
			if ( diff/nstep < 0.005) {
				break;
			}
			//std::cout << "--------------------------------------------------" << std::endl;
			//std::cout << "cycle=" << j << ", diff=" << diff << ", rho[nstep/2]=" << rho[nstep/2] << std::endl;
		}
		//for (i=0; i<nstep; i++){
		//	std::cout << "--------------------------------------------------" << std::endl;
		//	std::cout << "cycle=" << j << ", r[" << i << "]" << r[i] << " , -ext " << - phi_ext(r[i]) << ", rho[" << i << "]=" << rho[i] << std::endl;
		//}
		//
		//v_gamma = 0.0;
		//for (i=0; i<=(nstep-2)/2; i++){
			//std::cout << r[i] << ", " << rho[i] << std::endl;
		//	v_gamma = v_gamma + 2.0*rho[i]*dr;
		//}
		v_gamma = integral_simpson(rho, nstep-1, dr);
		v_gamma = v_gamma/(H-sigma_ss) - rho_b;
		if (v_gamma < 0) { v_gamma = 0.0; }
		//v_gamma = v_gamma * (0.8064/28.0134/1e21*6.02214e23)/rho_b;
		// N2(77K): 0.8064 g/mL, 0.8064/28.0134 mol/mL, 0.8064/28.0134/1e21 mol/nm3, 0.8064/28.0134/1e21*6.02214e23 molecules/nm3
		//std::cout << "V= " << v_gamma << std::endl;
		// press_hs(rho_b) from Carnahan-Starling (CS) equation of state
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		//std::cout << "P= " << press_b << std::endl;
		//std::cout << "P0= " << press_b0 << std::endl;
		pp0 = press_b/press_b0;
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov << pp0 << ", "<< v_gamma << std::endl;
		std::cout << pp0 << ", "<< v_gamma << std::endl;
	}
MPI::Finalize();
	return 0;
}
