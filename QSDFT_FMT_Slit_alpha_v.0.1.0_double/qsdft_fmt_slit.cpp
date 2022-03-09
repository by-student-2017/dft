#include <fstream>   // for file in and out
#include <iostream>  // for cout
#include <cmath>     // for log, exp
#include <sstream>   // for read parameters

//#include "maxwell_construction.h"

using namespace std;

// Fundamental Measure Theory (FMT)
// Quenched Solid Density Functional Theory (QSDFT)

//This QSDFT derived from NLDFT.
//  Smoothed Density Approximation (SDA)
//  Non-Local Density Functional Theory（NLDFT)
//  Reference: https://www.j-ad.org/adsorption_news/30_1.pdf

// Note
// This routine assumes that rho, etc is same value in x-y plane.
// Because of cut off (rc), it calculate circle in x-y plane.
// Units are fundamentaly [K] and [nm] in this routine.

// There are many imperfections, so I hope someone can make it better with a CC0 license. 
// It seems that this code is the first in the world at present (2021/7/5) to be released on CC0 even in NLDFT. 

// compiling: c++ qsdft_fmt_slit.cpp -O2
// usage: ./a.out

// debag mode
// compiling: c++ qsdft_fmt_slit.cpp -g -Wall -O0
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
// fluid-fluid
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
double rc;  // for fluid
double rcsf; // for solid-fluid
// ---------- ----------- ------------ ------------
//double rm = std::pow(2.0,1.0/6.0)*sigma_ff; //minimum position of LJ
//double rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
double rm;   // fluid
double rmsf; // solid-fluid
// ---------- ----------- ------------ ------------
// fluid-solid
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
double kb1 = 1.0;  // use all, except thermal de Broglie wavelength calculation.
// ---------- ----------- ------------ ------------
// thermal de Broglie wavelength
//extern double lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm], Maxwell_construction()
double kb = 1.38e-23; //[J/K] (8.61733262e-5 [eV/K])
// ----------
//double m = 14.0067*2.0/(6.02214076e23)/1000; // N2 = 4.65173e-26 [kg]
//double m = 4.65173e-26; //[kg] (N2) (e.g., Ar = 6.63e-26 [kg])
double m;  // for fluid
double ms; // for solid
//extern double T = 77.347; //[K]
// ----------
double T;
// ----------
double h = 6.63e-34; //[Js] (4.135667696e-15 [eVs])
// ----------
double lam;  // for fluid
double lams; // for solid
// Ref: https://www1.doshisha.ac.jp/~bukka/lecture/statistic/pdftext/std-07.pdf
// ---------- ----------- ------------ ------------
// alpha = integal phi_att_ff * -1.0
//extern double alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
//	( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
double alpha;
// ---------- ----------- ------------ ------------
// rho_b0 is related with P0
double rho_b0;
// ---------- ----------- ------------ ------------
// The edge position ze of the solid wall
double h0;
double ze;
// ---------- ----------- ------------ ------------
// Ris [nm] is the hard-sphere radius of solid (for QSDFT)
double Ris;
// ---------- ----------- ------------ ------------

double integral_trapezoidal(double *f, int n, double dx){
	double sum;
	sum = 0.0;
	int i;
	for(i=1; i<n; i++){
		sum += (f[i-1]+f[i])/2.0*dx;
	}
	return sum;
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

//Barker-Henderson (BH) perturbation theory
//double d_bh_calc(double epsilon, double sigma){
//	//double epsilon = 94.45;
//	//double sigma = 0.3575;
//	//Lstoskie et al.,
//	double d_bh_out;
//	double Ts = kb1*T/epsilon;
//	d_bh_out = (1.0+0.2977*Ts)/(1.0+0.331637*Ts+0.00104771*Ts*Ts)*sigma;
//	std::cout << "--------------------------------------------------" << std::endl;
//	std::cout << "d = d_hs = " << d_bh_out << " [nm] at " << T << " [K] from Barker-Henderson (BH) perturbation theory" << std::endl;
//	return d_bh_out;
//}

double rho_ssq(double z){
	double rho_ssq_out;
	//double h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	//double rho_ss = 114.0; // [molecules/nm3] (the density of bulk carbon)
	//double delta = 0.13; // [nm] (the roughness parameter) (the half-width of the density ramp)
	if ( 0.0 <= z && z < h0 ){
		rho_ssq_out = rho_ss;
	} else if ( h0 <= z && z < h0+2.0*delta ){
		rho_ssq_out = 0.75*rho_ss * (1.0 - (z - h0)/(2.0*delta));
	} else {
		rho_ssq_out = 0.0;
	}
	return rho_ssq_out;
}

// The edge position ze of the solid wall
double calc_ze(int ze_nstep){
	int i;
	double ze_out;
	//double h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	//double rho_ss = 114.0; // [molecules/nm3] (the density of bulk carbon)
	//double delta = 0.13; // [nm] (the roughness parameter) (the half-width of the density ramp)
	double rho_ssi[ze_nstep];
	double dss = (2.0*delta)/(ze_nstep-1);
	for (i=0; i<ze_nstep; i++){
		rho_ssi[i] = rho_ssq(h0+double(i)*dss);
	}
	//integral_trapezoidal(double *f, int n, double dx)
	//ze_out = integral_trapezoidal(rho_ssi, ze_nstep-1, dss)/rho_ss + h0;
	//integral_simpson(double *f, int n, double dx)
	ze_out = integral_simpson(rho_ssi, ze_nstep-1, dss)/rho_ss + h0;
	return ze_out;
}

void read_parameters(void){
	std::ifstream ifs("parameters.txt");
	std::string str;
	double num[25];
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
	std::cout << "--------------------------------------------------" << std::endl;
	//
	// ---------- ----------- ------------ ------------
	H = num[0]; //distace of slit [nm]
	// ---------- ----------- ------------ ------------
	sigma_ss = num[1]; // [nm]
	// ---------- ----------- ------------ ------------
	nstep = int(num[2]);
	// move below (ze)
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
	// move below(sigma_sf)
	// ---------- ----------- ------------ ------------
	nrmesh = int(num[9]);
	if ( nrmesh == 0 ) {
		nrmesh = int(rc/0.08 + 0.5);
		if ( nrmesh%2 == 0 ){
			nrmesh = nrmesh + 1;
		}
		std::cout << "autoset nrmesh = " << nrmesh << std::endl;
	}
	// ---------- ----------- ------------ ------------
	epsilon_sf = num[10]; // [K]
	// ---------- ----------- ------------ ------------
	sigma_sf = num[11]; // [nm]
	if ( rc == 0.0 ) { 
		rc = 5.0*sigma_ff;
		std::cout << "cut off, rc = " << rc << " [nm] (for fluid)" << std::endl;
		rcsf = 5.0*sigma_sf;
		std::cout << "cut off, rcsf = " << rcsf << " [nm] (for solid-fluid)" << std::endl;
		if ( rcsf > rc ) { rc = rcsf; }
		std::cout << "autoset (cut off) = " << rc << " [nm]" << std::endl;
	}
	std::cout << "--------------------------------------------------" << std::endl;
	// ---------- ----------- ------------ ------------
	delta = num[12]; // nm, delta < 0.3*sigma_ff for QSDFT
	// ---------- ----------- ------------ ------------
	rho_ss = num[13]; // [nm^-3], [mulecules/nm3]
	// ---------- ----------- ------------ ------------
	m = num[14]; // [kg], fluid
	// ---------- ----------- ------------ ------------
	T = num[15]; // [K]
	if ( d_hs == 0.0 ) { 
		d_hs = d_bh_calc(epsilon_ff, sigma_ff);
	} else {
		std::cout << "The hard-sphere radius of fluid =" << (d_hs/2.0) << " [nm]" << std::endl;
	}
	// ---------- ----------- ------------ ------------
	rho_b0 = num[16];
	// ---------- ----------- ------------ ------------
	h0 = num[17]; // [nm]
	// ---------- ----------- ------------ ------------
	Ris = num[18]; // [nm], the hard-sphere radius of solid
	std::cout << "The hard-sphere radius of solid =" << Ris << " [nm]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	// ---------- ----------- ------------ ------------
	ms = num[19]; // [kg], solid
	// ---------- ----------- ------------ ------------
	
	ze = calc_ze(2000);
	std::cout << "The edge position of the solid wall (one side), ze = " << ze << " [nm]" << std::endl;
	double ze_ssf;
	ze_ssf = (ze+sigma_sf);
	std::cout << "ze+sigma_sf = " << ze_ssf << " [nm]" << std::endl;
	std::cout << "(2.0*ze+sigma_sf) = " << (2.0*ze+sigma_sf) << " [nm]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	
	if ( nstep == 0 ) {
		//nstep = int((H-2.0*h0)/0.02 + 0.5);
		nstep = int((H-(2.0*ze+sigma_sf))/0.005 + 0.5);
		if ( nstep%2 == 1 ){
			nstep = nstep + 1;
		}
		std::cout << "autoset nstep = " << nstep << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
	}
	
	// ---------- ----------- ------------ ------------
	
	w_pw = (H-(2.0*ze+sigma_sf)); // pore width [nm]
	dr = (H-(2.0*ze+sigma_sf))/double(nstep-1);
	rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
	rmsf = 1.12246205*sigma_sf; // 2^(1/6)=1.12246205
	
	// ---------- ----------- ------------ ------------
	
	//ndmesh = int(2*d_hs*nrmesh/rc); // why ? this setting occures nan.
	//if ( ndmesh < 9 ) { 
	//	ndmesh = 9;
	//	std::cout << "autoset ndmesh = " << ndmesh << std::endl;
	//}
	//if ( ndmesh%2 == 0 ) { ndmesh = ndmesh + 1; }
	//dd = 2.0*d_hs/double(ndmesh-1); // rho_si(), xi()
	ndmesh = nrmesh;
	drc = rc/double(nrmesh-1); // xi()
	dd = drc;
	
	// ---------- ----------- ------------ ------------
	
	// thermal de Broglie wavelength of fluid
	//lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm], Maxwell_construction()
	lam = h/std::sqrt(2.0*M_PI*m*kb*T)*1e9; //[nm], Maxwell_construction()
	
	// thermal de Broglie wavelength of solid
	//double ms = 12.0107/(6.02214076e23)/1000; // 1.99442366e-26 [kg]
	//double lams = h/std::pow((2.0*M_PI*ms*kb*T),0.5)*1e9;
	lams = h/std::pow((2.0*M_PI*ms*kb*T),0.5)*1e9;
	
	// alpha = integal phi_att_ff * -1.0
	alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
		( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
	// rm = rm when the potential is split according to the WCA schem and rm = simga_ff when the LJ potential is split according to the BH decomposition.
	
	std::cout << "thermal de Broglie wavelength of fluid = " << lam << " [nm]" << std::endl;
	std::cout << "thermal de Broglie wavelength of solid = " << lams << " [nm]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "integal phi_att_ff * -1.0 = alpha = " << alpha << std::endl;
}

// The attractive potentials of fluid-fluid interactions.
double phi_att_ff(double r){
	double e;
	// WCA (Weeks-Chandler-Anderson) type
	if (r < rm){
		e = - epsilon_ff;
	}else if (rm <= r && r <= rc){
		// Lennard-Jones（LJ) potential
		//e = 4.0*epsilon_ff*( std::pow((sigma_ff/r),12.0) - std::pow((sigma_ff/r),6.0) );
		e = std::pow((sigma_ff/r),6.0);
		e = 4.0*epsilon_ff*( e*e - e );
	//}else {
	//	e = 0.0;
	}
	//std::cout << e << std::endl;
	return e;
}

// The attractive potentials of solid-fluid interactions.
double phi_att_sf(double r){
	double e;
	// WCA (Weeks-Chandler-Anderson) type
	if (r < rmsf){
		e = - epsilon_sf;
	}else if (rmsf <= r && r <= rcsf){
		// Lennard-Jones（LJ) potential
		//e = 4.0*epsilon_sf*( std::pow((sigma_sf/r),12.0) - std::pow((sigma_sf/r),6.0) );
		e = std::pow((sigma_sf/r),6.0);
		e = 4.0*epsilon_sf*( e*e - e );
	//}else {
	//	e = 0.0;
	}
	//std::cout << e << std::endl;
	return e;
}

// For NLDFT
// Percus-Yevick (PY) two-particle direct correlation function of the homogeneous hard-sphere fluid
//double wi(double r, int i){
//	double wi_out;
//	double rpdhs;
//	switch(i){
//		case 0:
//			if (r <= d_hs){
//				//wi_out = 3.0/(4.0*M_PI*std::pow(d_hs,3.0));
//				wi_out = 3.0/(4.0*M_PI*(d_hs*d_hs*d_hs));
//			} else {
//				wi_out = 0.0;
//			}
//			break;
//		case 1:
//			rpdhs = r/d_hs;
//			if (r <= d_hs) {
//				//wi_out = 0.475-0.648*(r/d_hs)+0.113*std::pow((r/d_hs),2.0);
//				//wi_out = 0.475-0.648*(r/d_hs)+0.113*((r/d_hs)*(r/d_hs));
//				wi_out = 0.475-0.648*(rpdhs)+0.113*(rpdhs*rpdhs);
//			} else if (d_hs < r && r <= 2.0*d_hs) {
//				//wi_out = 0.288*(d_hs/r)-0.924+0.764*(r/d_hs)-0.187*std::pow((r/d_hs),2.0);
//				//wi_out = 0.288*(d_hs/r)-0.924+0.764*(r/d_hs)-0.187*((r/d_hs)*(r/d_hs));
//				wi_out = 0.288/rpdhs-0.924+0.764*(rpdhs)-0.187*(rpdhs*rpdhs);
//			} else {
//				wi_out = 0.0;
//			}
//			break;
//		case 2:
//			rpdhs = r/d_hs;
//			if (r <= d_hs) {
				//wi_out = 5.0*M_PI*std::pow(d_hs,3.0)/144.0 * (6.0-12.0*(r/d_hs)+5.0*std::pow((r/d_hs),2.0));
				//wi_out = 5.0*M_PI*(d_hs*d_hs*d_hs)/144.0 * (6.0-12.0*(r/d_hs)+5.0*((r/d_hs)*(r/d_hs)));
//				wi_out = 5.0*M_PI*(d_hs*d_hs*d_hs)/144.0 * (6.0-12.0*(rpdhs)+5.0*(rpdhs*rpdhs));
//			} else {
//				wi_out = 0.0;
//			}
//			break;
//		default:
//			std::cout << "Error: " << i << std::endl;
//			break;
//	}
//	return wi_out;
//}

// Tarazona theory for NLDFT
//double rho_si(double *rho, double r1, double *r, int i){
//	int j,k;
//	double ra;
//	double raj;
//	double rak;
	//double ndmesh = 2*d_hs*nrmesh/rc;
	//double dd = 2.0*d_hs/double(ndmesh-1);
	//dd = drc;
//	double tpidd = 2.0*M_PI*dd;
//	double rho_si_out;
//	double rho_si_int_j[nstep];
//	double rho_si_int_k[nrmesh];
//	rho_si_int_k[0] = 0.0;
//	for (j=0; j<nstep; j++) {
//		raj = (r1-r[j]);
//		for (k=1; k<ndmesh; k++) {
//			rak = dd*double(k);
			//ra = std::pow((r1-r[j]),2.0) + std::pow((double(k)*dd),2.0);
			//ra = (r1-r[j])*(r1-r[j]) + (double(k)*dd)*(double(k)*dd);
//			ra = raj*raj + rak*rak;
			//ra = std::pow(ra,0.5);
//			ra = std::sqrt(ra);
			//std::cout << ra << std::endl;
			//
			//rho_si_int_k[k] = rho[j]*wi(ra,i)*(2.0*M_PI*(double(k)*dd)); // old ver.1.1.0
			//rho_si_int_k[k] = wi(ra,i)*(2.0*M_PI*(double(k)*dd));
//			rho_si_int_k[k] = wi(ra,i)*(tpidd*double(k));
//		}
		//integral_simpson(double *f, int n, double dx)
		//rho_si_int_j[j] = integral_simpson(rho_si_int_k, ndmesh, dd); // old ver.1.1.0
//		rho_si_int_j[j] = rho[j]*integral_simpson(rho_si_int_k, ndmesh, dd);
//	}
	//integral_simpson(double *f, int n, double dx)
//	rho_si_out = integral_simpson(rho_si_int_j, nstep-1, dr);
	//
//	return rho_si_out;
//}

// smoothed density approximation (SDA) for NLDFT
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

// smoothed density approximation (SDA), modified version for NLDFT v.0.9.0
//double rho_s(double *rho, double *r, double *rho_sj, double *rho_s0j, double *rho_s1j, double *rho_s2j){
//	int j;
//	double rho_den1j, rho_den2j;
//	for (j=0; j<nstep; j++) {
//		rho_s0j[j] = rho_si(rho, r[j], r, 0);
//		rho_s1j[j] = rho_si(rho, r[j], r, 1);
//		rho_s2j[j] = rho_si(rho, r[j], r, 2);
//		//rho_den1j = std::pow((1.0 - rho_s1j[j]),2.0);
//		rho_den1j = (1.0 - rho_s1j[j]);
//		rho_den1j = rho_den1j * rho_den1j;
//		//rho_den2j = std::pow((rho_den1j - 4.0*rho_s0j[j]*rho_s2j[j]),0.5);
//		//rho_den2j = std::sqrt(rho_den1j - 4.0*rho_s0j[j]*rho_s2j[j]);
//		rho_den2j = rho_den1j - 4.0*rho_s0j[j]*rho_s2j[j];
//		// to avoide nan
//		if ( rho_den2j > 0 ) {
//			rho_den2j = std::sqrt(rho_den2j);
//		} else {
//			rho_den2j = 0.0;
//		}
//		rho_sj[j] = 2.0*rho_s0j[j]/(1.0 - rho_s1j[j]+rho_den2j);
//		//std::cout << j << ", " << rho[j] << ", " << rho_sj[j] << ", " << rho_s0j[j] << ", " << rho_s1j[j] << ", " << rho_s2j[j] << std::endl;
//		//std::cout << rho_den1j << ", " << rho_den2j << std::endl;
//	}
//	return 0;
//}

// Steele 10-4-3 potential for NLDFT
//double phi_sf(double z){
//	double phi_sf_out;
//	double sigma_sf2 = sigma_sf*sigma_sf;
//	double sfpz = (sigma_sf/z);
//	double sfpz2 = sfpz*sfpz;
//	double dez = (0.61*delta+z);
	//phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*std::pow(sigma_sf,2.0)*delta*
	//			( (2.0/5.0)*std::pow((sigma_sf/z),10.0)-std::pow((sigma_sf/z),4.0)-std::pow(sigma_sf,4.0)/
	//			(3.0*delta*std::pow((0.61*delta+z),3.0)) );
//	phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*(sigma_sf2)*delta*
//				( (2.0/5.0)*std::pow(sfpz2,5.0)-(sfpz2*sfpz2)-(sigma_sf2*sigma_sf2)/
//				(3.0*delta*(dez*dez*dez)) );
//	return phi_sf_out;
//}

// e.g., wall potential (Carbon slit) for NLDFT
//double phi_ext(double z){
//	double phi_ext_out;
//	phi_ext_out = phi_sf(z) + phi_sf(H-z);
//	//std::cout << phi_ext_out << std::endl;
//	return phi_ext_out;
//}

// from Carnahan-Starling (CS) equation of state for NLDFT
//double mu_ex(double rho_b){
//	double y, mu_ex_out;
//	//y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
//	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
//	double den1y = (1.0-y);
//	//mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/std::pow((1.0-y),3.0);
//	//mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/((1.0-y)*(1.0-y)*(1.0-y));
//	mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/(den1y*den1y*den1y);
//	return mu_ex_out;
//}

// The excess hard sphere chemical potential (mu_ex) in the bulk fulid.
// mu_ex is calculated by the PY equation.
double mu_ex(double rho_b){
	double y, mu_ex_out;
	//eta = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	//eta = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
	//double den1e = (1.0-eta);
	double den1y = (1.0-y);
	//mu_ex_out = kb1*T*(-std::log(1-eta) + eta*(14.0 - 13.0*eta + 5.0*eta*eta)/(2.0*std::pow((1.0-eta),3.0)));
	//mu_ex_out = kb1*T*(-std::log(den1e) + eta*(14.0 - 13.0*eta + 5.0*eta*eta)/(2.0*(den1e*den1e*den1e)));
	mu_ex_out = kb1*T*(-std::log(den1y) + y*(14.0 - 13.0*y + 5.0*y*y)/(2.0*den1y*den1y*den1y));
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

// For NLDFT
//double f_ex(double rho_s){
//	double eta, f_ex_out;
	//eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
//	eta = M_PI*rho_s*(d_hs*d_hs*d_hs)/6.0;
//	double den1e = (1.0-eta);
	//f_ex_out = kb1*T*eta*(4.0-3.0*eta)/std::pow((1.0-eta),2.0);
	//f_ex_out = kb1*T*eta*(4.0-3.0*eta)/((1.0-eta)*(1.0-eta));
//	f_ex_out = kb1*T*eta*(4.0-3.0*eta)/(den1e*den1e);
//	return f_ex_out;
//}

// d(f_ex)/d(rho_s) for NLDFT
//double dfex_per_drhos(double rho_s){
//	double dfex_per_drhos_out;
//	double eta;
	//eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
//	eta = M_PI*rho_s*(d_hs*d_hs*d_hs)/6.0;
//	double den1e = (1.0-eta);
	//dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/std::pow((1.0-eta),3.0)*M_PI*std::pow(d_hs,3.0)/6.0;
	//dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/((1.0-eta)*(1.0-eta)*(1.0-eta))*M_PI*(d_hs*d_hs*d_hs)/6.0;
//	dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/(den1e*den1e*den1e)*(M_PI*(d_hs*d_hs*d_hs)/6.0);
//	return dfex_per_drhos_out;
//}

// d(rho_s)/d(rho) for NLDFT
//double drhos_per_drho(double *rho, double r1, double r2, double *r, double ra){
//	double w, drhos_per_drho_out;
//	// Percus-Yevick approximation, Tarazona theory
//	w = wi(ra,0) + wi(ra,1)*rho_s(rho,r1,r) + wi(ra,2)*std::pow(rho_s(rho,r1,r),2.0);
//	drhos_per_drho_out = w/(1.0-rho_si(rho,r2,r,1)-2.0*rho_si(rho,r2,r,2)*rho_s(rho,r2,r));
//	return drhos_per_drho_out;
//}

// d(rho_s)/d(rho), modified version for NLDFT v.0.9.0
//double drhos_per_drho_j(double ra, double rho_sj, double rho_s1j, double rho_s2j){
//	double w, drhos_per_drho_out;
//	// Percus-Yevick approximation, Tarazona theory
//	//w = wi(ra,0) + wi(ra,1)*rho_sj + wi(ra,2)*std::pow(rho_sj,2.0);
//	w = wi(ra,0) + wi(ra,1)*rho_sj + wi(ra,2)*(rho_sj*rho_sj);
//	drhos_per_drho_out = w/(1.0-rho_s1j-2.0*rho_s2j*rho_sj);
//	return drhos_per_drho_out;
//}

double ni_wall(double *r, double *n0_wall_i, double *n1_wall_i, double *n2_wall_i, double *n3_wall_i, double *nv1_wall_i, double *nv2_wall_i) {
	//
	int i;
	int w;
	double rai;
	double xs, xs2;
	double n0, n1, n2, n3, nv1, nv2;
	//
	int nwstep = 500;
	double dw = (ze+sigma_sf/2.0)/nwstep;
	//
	double n0_wall_w[nwstep];
	double n1_wall_w[nwstep];
	double n2_wall_w[nwstep];
	double n3_wall_w[nwstep];
	double nv1_wall_w[nwstep];
	double nv2_wall_w[nwstep];
	//
	// left
	for (i=0; i<nstep; i++) {
		for (w=0; w<nwstep; w++) {
			rai = (dw*double(w)-r[i]);
			//
			xs2 = (Ris*Ris-rai*rai);
			if ( xs2 >= 0.0 ){
				xs = std::sqrt(xs2);
			} else{
				xs = 0.0;
			}
			//
			n0_wall_w[w] = rho_ssq(dw*double(w))/(2.0*Ris*Ris)*xs;
			n1_wall_w[w] = rho_ssq(dw*double(w))/(2.0*Ris)*xs;
			n2_wall_w[w] = rho_ssq(dw*double(w))*(2.0*M_PI*xs);
			n3_wall_w[w] = rho_ssq(dw*double(w))*(M_PI*xs*xs);
			nv1_wall_w[w] = rho_ssq(dw*double(w))/(2.0*Ris)*(rai/Ris)*xs;
			nv2_wall_w[w] = rho_ssq(dw*double(w))*(rai/Ris)*(2.0*M_PI*xs);
		}
		n0_wall_i[i] = integral_simpson(n0_wall_w, nwstep-1, dw);
		n1_wall_i[i] = integral_simpson(n1_wall_w, nwstep-1, dw);
		n2_wall_i[i] = integral_simpson(n2_wall_w, nwstep-1, dw);
		n3_wall_i[i] = integral_simpson(n3_wall_w, nwstep-1, dw);
		nv1_wall_i[i] = integral_simpson(nv1_wall_w, nwstep-1, dw);
		nv2_wall_i[i] = integral_simpson(nv2_wall_w, nwstep-1, dw);
	}
	// right
	for (i=0; i<nstep; i++) {
		for (w=0; w<nwstep; w++) {
			rai = ((H-dw*double(w))-r[i]);
			//
			xs2 = (Ris*Ris-rai*rai);
			if ( xs2 >= 0.0 ){
				xs = std::sqrt(xs2);
			} else{
				xs = 0.0;
			}
			//
			n0_wall_w[w] = rho_ssq(H-dw*double(w))/(2.0*Ris*Ris)*xs;
			n1_wall_w[w] = rho_ssq(H-dw*double(w))/(2.0*Ris)*xs;
			n2_wall_w[w] = rho_ssq(H-dw*double(w))*(2.0*M_PI*xs);
			n3_wall_w[w] = rho_ssq(H-dw*double(w))*(M_PI*xs*xs);
			nv1_wall_w[w] = rho_ssq(H-dw*double(w))/(2.0*Ris)*(rai/Ris)*xs;
			nv2_wall_w[w] = rho_ssq(H-dw*double(w))*(rai/Ris)*(2.0*M_PI*xs);
		}
		n0_wall_i[i] = n0_wall_i[i] + integral_simpson(n0_wall_w, nwstep-1, dw);
		n1_wall_i[i] = n1_wall_i[i] + integral_simpson(n1_wall_w, nwstep-1, dw);
		n2_wall_i[i] = n2_wall_i[i] + integral_simpson(n2_wall_w, nwstep-1, dw);
		n3_wall_i[i] = n3_wall_i[i] + integral_simpson(n3_wall_w, nwstep-1, dw);
		nv1_wall_i[i] = nv1_wall_i[i] + integral_simpson(nv1_wall_w, nwstep-1, dw);
		nv2_wall_i[i] = nv2_wall_i[i] + integral_simpson(nv2_wall_w, nwstep-1, dw);
	}
	return 0;
}

double ni(double *rho, double *r, int i, double *n0_j, double *n1_j, double *n2_j, double *n3_j, double *nv1_j, double *nv2_j,
		  double *n0, double *n1, double *n2, double *n3, double *nv1, double *nv2,
		  double *n0_wall_i, double *n1_wall_i, double *n2_wall_i, double *n3_wall_i, double *nv1_wall_i, double *nv2_wall_i){
	int j;
	double raj;
	double xf, xs, xf2, xs2;
	double Rif;
	Rif = d_hs/2.0; // [nm], Rif is the hard-sphere radius of fluid
	//double Ris;
	//Ris = 0.2217/2.0; // [nm], Ris is the hard-sphere radius of solid (for QSDFT)
	//
	// Memo
	// x = y = sqrt(Ri^2-z^2)
	// z = r2 - r1 on z-axis
	//n0 = integral rho(r)/(4.0*M_PI*Ri*Ri) dr, r=Ri
	//n0 = integral rho(z)/(4.0*M_PI*Ri*Ri)*(2.0*M_PI*x) dz, z<Ri, x^2+z^2=Ri^2
	//n1 = integral rho(r)/(4.0*M_PI*Ri) dr, r=Ri
	//n1 = integral rho(z)/(4.0*M_PI*Ri)*(2.0*M_PI*x) dz, z<Ri, x^2+z^2=Ri^2
	//n2 = integral rho(r) dr, r=Ri
	//n2 = integral rho(z)*(2.0*M_PI*x) dz, z<Ri, x^2+z^2=Ri^2
	//n3 = integral rho(r) dr, r < Ri
	//n3 = integral rho(z)*(2.0*M_PI*x) dxdz, z<Ri, x^2+z^2<Ri^2
	//nv1 = integral rho(r)/(4.0*M_PI*Ri)*(Rvi/Ri) dr, r=Ri, Rvi = vector Ri
	//nv1 = integral rho(z)/(4.0*M_PI*Ri)*(z/Ri)*(2.0*M_PI*x) dxdr, z<Ri, x^2+z^2=Ri^2, Rvi = vector Ri
	//nv2 = integral rho(r)*(Rvi/Ri) dr, r=Ri, Rvi = vector Ri
	//nv2 = integral rho(z)*(z/Ri)*(2.0*M_PI*x) dxdr, z<Ri, x^2+z^2=Ri^2, Rvi = vector Ri
	//
	for (j=0; j<nstep; j++) {
		raj = (r[j]-r[i]);
		xf2 = (Rif*Rif-raj*raj);
		if ( xf2 >= 0.0 ){
			xf = std::sqrt(xf2);
		} else{
			xf = 0.0;
		}
		xs2 = (Ris*Ris-raj*raj);
		if ( xs2 >= 0.0 ){
			xs = std::sqrt(xs2);
		} else{
			xs = 0.0;
		}
		//
		//n0_j[j] = (rho[j])/(4.0*M_PI*Rif*Rif)*(2.0*M_PI*x);
		//n0_j[j] = (rho[j])/(2.0*Rif*Rif)*x;
		n0_j[j] = (rho[j])/(2.0*Rif*Rif)*xf + (rho_ssq(r[j])+rho_ssq(H-r[j]))/(2.0*Ris*Ris)*xs;
		//n0_j[j] = (rho[j])/(2.0*Rif*Rif)*xf;
		//
		//n1_j[j] = (rho[j])/(4.0*M_PI*Rif)*(2.0*M_PI*x);
		//n1_j[j] = (rho[j])/(2.0*Rif)*x;
		n1_j[j] = (rho[j])/(2.0*Rif)*xf + (rho_ssq(r[j])+rho_ssq(H-r[j]))/(2.0*Ris)*xs;
		//n1_j[j] = (rho[j])/(2.0*Rif)*xf;
		//
		n2_j[j] = (rho[j])*(2.0*M_PI*xf) + (rho_ssq(r[j])+rho_ssq(H-r[j]))*(2.0*M_PI*xs);
		//n2_j[j] = (rho[j])*(2.0*M_PI*xf);
		//
		n3_j[j] = (rho[j])*(M_PI*xf*xf) + (rho_ssq(r[j])+rho_ssq(H-r[j]))*(M_PI*xs*xs);
		//n3_j[j] = (rho[j])*(M_PI*xf*xf);
		//
		//nv1_j[j] = (rho[j])/(4.0*M_PI*Rif)*(raj/Rif)*(2.0*M_PI*x);
		nv1_j[j] = (rho[j])/(2.0*Rif)*(raj/Rif)*xf + (rho_ssq(r[j])+rho_ssq(H-r[j]))/(2.0*Ris)*(raj/Ris)*xs;
		//nv1_j[j] = (rho[j])/(2.0*Rif)*(raj/Rif)*xf;
		//
		nv2_j[j] = (rho[j])*(raj/Rif)*(2.0*M_PI*xf) + (rho_ssq(r[j])+rho_ssq(H-r[j]))*(raj/Ris)*(2.0*M_PI*xs);
		//nv2_j[j] = (rho[j])*(raj/Rif)*(2.0*M_PI*xf);
		
		//
		//std::cout << i << ", " << j << ", " << r[i] << ", " << r[j] << ", " << raj << ", " << x << std::endl;
		//std::cout << "i, j, rho[j], n0_j[j], n1_j[j], n2_j[j], n3_j[j], nv1_j[j], nv2_j[j]" << std::endl;
		//std::cout << i << ", " << j << ", " << rho[j] << ", " << n0_j[j] << ", " << n1_j[j] << ", " << n2_j[j] << ", " << n3_j[j] << ", " << nv1_j[j] << ", " << nv2_j[j] << std::endl;
	}
    //integral_trapezoidal(double *f, int n, double dx)
	//n0[i] = integral_trapezoidal(n0_j, nstep-1, dr) + n0_wall_i[i];
	//n1[i] = integral_trapezoidal(n1_j, nstep-1, dr) + n1_wall_i[i];
	//n2[i] = integral_trapezoidal(n2_j, nstep-1, dr) + n2_wall_i[i];
	//n3[i] = integral_trapezoidal(n3_j, nstep-1, dr) + n3_wall_i[i];
	//nv1[i] = integral_trapezoidal(nv1_j, nstep-1, dr) + nv1_wall_i[i];
	//nv2[i] = integral_trapezoidal(nv2_j, nstep-1, dr) + nv2_wall_i[i];
	//
	//integral_simpson(double *f, int n, double dx)
	n0[i] = integral_simpson(n0_j, nstep-1, dr) + n0_wall_i[i];
	n1[i] = integral_simpson(n1_j, nstep-1, dr) + n1_wall_i[i];
	n2[i] = integral_simpson(n2_j, nstep-1, dr) + n2_wall_i[i];
	n3[i] = integral_simpson(n3_j, nstep-1, dr) + n3_wall_i[i];
	nv1[i] = integral_simpson(nv1_j, nstep-1, dr) + nv1_wall_i[i];
	nv2[i] = integral_simpson(nv2_j, nstep-1, dr) + nv2_wall_i[i];
	//
	//double in2, inv2;
	//double rai;
	//rai = (r[i]-r[0]);
	//if ( rai <= Ris ) {
	//	// integral sin(x) dx = -cos(x) + C
	//	in2 = rho_ss*(2.0*M_PI*Ris)*-(rai/Ris-1.0)*Ris;
	//	xs = std::sqrt(Ris*Ris-rai*rai);
	//	n0[i] = n0[i] + in2/(4.0*M_PI*Ris*Ris);
	//	n1[i] = n1[i] + in2/(4.0*M_PI*Ris);
	//	n2[i] = n2[i] + in2;
	//	// integral sin(x)*sin(x) dx = (1/2)*x - (1/4)*sin(2x) + C
	//	// sin(2x) = 2*sin(x)*cos(x)
	//	n3[i] = n3[i] + rho_ss*(M_PI*Ris*Ris)*(0.5*std::asin(xs/Ris) - 0.25*2.0*(xs/Ris)*(rai/Ris))*Ris;
	//	// integral sin(x)*cos(x) dx = integral sin(2x)/2 dx = -(1/4)*cos(2x) + C
	//	// cos(2x) = 1 - 2*sin(x)*sin(x)
	//	inv2 = rho_ss*(2.0*M_PI*Ris)*(-0.25*(1.0-2.0*(xs/Ris)*(xs/Ris))+0.25)*Ris;
	//	nv1[i] = nv1[i] - inv2/(4.0*M_PI*Ris);
	//	nv2[i] = nv2[i] - inv2;
	//}
	//rai = (r[nstep-1]-r[i]);
	//if ( rai <= Ris ) {
	//	// integral sin(x) dx = -cos(x) + C
	//	in2 = rho_ss*(2.0*M_PI*Ris)*-(rai/Ris-1.0)*Ris;
	//	xs = std::sqrt(Ris*Ris-rai*rai);
	//	n0[i] = n0[i] + in2/(4.0*M_PI*Ris*Ris);
	//	n1[i] = n1[i] + in2/(4.0*M_PI*Ris);
	//	n2[i] = n2[i] + in2;
	//	// integral sin(x)*sin(x) dx = (1/2)*x - (1/4)*sin(2x) + C = (1/2)*x - (1/2)*sin(x)*cos(x) + C
	//	// sin(2x) = 2*sin(x)*cos(x)
	//	n3[i] = n3[i] + rho_ss*(M_PI*Ris*Ris)*(0.5*std::asin(xs/Ris) - 0.5*(xs/Ris)*(rai/Ris))*Ris;
	//	// integral sin(x)*cos(x) dx = integral sin(2x)/2 dx = -(1/4)*cos(2x) + C
	//	// cos(2x) = 1 - 2*sin(x)*sin(x)
	//	inv2 = rho_ss*(2.0*M_PI*Ris)*(-0.25*(1.0-2.0*(xs/Ris)*(xs/Ris))+0.25)*Ris;
	//	nv1[i] = nv1[i] + inv2/(4.0*M_PI*Ris);
	//	nv2[i] = nv2[i] + inv2;
	//}
	//
	// debug
	//std::cout << "-------------------------------------------------------------------------------------------------" << std::endl;
	//std::cout << "i=" << i << ", r[i]=" << r[i] << ", raj=" << raj << ", xs=" << xs << std::endl;
	//std::cout << "i=" << i << ", n0[i]=" << n0[i] << ", n1[i]=" << n1[i] << ", n2[i]=" <<  n2[i] << ", n3[i]=" << n3[i] << std::endl;
	//std::cout << "i=" << i << ", nv1[i]=" << nv1[i] << ", nv2[i]=" <<  nv2[i] << std::endl;
	return 0;
}

// c(1)(r) = (-k*T)*dF_HS(rho)/drho(r) = - integral f_ex dr' = dfex
// f_ex = sigma dphi/dn * w(r-r')
// phi = phi1 + phi2 + phi3 (RSLT2, the second modification)
// phi1 = (-n0*ln(1-n3)
// phi2 = (n1*n2-nv1*nv2)/(1-n3)
// phi3 = n2*n2*n2/(24*pi*(1-n3)*(1-n3)) * (1-3*|nv2/n2|*|nv2/n2|+2*|nv2/n2|*|nv2/n2|*|nv2/n2|)
// Only the z-axis components of nv1 and nv2 (with positive and negative) are remained due to symmetry.
double dfex(double *r, int i, double *n0, double *n1, double *n2, double *n3, double *nv1, double *nv2){
	int j;
	double raj;
	double x, x2;
	double dfex_out;
	double dphi_per_n0, dphi_per_n0_j[nstep];
	double dphi_per_n1, dphi_per_n1_j[nstep];
	double dphi_per_n2, dphi_per_n2_j[nstep];
	double dphi_per_n3, dphi_per_n3_j[nstep];
	double dphi_per_nv1, dphi_per_nv1_j[nstep];
	double dphi_per_nv2, dphi_per_nv2_j[nstep];
	double sxi;
	double sign;
	double Rif;
	Rif = d_hs/2.0; // [nm] Rif is the hard-sphere radius of fluid
	//double Ris;
	//Ris = 0.2217/2.0; // [nm] Ris is the hard-sphere radius of solid (for QSDFT)
	// 2.217e-10 [m], The hard sphere diameter of carbon atoms
	//
	// Memo
	// df(x)/dx = [d/dx1,...,d/dxn]t * [f1(x),...,fn(x)]
	// df(r)/dr = [d/dx,d/dy,d/dz]t * [fx(r),fy(r),fz(r)]
	// dr/dr = [d/dx,d/dy,d/dz]t * [x,y,z] = [dx/dx,dy/dy,dz/dz] = [1,1,1] = [I] = 1
	// r > 0: d|r|/dr = dr/dr = [I] = 1
	// r < 0: d|r|/dr = d(-r)/dr = -[I] = -1
	//
	for (j=0; j<nstep; j++) {
		raj = (r[j]-r[i]);
		x2 = (Rif*Rif-raj*raj);
		//
		if ( x2 >= 0.0 ){
			x = std::sqrt(x2);
		} else {
			x = 0.0;
		}
		//
		if ( n2[j] > 0.0 ) {
		sxi = std::abs(nv2[j]/n2[j]);
		} else {
			sxi = 0.0;
		}
		//std::cout << j << ", sxi = " << sxi << std::endl;
		//
		// dphi/dn0
		//dphi_per_n0[j] = -std::log(1.0-n3[j])/(4.0*M_PI*Rif*Rif)*(2.0*M_PI*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
		dphi_per_n0_j[j] = -std::log(1.0-n3[j])/(2.0*Rif*Rif)*x; // PHYSICAL REVIEW E, VOLUME 64, 011602
		//
		// dphi/dn1
		//dphi_per_n1_j[j] = ( n2[j]/(1.0-n3[j]) )/(4.0*M_PI*Rif)*(2.0*M_PI*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
		dphi_per_n1_j[j] = ( n2[j]/(1.0-n3[j]) )/(2.0*Rif)*x; // PHYSICAL REVIEW E, VOLUME 64, 011602
		//
		// dphi/dn2
		if ( nv2[j]/n2[j] >= 0.0 ){
			// nv2/n2 > 0  ->  nv2/n2 = sxi
			sign = 1.0;
		}else if ( nv2[j]/n2[j] < 0.0 ) {
			// nv2/n2 < 0  ->  -nv2/n2 = -sxi
			sign = -1.0;
		}
		// dphi/dn2 // Cite as: J. Chem. Phys. 98, 8126 (1993); https://doi.org/10.1063/1.464569
		//dphi_per_n2_j[j] = ( n1[j]/(1.0-n3[j])
		//	+ (1.0/(8.0*M_PI))*n2[j]*n2[j]/((1.0-n3[j])*(1.0-n3[j])) 
		//	- (1.0/(8.0*M_PI))*nv2[j]*nv2[j]/((1.0-n3[j])*(1.0-n3[j])) 
		//)*(2.0*M_PI*x);
		//
		// dphi/dn2, q=2 case, RSLT version // PHYSICAL REVIEW E 64 011602
		//dphi_per_n2_j[j] = ( n1[j]/(1.0-n3[j])
		//	+ 3.0*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
		//	+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//	* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(-nv2[j]/(n2[j]*n2[j])*sign)
		//)*(2.0*M_PI*x);
		//
		// dphi/dn2, q=3 case, RSLT version // PHYSICAL REVIEW E 64 011602
		//dphi_per_n2_j[j] = ( n1[j]/(1.0-n3[j])
		//	+ 3.0*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)*(1.0-sxi*sxi)
		//	+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//	* 3.0*(1.0-sxi*sxi)*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(-nv2[j]/(n2[j]*n2[j])*sign)
		//)*(2.0*M_PI*x);
		//
		// dphi/dn2, RSLT2 version // PHYSICAL REVIEW E 64 011602
		dphi_per_n2_j[j] = ( n1[j]/(1.0-n3[j])
			+ 3.0*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j]))*(1.0-3.0*sxi*sxi+2.0*sxi*sxi*sxi*sign)
			+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
				* (1.0-6.0*sxi*sign+6.0*sxi*sxi)*(-nv2[j]/(n2[j]*n2[j])*sign)
		)*(2.0*M_PI*x);
		//
		// dphi/dn2, The modified fundamental-measure theory (MFMT) // Langmuir 2008, 24, 12431-12439
		//dphi_per_n2_j[j] = ( n1[j]/(1.0-n3[j])
		//	+ 3.0*n2[j]*n2[j]*std::log(1.0-n3[j])/(36.0*M_PI*n3[j]*n3[j])
		//	+ 3.0*n2[j]*n2[j]/(36.0*M_PI*n3[j]*(1.0-n3[j])*(1.0-n3[j]))
		//	-nv2[j]*nv2[j]*std::log(1.0-n3[j])/(12.0*M_PI*n3[j]*n3[j])
		//	-nv2[j]*nv2[j]/(12.0*M_PI*n3[j]*(1.0-n3[j])*(1.0-n3[j]))
		//)*(2.0*M_PI*x);
		//
		// dphi/dn3 // Cite as: J. Chem. Phys. 98, 8126 (1993); https://doi.org/10.1063/1.464569
		//dphi_per_n3_j[j] = ( n0[j]/(1.0-n3[j])
		//	+ n1[j]*n2[j]/((1.0-n3[j])*(1.0-n3[j]))
		//	+ (1.0/(12.0*M_PI))*n2[j]*n2[j]*n2[j]/((1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//	- nv1[j]*nv2[j]//((1.0-n3[j])*(1.0-n3[j]))
		//	- (1.0/(12.0*M_PI))*n2[j]*nv2[j]*nv2[j]/((1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//)*(M_PI*x*x);
		//
		// dphi/dn3, q=2 case, RSLT version // PHYSICAL REVIEW E, VOLUME 64, 011602
		//dphi_per_n3_j[j] = ( n0[j]/(1.0-n3[j])
		//	+ (n1[j]*n2[j] - nv1[j]*nv2[j])/((1.0-n3[j])*(1.0-n3[j])) 
		//	+ 2.0*n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
		//)*(M_PI*x*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
		//
		// dphi/dn3, q=3 case, RSLT version // PHYSICAL REVIEW E, VOLUME 64, 011602
		//dphi_per_n3_j[j] = ( n0[j]/(1.0-n3[j])
		//	+ (n1[j]*n2[j] - nv1[j]*nv2[j])/((1.0-n3[j])*(1.0-n3[j])) 
		//	+ 2.0*n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)*(1.0-sxi*sxi)
		//)*(M_PI*x*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
		//
		// dphi/dn3, RSLT2 version // PHYSICAL REVIEW E, VOLUME 64, 011602
		dphi_per_n3_j[j] = ( n0[j]/(1.0-n3[j])
			+ (n1[j]*n2[j] - nv1[j]*nv2[j])/((1.0-n3[j])*(1.0-n3[j])) 
			+ 2.0*n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))*(1.0-3.0*sxi*sxi+2.0*sxi*sxi*sxi)
		)*(M_PI*x*x);
		//
		// dphi/dn3, The modified fundamental-measure theory (MFMT) // Langmuir 2008, 24, 12431-12439
		//dphi_per_n3_j[j] = ( n0[j]/(1.0-n3[j])
		//	//
		//	+ n1[j]*n2[j]/((1.0-n3[j])*(1.0-n3[j]))
		//	//
		//	- n2[j]*n2[j]*n2[j]/(36.0*M_PI*n3[j]*n3[j]*(1.0-n3[j]))
		//	-2.0*n2[j]*n2[j]*n2[j]*std::log(1.0-n3[j])/(36.0*M_PI*n3[j]*n3[j]*n3[j])
		//	//
		//	-n2[j]*n2[j]*n2[j]/(36.0*M_PI*n3[j]*n3[j]*(1.0-n3[j])*(1.0-n3[j]))
		//	+2.0*n2[j]*n2[j]*n2[j]/(36.0*M_PI*n3[j]*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//	//
		//	-nv1[j]*nv2[j]/((1.0-n3[j])*(1.0-n3[j]))
		//	//
		//	+n2[j]*nv2[j]*nv2[j]/(12.0*M_PI*n3[j]*n3[j]*(1.0-n3[j]))
		//	+2.0*n2[j]*nv2[j]*nv2[j]*std::log(1.0-n3[j])/(12.0*M_PI*n3[j]*n3[j]*n3[j])
		//	//
		//	+n2[j]*nv2[j]*nv2[j]/(12.0*M_PI*n3[j]*n3[j]*(1.0-n3[j])*(1.0-n3[j]))
		//	-2.0*n2[j]*nv2[j]*nv2[j]/(12.0*M_PI*n3[j]*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//)*(M_PI*x*x);
		//
		// dphi/dnv1
		//dphi_per_nv1_j[j] = ( -nv2[j]/(1.0-n3[j]) )/(4.0*M_PI*Ri)*(raj/Ri)*(2.0*M_PI*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
		dphi_per_nv1_j[j] = ( -nv2[j]/(1.0-n3[j]) )/(2.0*Rif)*(raj/Rif)*x; // PHYSICAL REVIEW E, VOLUME 64, 011602
		//
		// dphi/dnv2 // Cite as: J. Chem. Phys. 98, 8126 (1993); https://doi.org/10.1063/1.464569
		//dphi_per_nv2_j[j] = ( -nv1[j]/(1.0-n3[j])
		//	- (1.0/(4.0*M_PI))*n2[j]*nv2[j]/((1.0-n3[j])*(1.0-n3[j]))
		//)*(raj/Rif)*(2.0*M_PI*x);
		//
		// dphi/dnv2, q=2 case, RSLT version // PHYSICAL REVIEW E, VOLUME 64, 011602
		//dphi_per_nv2_j[j] = ( -nv1[j]/(1.0-n3[j])
		//	+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//	* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(1.0/n2[j])
		//)*(raj/Rif)*(2.0*M_PI*x);
		//
		// dphi/dnv2, q=3 case, RSLT version // PHYSICAL REVIEW E, VOLUME 64, 011602
		//dphi_per_nv2_j[j] = ( -nv1[j]/(1.0-n3[j])
		//	+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//	* 3.0*(1.0-sxi*sxi)*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(1.0/n2[j])
		//)*(raj/Rif)*(2.0*M_PI*x);
		//
		// dphi/dnv2, RSLT2 version // PHYSICAL REVIEW E, VOLUME 64, 011602
		dphi_per_nv2_j[j] = ( -nv1[j]/(1.0-n3[j])
			+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
				* (1.0-6.0*sxi*sign+6.0*sxi*sxi)*(sign*1.0/n2[j])
		)*(raj/Rif)*(2.0*M_PI*x);
		//
		// dphi/dnv2, The modified fundamental-measure theory (MFMT) // Langmuir 2008, 24, 12431-12439
		//dphi_per_nv2_j[j] = ( -nv1[j]/(1.0-n3[j])
		//	-2.0*n2[j]*nv2[j]*std::log(1.0-n3[j])/(12.0*M_PI*n3[j]*n3[j])
		//	-2.0*n2[j]*nv2[j]/(12.0*M_PI*n3[j]*(1.0-n3[j])*(1.0-n3[j]))
		//)*(raj/Rif)*(2.0*M_PI*x);
	}
	//
    //integral_trapezoidal(double *f, int n, double dx)
	//dphi_per_n0  = integral_trapezoidal(dphi_per_n0_j, nstep-1, dr);
	//dphi_per_n1  = integral_trapezoidal(dphi_per_n1_j, nstep-1, dr);
	//dphi_per_n2  = integral_trapezoidal(dphi_per_n2_j, nstep-1, dr);
	//dphi_per_n3  = integral_trapezoidal(dphi_per_n3_j, nstep-1, dr);
	//dphi_per_nv1 = integral_trapezoidal(dphi_per_nv1_j, nstep-1, dr);
	//dphi_per_nv2 = integral_trapezoidal(dphi_per_nv2_j, nstep-1, dr);
	//
	//integral_simpson(double *f, int n, double dx)
	dphi_per_n0 = integral_simpson(dphi_per_n0_j, nstep-1, dr);
	dphi_per_n1 = integral_simpson(dphi_per_n1_j, nstep-1, dr);
	dphi_per_n2 = integral_simpson(dphi_per_n2_j, nstep-1, dr);
	dphi_per_n3 = integral_simpson(dphi_per_n3_j, nstep-1, dr);
	dphi_per_nv1 = integral_simpson(dphi_per_nv1_j, nstep-1, dr);
	dphi_per_nv2 = integral_simpson(dphi_per_nv2_j, nstep-1, dr);
	//
	//std::cout << "i, dphi_per_n0, dphi_per_n1, dphi_per_n2, dphi_per_n3, dphi_per_nv1, dphi_per_nv2" << std::endl;
	//std::cout << i << ", " << dphi_per_n0 << "," << dphi_per_n1 << "," << dphi_per_n2 << "," << dphi_per_n3 << "," << dphi_per_nv1 << "," << dphi_per_nv2 << "," << std::endl;
	//
	dfex_out = -(dphi_per_n0 + dphi_per_n1 + dphi_per_n2 + dphi_per_n3 + dphi_per_nv1 + dphi_per_nv2);
	//std::cout << dfex_out << std::endl;
	return dfex_out;
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
	//
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
				//alpha_int_k[k]  = -phi_att_ff(ra)*(2.0*M_PI*(double(k)*drc));
				alpha_int_k[k]  = -phi_att_ff(ra)*(tpidrc*double(k));
			}
			//integral_simpson(double *f, int n, double dx)
			alpha_int_j[j]  = integral_simpson(alpha_int_k, nrmesh-1, drc);
		}
		//integral_simpson(double *f, int n, double dx)
		//alpha_other_method  = alpha_other_method + integral_simpson(alpha_int_j, nstep-1, dr)*2.0*dr;
		alpha_other_method  = alpha_other_method + integral_simpson(alpha_int_j, nstep-1, dr);
	}
	alpha_other_method  = alpha_other_method * 2.0 * dr / (H-(2.0*ze+sigma_sf));
	//std::cout << "--------------------------------------------------" << std::endl;
	//std::cout << "average alpha of other method = " << alpha_other_method << " in (carbon) slit" << std::endl;
	return alpha_other_method;
}

// fluid-fluid
double phi_att_ff_int(double *r, double *phi_att_ff_int_ij){
	int i,j,k;
	double ra;
	double raj;
	double rak;
	//double drc = rc/double(nrmesh-1);
	//dd = drc;
	double phi_ff_int_k[nrmesh];
	double tpidrc = 2.0*M_PI*drc;
	phi_ff_int_k[0] = 0.0;
	//
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
				//rho_phi_ff_int_k[k]  = rho[j]*phi_att_ff(ra)*(2.0*M_PI*(double(k)*drc)); // old ver.1.1.0
				//phi_ff_int_k[k]  = phi_att_ff(ra)*(2.0*M_PI*(double(k)*drc));
				phi_ff_int_k[k]  = phi_att_ff(ra)*(tpidrc*double(k));
			}
		phi_att_ff_int_ij[i*nstep+j] = integral_simpson(phi_ff_int_k, nrmesh-1, drc);
		}
	}
	return 0;
}

// solid-fluid
double phi_att_sf_int(double *r, double *rhos_phi_sf_int_i){
	int i,j,k;
	double ra_left;
	double ra_right;
	double raj_left;
	double raj_right;
	double rak;
	//dd = drc = rc/double(nrmesh-1);
	//
	int sfmesh = 500;
	double dsf = (h0+2.0*delta)/(sfmesh-1);
	double rhos_phi_sf_int_j[sfmesh];
	//
	int sfnrmesh = 500;
	double drcsf = rcsf/(sfnrmesh-1);
	double phi_sf_int_k[sfnrmesh];
	//
	double tpi = 2.0*M_PI;
	phi_sf_int_k[0] = 0.0;
	//
	//for (i=0; i<nstep; i++) {
	for (i=0; i<=(nstep-2)/2; i++){
		//
		for (j=0; j<sfmesh; j++) {
			raj_left  = (double(j)*dsf - r[i]);
			raj_right = ((H-double(j)*dsf) - r[i]);
			for (k=1; k<sfnrmesh; k++) {
				rak = drcsf*double(k);
				//
				ra_left  = rak*rak + raj_left*raj_left;
				ra_left  = std::sqrt(ra_left);
				//
				ra_right = rak*rak + raj_right*raj_right;
				ra_right = std::sqrt(ra_right);
				//
				phi_sf_int_k[k]  = ( phi_att_sf(ra_left) + phi_att_sf(ra_right) ) * (tpi*rak);
			}
			rhos_phi_sf_int_j[j] = rho_ssq(double(j)*dsf)*integral_simpson(phi_sf_int_k, sfnrmesh-1, drcsf);
			//rhos_phi_sf_int_j[j] = rho_ssq(double(j)*dsf)*integral_trapezoidal(phi_sf_int_k, sfnrmesh-1, drcsf);
		}
		//
		rhos_phi_sf_int_i[i] = integral_simpson(rhos_phi_sf_int_j, sfmesh-1, dsf);
		//rhos_phi_sf_int_i[i] = integral_trapezoidal(rhos_phi_sf_int_j, sfmesh-1, dsf);
		rhos_phi_sf_int_i[(nstep-1)-i] = rhos_phi_sf_int_i[i];
	}
	//for (i=0; i<nstep; i++) {
	//	std::cout << i << ", " << r[i] << ", " << rhos_phi_sf_int_i[i] << std::endl;
	//}
	return 0;
}

// xi include kb1*T*(std::log(rho_b)) type.
// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
//double xi(double *rho, double *r, int i, double rho_b, double *rho_sj, double *rho_s0j, double *rho_s1j, double *rho_s2j, double *phi_att_ff_int_ij){
double xi(double *rho, double *r, int i, double rho_b, double *phi_att_ff_int_ij, double *rho_phi_ff_int_i){
	int j;
	double rho_phi_ff_int_j[nstep];
	for (j=0; j<nstep; j++) {
		rho_phi_ff_int_j[j]  = rho[j]*phi_att_ff_int_ij[i*nstep+j];
	}
	//integral_simpson(double *f, int n, double dx)
	rho_phi_ff_int_i[i]  = integral_simpson(rho_phi_ff_int_j, nstep-1, dr);
	//
	double xi_out;
	xi_out = ( - rho_b*alpha ) + ( mu_ex(rho_b) - rho_phi_ff_int_i[i] ) + ( kb1*T*std::log(rho_b) );
	//
	return xi_out;
}

// For SDA, from Carnahan-Starling (CS) equation
//double press_hs(double rho_b){
//	double y, press_hs_out;
//	//y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
//	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
//	double den1y = (1.0-y);
//	//press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/std::pow((1.0-y),3.0);
//	//press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/((1.0-y)*(1.0-y)*(1.0-y));
//	press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/(den1y*den1y*den1y);
//	return press_hs_out;
//}

// For FMT, from Percus Yevick (PY) equation
// http://www.sklogwiki.org/SklogWiki/index.php/Exact_solution_of_the_Percus_Yevick_integral_equation_for_hard_spheres
double press_hs(double rho_b){
	double eta, press_hs_out;
	//eta = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	eta = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
	double den1e = (1.0-eta);
	//press_hs_out = rho_b*kb1*T* (1.0 + eta + eta*eta)/std::pow((1.0-eta),3.0);
	//press_hs_out = rho_b*kb1*T* (1.0 + eta + eta*eta)/((1.0-eta)*(1.0-eta)*(1.0-eta));
	press_hs_out = rho_b*kb1*T* (1.0 + eta + eta*eta)/(den1e*den1e*den1e);
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
	double rho_b0_out;
	double rho_b0_gas, rho_b0_metastable, rho_b0_liquid;
	double press_b0;
	//
	// rho_b vs. mu_b/epsilon_ff
	std::ofstream ofs("./Maxwell_construction_data.txt");
	ofs << "Chemical_potential(mu_b/epsilon_ff), Density(rho_b*d_hs^3)" << std::endl;
	for (i=0; i<iter_max_drhob0; i++){
		rho_b0_out = drhob0*double(i+1.0);
		mu_b_per_epsilon_ff[i] = mu_b(rho_b0_out)/epsilon_ff;
		//ofs << mu_b_per_epsilon_ff[i] << ", " << rho_b0_out*std::pow(d_hs,3.0) << std::endl;
		ofs << mu_b_per_epsilon_ff[i] << ", " << rho_b0_out*(d_hs*d_hs*d_hs) << std::endl;
		//std::cout << "rho_b0 = "<< rho_b0_out << ", mu_b/epsilon_ff = " << mu_b_per_epsilon_ff[i] << std::endl;
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
		rho_b0_out = drhob0*double(j+1.0);
		if (std::abs(diff) <= threshold_diff) {
			//std::cout << "mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << ", diff = " << diff << std::endl;
			break;
		}
	}
	// find rho_b0
	flag = 0;
	for (i=0; i<iter_max_drhob0; i++){
		rho_b0_out = drhob0*double(i+1.0);
		//if ( std::abs(mu_b(rho_b0_out)/epsilon_ff - mu_e_per_epsilon_ff) <= threshold_find &&
		//	 0.05 <= rho_b0_out*std::pow(d_hs,3.0) &&  rho_b0_out*std::pow(d_hs,3.0) <= 0.75) {
		if ( std::abs(mu_b(rho_b0_out)/epsilon_ff - mu_e_per_epsilon_ff) <= threshold_find ) {
			//std::cout << "rho_b0 = " << rho_b0_out << ", rho_b0*d_hs^3 = " << rho_b0_out*std::pow(d_hs,3.0) << std::endl;
			if ( flag == 0 ){
				rho_b0_gas = rho_b0_out;
				flag = 1;
			} else if ( flag == 1 && 5.0*rho_b0_gas < rho_b0_out ){
				rho_b0_metastable = rho_b0_out;
				flag = 2;
			} else if ( flag == 2 && 1.5*rho_b0_metastable < rho_b0_out ){
				rho_b0_liquid = rho_b0_out;
				break;
			}
		}
	}	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "Maxwell construction (Maxwell equal area rule)" << std::endl;
	std::cout << "chemical potential, mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << std::endl;
	rho_b0_out = rho_b0_gas;
	//std::cout << "density, rho_b0*d_hs^3 = " << rho_b0_out*std::pow(d_hs,3.0) << std::endl;
	std::cout << "density, rho_b0*d_hs^3 = " << rho_b0_out*(d_hs*d_hs*d_hs) << std::endl;
	//press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
	press_b0 = press_hs(rho_b0_out) - 0.5*(rho_b0_out*rho_b0_out)*alpha;
	std::cout << "Bulk pressure, P0 = " << press_b0 << " (rho_b0 = " << rho_b0_out << ")" <<std::endl;
	std::cout << std::endl;
	std::cout << "gas phase   : rho_b0_gas        = " << rho_b0_gas        << ", rho_b0_gas*d_hs^3        = " << rho_b0_gas*std::pow(d_hs,3.0)        << std::endl;
	std::cout << "metastable  : rho_b0_metastable = " << rho_b0_metastable << ", rho_b0_metastable*d_hs^3 = " << rho_b0_metastable*std::pow(d_hs,3.0) << std::endl;
	std::cout << "liquid phase: rho_b0_liquid     = " << rho_b0_liquid     << ", rho_b0_liquid*d_hs^3     = " << rho_b0_liquid*std::pow(d_hs,3.0)     << std::endl;
	return rho_b0_out;
}

double fex(int i, double *n0, double *n1, double *n2, double *n3, double *nv1, double *nv2){
	double fex_out;
	double sxi = std::abs(nv2[i]/n2[i]);
	double phi1, phi2, phi3;
	phi1 = -n0[i]*std::log(1.0-n3[i]);
	phi2 = (n1[i]*n2[i] - nv1[i]*nv2[i])/(1.0-n3[i]);
	//
	//phi3 = ((1.0/3.0)*n2[i]*n2[i]*n2[i] - n2[i]*(nv2[i]*nv2[i]))/(8.0*M_PI*(1.0-n3[i])*(1.0-n3[i]));
	//
	// RSLT1, q=2
	//phi3 = n2[i]*n2[i]*n2[i]/(24.0*M_PI*(1.0-n3[i])*(1.0-n3[i]))*(1.0-sxi*sxi)*(1.0-sxi*sxi);
	//
	// RSLT1, q=3
	//phi3 = n2[i]*n2[i]*n2[i]/(24.0*M_PI*(1.0-n3[i])*(1.0-n3[i]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)*(1.0-sxi*sxi);
	//
	// RSLT2 version // PHYSICAL REVIEW E 64 011602
	if ( nv2[i]/n2[i] < 0.0 ){
		sxi = sxi*-1.0;
	}
	phi3 = n2[i]*n2[i]*n2[i]/(24.0*M_PI*(1.0-n3[i])*(1.0-n3[i]))*(1.0-3.0*sxi*sxi+2.0*sxi*sxi*sxi);
	//
	fex_out = phi1 + phi2 + phi3;
	return fex_out;
}

// grand potential for FMT
// Omega(rho of fluid) + integral 0.5*rho*rhos*phi_att_sf(r-r') drdr'
double omega(double *rho, double *r, double *fex_i, double *rho_phi_ff_int_i, double *rhos_phi_sf_int_i, double rho_b){
	double omega_out;
	omega_out = 1.0;
	double omega1, omega2, omega3_ff, omega3_sf, omega4;
	int i;
	//int i,j;
	double fidf[nstep];
	double rho_x_rho_phi_ff_int[nstep];
	double rho_x_rhos_phi_sf_int[nstep];
	double rho_x_muf[nstep];
	double muf = (kb1*T)*std::log(rho_b*lam*lam*lam) + mu_ex(rho_b) - rho_b*alpha;
	//
	for (i=0; i<nstep; i++){
		fidf[i] = rho[i]*(std::log(rho[i]*lam*lam*lam)-1.0);
		rho_x_rho_phi_ff_int[i] = rho[i] * rho_phi_ff_int_i[i];
		rho_x_rhos_phi_sf_int[i] = rho[i] * rhos_phi_sf_int_i[i];
		rho_x_muf[i] = rho[i] * - muf;
	}
	omega1 = (kb1*T) * integral_simpson(fidf, nstep-1, dr);
	//
	// solid
	//int sfmesh = 2000;
	//double dsf = H/(sfmesh-1);
	//double rhos[sfmesh];
	//double fids[sfmesh];
	//double ms = 12.0107/(6.02214076e23)/1000;
	//double lams = h/std::pow((2.0*M_PI*ms*kb*T),0.5)*1e9;
	//for (j=0; j<sfmesh; j++) {
	//	rhos[i] = rho_ssq(double(j)*dsf)+rho_ssq(H-double(j)*dsf);
	//	fids[i] = rhos[i]*(std::log(rhos[i]*lam*lam*lam)-1.0);
	//}
	//
	omega2 = (kb1*T) * integral_simpson(fex_i, nstep-1, dr); // Fex (fluid + solid)
	//
	omega3_ff = 0.5 * integral_simpson(rho_x_rho_phi_ff_int, nstep-1, dr);
	omega3_sf = 0.5 * integral_simpson(rho_x_rhos_phi_sf_int, nstep-1, dr);
	//
	omega4 = integral_simpson(rho_x_muf, nstep-1, dr);
	//
	omega_out = (omega1 + omega2 + omega3_ff + omega3_sf + omega4) / epsilon_ff;
	//std::cout << omega1 << ", " << omega2 << ", " << omega3_ff << ", " << omega3_sf << ", " << omega4 << std::endl;
	return omega_out;
}

int main(){
	int i,j,k;
	double v_gamma;
	double press_b, press_b0, pp0;
	double rho_b;
	double v_mmol_per_cm3;
	double v_cm3STP_per_cm3;
	double grand_potential;
	//
	read_parameters();
	double r[nstep];
	double rho[nstep], rho_new[nstep];
	//
	for (i=0; i<nstep; i++){
		r[i] = (2.0*ze+sigma_sf)/2.0 + dr*double(i);
		//std::cout << i << ", " << r[i] << std::endl;
	}
	
	// show alpha
	//calc_alpha(r);
	// alpha = calc_alpha(r);
	
	// set rho_b0
	if ( rho_b0 != 0.0 ){
		std::cout << "rho_b0 = " << rho_b0 << std::endl;
	} else {
		rho_b0 = Maxwell_construction();
	}
	
	//std::cout << rho_b0 << std::endl;
	// initialization
	for (i=0; i<nstep; i++){
		rho[i] = rho_b0/(nstep*dr);
		rho_new[i] = 0.0;
	}
	// P/P0, V[molecules/nm^3], Omega/epsilon_ff[nm^-2]
	std::ofstream ofsppov_vs("./PP0_vs_Vgamma_data_vs.txt");
	ofsppov_vs << "# w = (H-(2.0*ze+sigma_sf)) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_vs << "# P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega(ff+sf part)/epsilon_ff[1/nm2](now dummy)" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "w = (H-(2.0*ze+sigma_sf)) = pore width = " << w_pw << " [nm]" << std::endl;
	std::cout << "P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega(ff+sf part)/epsilon_ff[1/nm2](now dummy)" << std::endl;
	//double rho_sj[nstep];
	//double rho_s0j[nstep];
	//double rho_s1j[nstep];
	//double rho_s2j[nstep];
	//double phi_att_ff_int_ij[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	double *phi_att_ff_int_ij = (double *)malloc(sizeof(double)*((nstep+1)*nstep));
	if (phi_att_ff_int_ij == NULL) {
		printf("Memory cannot be allocated.");
		std::exit(1);
	} else {
		printf("Memory has been allocated. The address is %p\n", phi_att_ff_int_ij);
	}
	phi_att_ff_int(r, phi_att_ff_int_ij); // calculate integral phi_att_ff at r[i]
	//double phi_att_sf_int_i[nstep];
	double rho_phi_ff_int_i[nstep];
	double rhos_phi_sf_int_i[nstep];
	phi_att_sf_int(r, rhos_phi_sf_int_i); // calculate integral phi_att_sf at r[i] -> rhos * phi_att_sf
	//
	double n0_wall_i[nstep];
	double n1_wall_i[nstep];
	double n2_wall_i[nstep];
	double n3_wall_i[nstep];
	double nv1_wall_i[nstep];
	double nv2_wall_i[nstep];
	ni_wall(r, n0_wall_i, n1_wall_i, n2_wall_i, n3_wall_i, nv1_wall_i, nv2_wall_i);
	//
	double n0_j[nstep], n0[nstep];
	double n1_j[nstep], n1[nstep];
	double n2_j[nstep], n2[nstep];
	double n3_j[nstep], n3[nstep];
	double nv1_j[nstep], nv1[nstep];
	double nv2_j[nstep], nv2[nstep];
	double c1;
	double fex_i[nstep]; // For grand potential, Omega
	double diff = 1.0;
	double old_diff;
	double diff0, diff1;
	double mixing;
	for (k=0; k<100; k++){
		rho_b = rho_b0 * std::exp(-(20.0-2.0*double(k+1.0)/10.0));
		//rho_b = rho_b0 * std::exp(-(20.0-2.0*double(99.0-k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		//double check_data;
		for (j=0; j<cycle_max; j++){
			for (i=0; i<nstep; i++){
				ni(rho, r, i, n0_j, n1_j, n2_j, n3_j, nv1_j, nv2_j, n0, n1, n2, n3, nv1, nv2, n0_wall_i, n1_wall_i, n2_wall_i, n3_wall_i, nv1_wall_i, nv2_wall_i);
			}
			for (i=0; i<nstep; i++){
				c1 = dfex(r, i, n0, n1, n2, n3, nv1, nv2);
				rho_new[i] = std::exp(c1+xi(rho,r,i,rho_b, phi_att_ff_int_ij,rho_phi_ff_int_i)/(kb1*T)-rhos_phi_sf_int_i[i]/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				//
				// overflow about std::exp(730)
				// to avoid overflow
				if (rho_new[i] > 1e6){
					rho_new[i] = rho[i] * 10.0;
					//std::cout << "rho[i] > 1e6" << std::endl;
					//std::exit(1);
				}
				// to avoid -inf or int
				if (rho_new[i] < 1e-6 && rho[i] < 1e-6){
					rho_new[i] = 1e-6;
					rho[i] = 1e-6;
				}
			}
			diff = 0.0;
			for (i=0; i<=(nstep-2)/2; i++){
				diff0 = std::abs((rho_new[i]-rho[i])/rho[i]);
				diff = diff + 2.0*diff0;
				mixing = wmixing + wmixing/(0.5+diff0);
				//std::cout << i << ", " << mixing << std::endl;
				rho[i] = mixing*rho_new[i] + (1.0-mixing)*rho[i];
				rho[(nstep-1)-i] = rho[i]; // The rest is filled with mirror symmetry. 
			}
			if ( diff/nstep < 0.005 && j >= 100) {
				break;
			}
		}
		//for (i=0; i<nstep; i++){
		//	std::cout << "i=" << i << ", r=" << r[i] << ", rho_new=" << rho_new[i] << ", rho=" << rho[i] << ", mixing=" << mixing << ", (diff/nstep*100.0)=" << (diff/nstep*100.0) << std::endl;
		//}
		//
		v_gamma = integral_simpson(rho, nstep-1, dr);
		//v_gamma = v_gamma/(H-sigma_ss) - rho_b; // for NLDFT
		v_gamma = v_gamma/(H-(2.0*ze+sigma_sf)) - rho_b;
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076 ) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10; // [mmol/cm3]
		v_cm3STP_per_cm3 = v_mmol_per_cm3 * 22.414;
		if (v_gamma < 0) { v_gamma = 0.0; }
		//v_gamma = v_gamma * (0.8064/28.0134/1e21*6.02214e23)/rho_b;
		// N2(77K): 0.8064 g/mL, 0.8064/28.0134 mol/mL, 0.8064/28.0134/1e21 mol/nm3, 0.8064/28.0134/1e21*6.02214e23 molecules/nm3
		//std::cout << "V= " << v_gamma << std::endl;
		//
		// press_hs(rho_b) from Carnahan-Starling (CS) equation of state
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		//std::cout << "P= " << press_b << std::endl;
		//std::cout << "P0= " << press_b0 << std::endl;
		pp0 = press_b/press_b0;
		// grand ppotential, Omega(ff+sf part)
		for (i=0; i<nstep; i++){
			fex_i[i] = fex(i, n0, n1, n2, n3, nv1, nv2);
		}
		grand_potential = omega(rho, r, fex_i, rho_phi_ff_int_i, rhos_phi_sf_int_i, rho_b);
		//grand_potential = 1.0;
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_vs << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
	}
	// reverse
	// P/P0, V[molecules/nm^3], Omega/epsilon_ff[nm^-2]
	std::ofstream ofsppov_ls("./PP0_vs_Vgamma_data_ls.txt");
	ofsppov_ls << "# w = (H-(2.0*ze+sigma_sf)) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_ls << "# P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	//std::cout << "w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	//std::cout << "P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega/epsilon_ff[1/nm2]" << std::endl;
	for (k=0; k<100; k++){
		//rho_b = rho_b0 * std::exp(-(20.0-2.0*double(k+1.0)/10.0));
		rho_b = rho_b0 * std::exp(-(20.0-2.0*double(99.0-k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		//double check_data;
		for (j=0; j<cycle_max; j++){
			for (i=0; i<nstep; i++){
				ni(rho, r, i, n0_j, n1_j, n2_j, n3_j, nv1_j, nv2_j, n0, n1, n2, n3, nv1, nv2, n0_wall_i, n1_wall_i, n2_wall_i, n3_wall_i, nv1_wall_i, nv2_wall_i);
			}
			for (i=0; i<nstep; i++){
				c1 = dfex(r, i, n0, n1, n2, n3, nv1, nv2);
				rho_new[i] = std::exp(c1+xi(rho,r,i,rho_b, phi_att_ff_int_ij,rho_phi_ff_int_i)/(kb1*T)-rhos_phi_sf_int_i[i]/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				//
				// overflow about std::exp(730)
				// to avoid overflow
				if (rho_new[i] > 1e6){
					rho_new[i] = rho[i] * 10.0;
					//std::cout << "rho[i] > 1e6" << std::endl;
					//std::exit(1);
				}
				// to avoid -inf or int
				if (rho_new[i] < 1e-6 && rho[i] < 1e-6){
					rho_new[i] = 1e-6;
					rho[i] = 1e-6;
				}
			}
			diff = 0.0;
			for (i=0; i<=(nstep-2)/2; i++){
				diff0 = std::abs((rho_new[i]-rho[i])/rho[i]);
				diff = diff + 2.0*diff0;
				mixing = wmixing + wmixing/(0.5+diff0);
				//std::cout << i << ", " << mixing << std::endl;
				rho[i] = mixing*rho_new[i] + (1.0-mixing)*rho[i];
				rho[(nstep-1)-i] = rho[i]; // The rest is filled with mirror symmetry. 
			}
			if ( diff/nstep < 0.005 && j >= 100) {
				break;
			}
		}
		//
		v_gamma = integral_simpson(rho, nstep-1, dr);
		//v_gamma = v_gamma/(H-sigma_ss) - rho_b; // for NLDFT
		v_gamma = v_gamma/(H-(2.0*ze+sigma_sf)) - rho_b;
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076 ) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10; // [mmol/cm3]
		v_cm3STP_per_cm3 = v_mmol_per_cm3 * 22.414;
		if (v_gamma < 0) { v_gamma = 0.0; }
		//v_gamma = v_gamma * (0.8064/28.0134/1e21*6.02214e23)/rho_b;
		// N2(77K): 0.8064 g/mL, 0.8064/28.0134 mol/mL, 0.8064/28.0134/1e21 mol/nm3, 0.8064/28.0134/1e21*6.02214e23 molecules/nm3
		//std::cout << "V= " << v_gamma << std::endl;
		//
		// press_hs(rho_b) from Carnahan-Starling (CS) equation of state
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		//std::cout << "P= " << press_b << std::endl;
		//std::cout << "P0= " << press_b0 << std::endl;
		pp0 = press_b/press_b0;
		// grand ppotential, Omega(ff+sf part)
		for (i=0; i<nstep; i++){
			fex_i[i] = fex(i, n0, n1, n2, n3, nv1, nv2);
		}
		grand_potential = omega(rho, r, fex_i, rho_phi_ff_int_i, rhos_phi_sf_int_i, rho_b);
		//grand_potential = 1.0;
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_ls << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
	}
	free(phi_att_ff_int_ij);
	return 0;
}
