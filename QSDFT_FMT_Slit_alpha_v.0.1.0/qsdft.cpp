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

// compiling: c++ nldft.cpp
// usage: ./a.out

// debag mode
// compiling: c++ nldft.cpp -g3 -ggdb
// run: gdb ./a.out
//      (gdb) r
//      (gdb) backtrace

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
// alpha = integal phi_att_ff * -1.0
//extern double alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
//	( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
double alpha;
// ---------- ----------- ------------ ------------
// The edge position ze of the solid wall
double ze;
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

double rho_ssq(double z){
	double rho_ssq_out;
	double h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	double rho_ss0 = 114.0; // [molecules/nm3] (the density of bulk carbon)
	double delta_ss = 0.13; // [nm] (the roughness parameter) (the half-width of the density ramp)
	if ( 0.0 <= z && z < h0 ){
		rho_ssq_out = rho_ss0;
	} else if ( h0 <= z && z < h0+2.0*delta_ss ){
		rho_ssq_out = 0.75*rho_ss0 * (1.0 - (z - h0)/(2.0*delta_ss));
	} else {
		rho_ssq_out = 0.0;
	}
	return rho_ssq_out;
}

// The edge position ze of the solid wall
double calc_ze(void){
	int i;
	double ze_out;
	double h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	double rho_ss0 = 114.0; // [molecules/nm3] (the density of bulk carbon)
	double delta_ss = 0.13; // [nm] (the roughness parameter) (the half-width of the density ramp)
	double rho_ssi[nstep];
	double dss = (2.0*delta_ss)/(nstep-1);
	for (i=0; i<nstep; i++){
		rho_ssi[i] = rho_ssq(h0+double(i)*dss);
	}
	//integral_trapezoidal(double *f, int n, double dx)
	//ze_out = integral_trapezoidal(rho_ssi, nstep, dss)/rho_ss0 + h0;
	//integral_simpson(double *f, int n, double dx)
	ze_out = integral_simpson(rho_ssi, nstep, dss)/rho_ss0 + h0;
	return ze_out;
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
		nstep = int((H-2.0*(ze+sigma_sf))/0.02 + 0.5);
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
		if ( nrmesh%2 == 0 ){
			nrmesh = nrmesh + 1;
		}
		std::cout << "autoset nrmesh = 5.0*sigma_ff " << nrmesh << std::endl;
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
	
	ze = calc_ze();
	std::cout << "ze = " << ze << " [nm]" << std::endl;
	double ze_ssf;
	ze_ssf = (ze+sigma_sf);
	std::cout << "ze+sigma_sf = " << ze_ssf << " [nm]" << std::endl;
	std::cout << "2.0*(ze+sigma_sf) = " << 2.0*ze_ssf << " [nm]" << std::endl;
	
	w_pw = (H-2.0*(ze+sigma_sf)); // pore width [nm]
	//dr = H/double(nstep-1);
	dr = (H-2.0*(ze+sigma_sf))/double(nstep-1);
	//w_pw = (H-sigma_ss); // pore width [nm]
	//dr = (H-sigma_ss)/double(nstep+1);
	rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
	
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
	
	// thermal de Broglie wavelength
	//lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm], Maxwell_construction()
	lam = h/std::sqrt(2.0*M_PI*m*kb*T)*1e9; //[nm], Maxwell_construction()
	
	// alpha = integal phi_att_ff * -1.0
	alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
		( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
	// rm = rm when the potential is split according to the WCA schem and rm = simga_ff when the LJ potential is split according to the BH decomposition.
	
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "thermal de Broglie wavelength = " << lam << " [nm]" << std::endl;
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
	}else {
		e = 0.0;
	}
	//std::cout << e << std::endl;
	return e;
}

// The attractive potentials of solid-fluid interactions.
double phi_att_sf(double r){
	double e;
	// WCA (Weeks-Chandler-Anderson) type
	if (r < rm){
		e = - epsilon_sf;
	}else if (rm <= r && r <= rc){
		// Lennard-Jones（LJ) potential
		//e = 4.0*epsilon_sf*( std::pow((sigma_sf/r),12.0) - std::pow((sigma_sf/r),6.0) );
		e = std::pow((sigma_sf/r),6.0);
		e = 4.0*epsilon_sf*( e*e - e );
	}else {
		e = 0.0;
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
//	rho_si_out = integral_simpson(rho_si_int_j, nstep, dr);
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

double ni(double *rho, double *r, int i, double *n0_j, double *n1_j, double *n2_j, double *n3_j, double *nv1_j, double *nv2_j,
		  double *n0, double *n1, double *n2, double *n3, double *nv1, double *nv2){
	int j,k;
	double ra;
	double raj;
	double xf, xs, xf2, xs2;
	double Ri;
	Ri = 0.2217; // Ri is the hard-sphere radius
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
	//nv1 = integral rho(r)/(4.0*M_PI)*(Rvi/Ri) dr, r=Ri, Rvi = vector Ri
	//nv1 = integral rho(z)/(4.0*M_PI)*(z/Ri)*(2.0*M_PI*x) dxdr, z<Ri, x^2+z^2=Ri^2, Rvi = vector Ri
	//nv2 = integral rho(r)*(Rvi/Ri) dr, r=Ri, Rvi = vector Ri
	//nv2 = integral rho(z)*(z/Ri)*(2.0*M_PI*x) dxdr, z<Ri, x^2+z^2=Ri^2, Rvi = vector Ri
	//
	//n0[i] = 0.0;
	//n1[i] = 0.0;
	//n2[i] = 0.0;
	//n3[i] = 0.0;
	//nv1[i] = 0.0;
	//nv2[i] = 0.0;
	for (j=0; j<nstep; j++) {
		raj = (r[j]-r[i]);
		xf2 = (d_hs*d_hs-raj*raj);
		if ( xf2 >= 0.0 ){
			xf = std::sqrt(xf2);
		} else{
			xf = 0.0;
		}
		xs2 = (Ri*Ri-raj*raj);
		if ( xs2 >= 0.0 ){
			xs = std::sqrt(xs2);
		} else{
			xs = 0.0;
		}
		//n0_j[j] = (rho[j])/(4.0*M_PI*Ri*Ri)*(2.0*M_PI*x);
		//n0_j[j] = (rho[j])/(2.0*Ri*Ri)*x;
		n0_j[j] = (rho[j])/(2.0*d_hs*d_hs)*xf + (rho_ssq(r[j])+rho_ssq(H-r[j]))/(2.0*Ri*Ri)*xs;
		//n1_j[j] = (rho[j])/(4.0*M_PI*Ri)*(2.0*M_PI*x);
		//n1_j[j] = (rho[j])/(2.0*Ri)*x;
		n1_j[j] = (rho[j])/(2.0*d_hs)*xf + (rho_ssq(r[j])+rho_ssq(H-r[j]))/(2.0*Ri)*xs;
		n2_j[j] = (rho[j])*(2.0*M_PI*xf) + (rho_ssq(r[j])+rho_ssq(H-r[j]))*(2.0*M_PI*xs);
		n3_j[j] = (rho[j])*(M_PI*xf*xf) + (rho_ssq(r[j])+rho_ssq(H-r[j]))*(M_PI*xs*xs);
		//nv1_j[j] = (rho[j])/(4.0*M_PI)*(raj/Ri)*(2.0*M_PI*x);
		nv1_j[j] = (rho[j])/(2.0)*(raj/d_hs)*xf + (rho_ssq(r[j])+rho_ssq(H-r[j]))/(2.0)*(raj/Ri)*xs;
		nv2_j[j] = (rho[j])*(raj/d_hs)*(2.0*M_PI*xf) + (rho_ssq(r[j])+rho_ssq(H-r[j]))*(raj/Ri)*(2.0*M_PI*xs);
		
		//
		//std::cout << i << ", " << j << ", " << r[i] << ", " << r[j] << ", " << raj << ", " << x << std::endl;
		//std::cout << rho[j] << ", " << rho_ssq(r[j]) << ", " << rho_ssq(H-r[j]) << ", " << n0_j[j] << std::endl;
		//
		//n0[i] = n0[i] + n0_j[j]*dr;
		//n1[i] = n1[i] + n1_j[j]*dr;
		//n2[i] = n2[i] + n2_j[j]*dr;
		//n3[i] = n3[i] + n3_j[j]*dr;
		//nv1[i] = nv1[i] + nv1_j[j]*dr;
		//nv2[i] = nv2[i] + nv2_j[j]*dr;
	}
    //integral_trapezoidal(double *f, int n, double dx)
	n0[i] = integral_trapezoidal(n0_j, nstep, dr);
	n1[i] = integral_trapezoidal(n1_j, nstep, dr);
	n2[i] = integral_trapezoidal(n2_j, nstep, dr);
	n3[i] = integral_trapezoidal(n3_j, nstep, dr);
	nv1[i] = integral_trapezoidal(nv1_j, nstep, dr);
	nv2[i] = integral_trapezoidal(nv2_j, nstep, dr);
	//
	//integral_simpson(double *f, int n, double dx)
	//n0[i] = integral_simpson(n0_j, nstep, dr);
	//n1[i] = integral_simpson(n1_j, nstep, dr);
	//n2[i] = integral_simpson(n2_j, nstep, dr);
	//n3[i] = integral_simpson(n3_j, nstep, dr);
	//nv1[i] = integral_simpson(nv1_j, nstep, dr);
	//nv2[i] = integral_simpson(nv2_j, nstep, dr);
	//
	//std::cout << "i, r[i], j, r[j], raj, xs, n0[i], n1[i], n2[i], n3[i], nv1[i], nv2[i]" << std::endl;
	//std::cout << i << ", " << r[i] << ", " << j-1 << ", " << r[j-1] << ", " << raj << ", " << xs << ", " << n0[i] << ", " << n1[i] << ", " << n2[i] << ", " << n3[i] << ", " << nv1[i] << ", " << nv2[i] << ", " << std::endl;
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
	int j,k;
	double ra;
	double raj;
	double rak;
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
	double tpidrc = 2.0*M_PI*drc;
	double Ri;
	Ri = d_hs;
	//Ri = 0.2217; // [nm] Ri is the hard-sphere radius of carbon
	//
	// Memo
	// df(x)/dx = [d/dx1,...,d/dxn]t * [f1(x),...,fn(x)]
	// df(r)/dr = [d/dx,d/dy,d/dz]t * [fx(r),fy(r),fz(r)]
	// dr/dr = [d/dx,d/dy,d/dz]t * [x,y,z] = [dx/dx,dy/dy,dz/dz] = [1,1,1] = [I] = 1
	// r > 0: d|r|/dr = dr/dr = [I] = 1
	// r < 0: d|r|/dr = d(-r)/dr = -[I] = -1
	//
	dphi_per_n0 = 0.0;
	dphi_per_n1 = 0.0;
	dphi_per_n2 = 0.0;
	dphi_per_n3 = 0.0;
	dphi_per_nv1 = 0.0;
	dphi_per_nv2 = 0.0;
	for (j=0; j<nstep; j++) {
		raj = (r[j]-r[i]);
		x2 = (Ri*Ri-raj*raj);
		if ( x2 >= 0.0 ){
			x = std::sqrt(x2);
		} else {
			x = 0.0;
		}
		sxi = std::abs(nv2[j]/n2[j]);
		//std::cout << j << ", sxi = " << sxi << std::endl;
		//
		// dphi/dn0
		//dphi_per_n0[j] = -std::log(1.0-n3[j])/(4.0*M_PI*Ri*Ri)*(2.0*M_PI*x);
		dphi_per_n0_j[j] = -std::log(1.0-n3[j])/(2.0*Ri*Ri)*x;
		//
		// dphi/dn1
		//dphi_per_n1_j[j] = ( n2[j]/(1.0-n3[j]) )/(4.0*M_PI*Ri)*(2.0*M_PI*x);
		dphi_per_n1_j[j] = ( n2[j]/(1.0-n3[j]) )/(2.0*Ri)*x;
		//
		// dphi/dn2
		if ( nv2[j]/n2[j] >= 0.0 ){
			// nv2/n2 > 0  ->  nv2/n2 = sxi
			sign = 1.0;
		}else if ( nv2[j]/n2[j] < 0.0 ) {
			// nv2/n2 < 0  ->  -nv2/n2 = -sxi
			sign = -1.0;
		}
		dphi_per_n2_j[j] = ( n1[j]/(1.0-n3[j])
			+ 3.0*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j]))*(1.0-3.0*sxi*sxi+2.0*sxi*sxi*sxi*sign)
			+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
				* (1.0-6.0*sxi*(-nv2[j]/(n2[j]*n2[j]))+6.0*sxi*sxi*(-nv2[j]/(n2[j]*n2[j]))*sign)
		)*(2.0*M_PI*x);
		//
		// dphi/dn3
		dphi_per_n3_j[j] = ( n0[j]/(1.0-n3[j]) 
			+ (n1[j]*n2[j] - nv1[j]*nv2[j])/((1.0-n3[j])*(1.0-n3[j])) 
			+ 2.0*n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))*(1.0-3.0*sxi*sxi+2.0*sxi*sxi*sxi)
		)*(M_PI*x*x);
		//
		// dphi/dnv1
		//dphi_per_nv1_j[j] = ( -nv2[j]/(1.0-n3[j]) )/(4.0*M_PI)*(raj/Ri)*(2.0*M_PI*x);
		dphi_per_nv1_j[j] = ( -nv2[j]/(1.0-n3[j]) )/(2.0)*(raj/Ri)*x;
		//
		// dphi/dnv2
		dphi_per_nv2_j[j] = ( -nv1[j]/(1.0-n3[j])
			+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
				* (1.0-6.0*sxi*(1.0/n2[j])+6.0*sxi*sxi*(1.0/n2[j])*sign)
		)*(raj/Ri)*(2.0*M_PI*x);
		//
		//dphi_per_n0 = dphi_per_n0 + dphi_per_n0_j[j]*dr;
		//dphi_per_n1 = dphi_per_n1 + dphi_per_n1_j[j]*dr;
		//dphi_per_n2 = dphi_per_n2 + dphi_per_n2_j[j]*dr;
		//dphi_per_n3 = dphi_per_n3 + dphi_per_n3_j[j]*dr;
		//dphi_per_nv1 = dphi_per_nv1 + dphi_per_nv1_j[j]*dr;
		//dphi_per_nv2 = dphi_per_nv2 + dphi_per_nv2_j[j]*dr;
	}
    //integral_trapezoidal(double *f, int n, double dx)
	dphi_per_n0  = -integral_trapezoidal(dphi_per_n0_j, nstep, dr);
	dphi_per_n1  = -integral_trapezoidal(dphi_per_n1_j, nstep, dr);
	dphi_per_n2  = -integral_trapezoidal(dphi_per_n2_j, nstep, dr);
	dphi_per_n3  = -integral_trapezoidal(dphi_per_n3_j, nstep, dr);
	dphi_per_nv1 = -integral_trapezoidal(dphi_per_nv1_j, nstep, dr);
	dphi_per_nv2 = -integral_trapezoidal(dphi_per_nv2_j, nstep, dr);
	//
	//integral_simpson(double *f, int n, double dx)
	//dphi_per_n0 = -integral_simpson(dphi_per_n0_j, nstep, dr);
	//dphi_per_n1 = -integral_simpson(dphi_per_n1_j, nstep, dr);
	//dphi_per_n2 = -integral_simpson(dphi_per_n2_j, nstep, dr);
	//dphi_per_n3 = -integral_simpson(dphi_per_n3_j, nstep, dr);
	//dphi_per_nv1 = -integral_simpson(dphi_per_nv1_j, nstep, dr);
	//dphi_per_nv2 = -integral_simpson(dphi_per_nv2_j, nstep, dr);
	//
	//std::cout << "i, dphi_per_n0, dphi_per_n1, dphi_per_n2, dphi_per_n3, dphi_per_nv1, dphi_per_nv2" << std::endl;
	//std::cout << i << ", " << dphi_per_n0 << "," << dphi_per_n1 << "," << dphi_per_n2 << "," << dphi_per_n3 << "," << dphi_per_nv1 << "," << dphi_per_nv2 << "," << std::endl;
	//
	dfex_out = dphi_per_n0 + dphi_per_n1 + dphi_per_n2 + dphi_per_n3 + dphi_per_nv1 + dphi_per_nv2;
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
			alpha_int_j[j]  = integral_simpson(alpha_int_k, nrmesh, drc);
		}
		//integral_simpson(double *f, int n, double dx)
		//alpha_other_method  = alpha_other_method + integral_simpson(alpha_int_j, nstep, dr)*2.0*dr;
		alpha_other_method  = alpha_other_method + integral_simpson(alpha_int_j, nstep, dr);
	}
	alpha_other_method  = alpha_other_method * 2.0 * dr / (H-sigma_ss);
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
		phi_att_ff_int_ij[i*nstep+j] = integral_simpson(phi_ff_int_k, nrmesh, drc);
		}
	}
	return 0;
}

// solid-fluid
double phi_att_sf_int(double *r, double *phi_att_sf_int_i){
	int i,j,k;
	double ra;
	double raj;
	double rak;
	//double drc = rc/double(nrmesh-1);
	//dd = drc;
	int sfmesh = 2000;
	double dsf = H/(sfmesh-1);
	double phi_sf_int_k[nrmesh];
	double phi_sf_int_j[sfmesh];
	double tpidrc = 2.0*M_PI*drc;
	phi_sf_int_k[0] = 0.0;
	for (i=0; i<nstep; i++) {
		for (j=0; j<sfmesh; j++) {
			//raj = (r[i]-r[j]);
			raj = (r[i]-double(j)*dsf);
			for (k=1; k<nrmesh; k++) {
				rak = drc*double(k);
				//ra = std::pow((r[i]-r[j]),2.0) + std::pow((double(k)*drc),2.0);
				//ra = (r[i]-r[j])*(r[i]-r[j]) + (double(k)*drc)*(double(k)*drc);
				ra = raj*raj + rak*rak;
				//ra = std::pow(ra,0.5);
				ra = std::sqrt(ra);
				//std::cout << ra << std::endl;
				//rho_phi_sf_int_k[k]  = rho_ss[j]*phi_att_sf(ra)*(2.0*M_PI*(double(k)*drc)); // old ver.1.1.0
				//phi_sf_int_k[k]  = phi_att_sf(ra)*(2.0*M_PI*(double(k)*drc));
				phi_sf_int_k[k]  = phi_att_sf(ra)*(tpidrc*double(k));
			}
			phi_sf_int_j[j] = (rho_ssq(double(j)*dsf)+rho_ssq(H-double(j)*dsf))*integral_simpson(phi_sf_int_k, nrmesh, drc);
		}
		phi_att_sf_int_i[i] = integral_simpson(phi_sf_int_j, sfmesh, dsf);
		//std::cout << rho_phi_sf_int << std::endl;
	}
	return 0;
}

// xi include kb1*T*(std::log(rho_b)) type.
// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
//double xi(double *rho, double *r, int i, double rho_b, double *rho_sj, double *rho_s0j, double *rho_s1j, double *rho_s2j, double *phi_att_ff_int_ij){
double xi(double *rho, double *r, int i, double rho_b, double *phi_att_ff_int_ij){
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
	double rho_phi_ff_int_j[nstep];
	double rho_phi_sf_int_j[nstep];
	double rho_dfex_int_k[nrmesh];
	double rho_phi_ff_int = 0.0;
	double rho_phi_sf_int = 0.0;
	double dz = H/(nstep-1);
	for (j=0; j<nstep; j++) {
		rho_phi_ff_int_j[j]  = rho[j]*phi_att_ff_int_ij[i*nstep+j];
		//rho_phi_sf_int_j[j]  = (rho_ssq(double(j)*dz)+rho_ssq(H-double(j)*dz))*phi_att_sf_int_ij[i*nstep+j];
		//std::cout << rho_phi_sf_int_j[j] << std::endl;
		//rho_phi_sf_int = rho_phi_sf_int + rho_phi_sf_int_j[j]*dr;
	}
	//integral_simpson(double *f, int n, double dx)
	rho_phi_ff_int  = integral_simpson(rho_phi_ff_int_j, nstep, dr);
	//rho_phi_sf_int  = integral_simpson(rho_phi_sf_int_j, nstep, dr);
	//std::cout << rho_phi_sf_int << std::endl;
	//
	double xi_out;
	//xi_out = kb1*T*std::log(rho_b) + mu_ex(rho_b) - rho_b*alpha - phi_ext(r[i]) - f_ex(rho_sj[i]) - rho_dfex_int - rho_phi_int; // old ver.1.1.1
	//xi_out = ( - rho_b*alpha - rho_dfex_int - f_ex(rho_sj[i]) ) + ( mu_ex(rho_b) - rho_phi_int ) + ( kb1*T*std::log(rho_b) - phi_ext(r[i]) );
	xi_out = ( - rho_b*alpha ) + ( mu_ex(rho_b) - rho_phi_ff_int ) + ( kb1*T*std::log(rho_b) );
	// debug
	//std::cout << i << ", " << xi_out << ", " << kb1*T*std::log(rho_b) << ", " << mu_ex(rho_b) << ", " << -rho_b*alpha << ", " << -rho_phi_ff_int << std::endl;
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

int main(){
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
	for (i=0; i<nstep; i++){
		// r[i] = sigma_ss/2.0 + (H-sigma_ss)/double(nstep)*double(i);
		// 1.72 times is escape nan, etc from positive value of wall potential
		//r[i] = sigma_ss*1.74/2.0 + dr*double(i) + dr/2.0; // dr = (H-sigma_ss*1.74)/double(nstep+1);
		//r[i] = sigma_ss/2.0 + dr*double(i); // dr = (H-sigma_ss)/double(nstep+1);
		//r[i] = dr*double(i);
		r[i] = (ze+sigma_sf) + dr*double(i); // dr = (H-2.0*(ze+sigma_sf))/double(nstep-1);
		//std::cout << i << ", " << r[i] << std::endl;
	}
	
	// show alpha
	calc_alpha(r);
	// alpha = calc_alpha(r);
	
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
	ofsppov << "# w = (H-2.0*(ze-sigma_sf)) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov << "# P/P0, V[molecules/nm^3]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "w = (H-2.0*(ze-sigma_sf)) = pore width = " << w_pw << " [nm]" << std::endl;
	std::cout << "P/P0, V[molecules/nm^3]" << std::endl;
	double rho_sj[nstep];
	double rho_s0j[nstep];
	double rho_s1j[nstep];
	double rho_s2j[nstep];
	double phi_att_ff_int_ij[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	double phi_att_sf_int_i[nstep];
	//
	double n0_j[nstep], n0[nstep];
	double n1_j[nstep], n1[nstep];
	double n2_j[nstep], n2[nstep];
	double n3_j[nstep], n3[nstep];
	double nv1_j[nstep], nv1[nstep];
	double nv2_j[nstep], nv2[nstep];
	double c1;
	phi_att_ff_int(r, phi_att_ff_int_ij); // calculate integral phi_att_ff at r[i]
	phi_att_sf_int(r, phi_att_sf_int_i); // calculate integral phi_att_sf at r[i]
	for (k=0; k<100; k++){
		rho_b = rho_b0 * std::exp(-(20.0-2.0*double(k+1.0)/10.0));
		//rho_b = rho_b0 * std::exp(-(20.0-2.0*double(99.0-k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		double check_data;
		for (j=0; j<cycle_max; j++){
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dr*nstep/2. 
			//rho_s(rho, r, rho_sj, rho_s0j, rho_s1j, rho_s2j); // for NLDFT
			for (i=0; i<nstep; i++){
				ni(rho, r, i, n0_j, n1_j, n2_j, n3_j, nv1_j, nv2_j, n0, n1, n2, n3, nv1, nv2);
			}
			for (i=0; i<=(nstep-2)/2; i++){
				c1 = dfex(r, i, n0, n1, n2, n3, nv1, nv2);
				//std::cout << "c1*(kb1*T) = " << c1*(kb1*T) << std::endl;
				//std::cout << "-phi_att_sf_i[i] = " << -phi_att_sf_int_i[i] << std::endl;
				//rho_new[i] = rho_b*std::exp(xi(rho,r[i],rho_b,r)/(kb1*T)); // this equation occure inf.
				//rho_new[i] = std::exp(c1+xi(rho,r,i,rho_b, rho_sj, rho_s0j, rho_s1j, rho_s2j, phi_att_ff_int_ij)/(kb1*T)-phi_att_sf_int_i[i]/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				rho_new[i] = std::exp(c1+xi(rho,r,i,rho_b, phi_att_ff_int_ij)/(kb1*T)-phi_att_sf_int_i[i]/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				//check_data = c1*(kb1*T)+xi(rho,r,i,rho_b, phi_att_ff_int_ij)-phi_att_sf_int_i[i];
				//std::cout << i << ", check_data = " << check_data << std::endl;
				//std::cout << "num of cycle i, r[i], rho_new[i], rho[i]" << std::endl;
				//std::cout << i << ", " << r[i] << ", "<< rho_new[i] << ", " << rho[i] << std::endl;
				//std::cout << i << ", " << rho[i] << ", " << rho_sj[i] << ", " << rho_s0j[i] << ", " << rho_s1j[i] << ", " << rho_s2j[i] << std::endl;
			}
			diff = 0.0;
			for (i=0; i<=(nstep-2)/2; i++){
				rho[i] = wmixing*rho_new[i] + (1.0-wmixing)*rho[i];
				rho[nstep-i] = rho[i]; // The rest is filled with mirror symmetry. 
				diff = diff + 2.0*std::abs((rho_new[i]-rho[i])/rho[i]);
			}
			if ( (diff/nstep*100.0) < 5.0) {
				break;
			}
			//std::cout << "--------------------------------------------------" << std::endl;
			//std::cout << "cycle=" << j << ", diff=" << diff << ", rho[nstep/2]=" << rho[nstep/2] << std::endl;
		}
		//
		//v_gamma = 0.0;
		//for (i=0; i<=(nstep-2)/2; i++){
			//std::cout << i << ", " << r[i] << ", " << rho[i] << std::endl;
			//v_gamma = v_gamma + 2.0*rho[i]*dr;
		//}
		v_gamma = integral_simpson(rho, nstep, dr);
		//v_gamma = v_gamma/(H-sigma_ss) - rho_b; // for NLDFT
		v_gamma = v_gamma/(H-2.0*(ze+sigma_sf)) - rho_b;
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
	return 0;
}
