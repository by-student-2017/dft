#include <fstream>   // for file in and out
#include <iostream>  // for cout
#include <cmath>     // for log, exp
#include <sstream>   // for read parameters

//#include "maxwell_construction.h"

using namespace std;

//non-local smoothed density approximation：SDA
//non-local density functional theory（NLDFT)
//Reference: https://www.j-ad.org/adsorption_news/30_1.pdf

// slit -> cylinder
// x-y plane -> nstep
// z axis -> nhmesh
// H = 5*sigma_ff (cut off radius)

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
//double H;
double Rcc;
double Dcc;
double sigma_ss;
//#define nstep=1001;
//constexpr int nstep = 1001;
int nstep;
double w_pw;
double dr;
// ---------- ----------- ------------ ------------
// assume rho is same value in x-y plane.
// cylinder and normalization, because of cut off (rc).
//int nhmesh = 20; //rho_si and xi function
int nhmesh;
//int ndmesh = d_hs*nhmesh/rc
// nrmesh is on theta range.
int nrmesh;
//double drc = rc/double(nhmesh-1);
double drc;
//double dd = 2.0*d_hs/double(ndmesh-1);
//double dd;
double dh;
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
// rho_b0 is related with P0
double rho_b0;
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
	Dcc = num[0];
	Rcc = Dcc/2.0; // The radial coordinate of the adsorption centers [nm]
	// ---------- ----------- ------------ ------------
	sigma_ss = num[1]; // [nm]
	// ---------- ----------- ------------ ------------
	nstep = int(num[2]);
	if ( nstep == 0 ) {
		nstep = int((Dcc-sigma_ss)/0.02 + 0.5);
		if ( nstep%2 == 1 ){
			nstep = nstep + 1;
		}
		std::cout << "--------------------------------------------------" << std::endl;
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
	nhmesh = int(num[9]);
	if ( nhmesh == 0 ) {
		nhmesh = int(rc/0.08 + 0.5);
		if ( nhmesh%2 == 1 ){
			nhmesh = nhmesh + 1;
		}
		std::cout << "autoset nhmesh = " << nhmesh << std::endl;
	}
	// ---------- ----------- ------------ ------------
	epsilon_sf = num[10]; // [K]
	// ---------- ----------- ------------ ------------
	sigma_sf = num[11]; // [nm]
	// ---------- ----------- ------------ ------------
	delta = num[12]; // nm
	// ---------- ----------- ------------ ------------
	rho_ss = num[13]; // The surface number density of the adsorption centers [nm^-2]
	// ---------- ----------- ------------ ------------
	m = num[14]; // [kg]
	// ---------- ----------- ------------ ------------
	T = num[15]; // [K]
	if ( d_hs == 0.0 ) { d_hs = d_bh_calc(epsilon_ff, sigma_ff); }
	// ---------- ----------- ------------ ------------
	rho_b0 = num[16];
	// ---------- ----------- ------------ ------------
	nrmesh = num[17];
	//nrmesh = 180/9;
	// ---------- ----------- ------------ ------------
	
	w_pw = (Dcc-sigma_ss); // pore width [nm]
	//dr = (Dcc-sigma_ss)/2.0/double(nstep-1);
	dr = (Dcc-sigma_ss)/2.0/double(nstep);
	rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
	
	// ---------- ----------- ------------ ------------
	
	dh = 2.0*d_hs/double(nhmesh-1);
	
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
	int j,k,t;
	double ra;
	double raj;
	double rak;
	double rho_si_out;
	double rho_si_int_j[nstep];
	double rho_si_int_k[nhmesh];
	double rho_si_int_t[nrmesh];
	double x,y;
	double drad = M_PI/nrmesh;
	for (j=0; j<nstep; j++) {
		for (k=0; k<nhmesh; k++) {
			rak = dh*double(k);
			//rho_si_int_k[k] = 0.0;
			for (t=0; t<nrmesh; t++) {
				x = r[j]*std::cos(drad*double(t));
				y = r[j]*std::sin(drad*double(t));
				raj = x - r1;
				ra = raj*raj + y*y + rak*rak;
				ra = std::sqrt(ra);
				rho_si_int_t[t] = wi(ra,i);
				//rho_si_int_k[k] += wi(ra,i);
			}
			//integral_simpson(double *f, int n, double dx)
			rho_si_int_k[k] = integral_simpson(rho_si_int_t, nrmesh-1, drad);
			//integral_trapezoidal(double *f, int n, double dx)
			//rho_si_int_k[k] = integral_trapezoidal(rho_si_int_t, nrmesh-1, drad);
		}
		//integral_simpson(double *f, int n, double dx)
		rho_si_int_j[j] = 2.0*r[j]*rho[j]*integral_simpson(rho_si_int_k, nhmesh-1, dh)*2.0;
		//rho_si_int_j[j] = 2.0*drad*r[j]*rho[j]*integral_simpson(rho_si_int_k, nhmesh-1, dh)*2.0;
		//integral_trapezoidal(double *f, int n, double dx)
		//rho_si_int_j[j] = 2.0*r[j]*rho[j]*integral_trapezoidal(rho_si_int_k, nhmesh-1, dh)*2.0;
		//rho_si_int_j[j] = 2.0*drad*r[j]*rho[j]*integral_trapezoidal(rho_si_int_k, nhmesh-1, dh)*2.0;
	}
	//integral_simpson(double *f, int n, double dx)
	rho_si_out = integral_simpson(rho_si_int_j, nstep-1, dr) + rho_si_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	//integral_trapezoidal(double *f, int n, double dx)
	//rho_si_out = integral_trapezoidal(rho_si_int_j, nstep-1, dr) + rho_si_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
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
//double phi_sf(double z){
//	double phi_sf_out;
//	double sigma_sf2 = sigma_sf*sigma_sf;
//	double sfpz = (sigma_sf/z);
//	double sfpz2 = sfpz*sfpz;
//	double dez = (0.61*delta+z);
//	//phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*std::pow(sigma_sf,2.0)*delta*
//	//			( (2.0/5.0)*std::pow((sigma_sf/z),10.0)-std::pow((sigma_sf/z),4.0)-std::pow(sigma_sf,4.0)/
//	//			(3.0*delta*std::pow((0.61*delta+z),3.0)) );
//	phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*(sigma_sf2)*delta*
//				( (2.0/5.0)*std::pow(sfpz2,5.0)-(sfpz2*sfpz2)-(sigma_sf2*sigma_sf2)/
//				(3.0*delta*(dez*dez*dez)) );
//	return phi_sf_out;
//}

// e.g., wall potential (Carbon slit)
//double phi_ext(double z){
//	double phi_ext_out;
//	phi_ext_out = phi_sf(z) + phi_sf(H-z);
//	//std::cout << phi_ext_out << std::endl;
//	return phi_ext_out;
//}

// if either a or b is a nonpositive integer
//double F(double a, double b, double c, double z){
//	int i;
//	double F_out;
//	double bi = 1.0;
//	double ci = 1.0;
//	double ni = 1.0; // factorial
//	double bc = 1.0;
//	F_out = 1.0;
//	for (i=1; i<-a; i++){
//		bc = bc*double((-a+1-i)/i); // binomial coefficients
//		bi = bi*(b+double(i-1));
//		ci = ci*(c+double(i-1));
//		F_out = F_out + std::pow(-1.0,i)*bc*(bi/ci)*std::pow(z,i);
//	}
//	return F_out;
//}

// the hypergeometric function
double Fh(double a, double b, double c, double z){
	int i;
	double Fh_tmp;
	double Fh_out;
	double ai = 1.0;
	double bi = 1.0;
	double ci = 1.0;
	double ni = 1.0; // factorial
	Fh_out = 1.0;
	for (i=1; i<5000; i++){
		ai = ai*(a+double(i-1));
		bi = bi*(b+double(i-1));
		ci = ci*(c+double(i-1));
		ni = ni*double(i);
		Fh_tmp = (ai*bi/ci)*(std::pow(z,i)/ni);
		//std::cout << "F[" << i << "] = " << Fh_tmp << std::endl;
		Fh_out = Fh_out + Fh_tmp;
		if ( std::abs(Fh_tmp) <= 1e-12 ){ break; }
	}
	return Fh_out;
}

// e.g., cylindrial layer
double phi_sf(double r){
	//double Dcc;
	//double Rcc = Dcc/2.0;
	double phi_ext_out;
	phi_ext_out = (M_PI*M_PI)*rho_ss*epsilon_sf*(sigma_sf*sigma_sf) * (
		(63.0/32.0)*std::pow(((Rcc-r)/sigma_sf)*(1.0+(r/Rcc)),-10.0) * Fh(-9.0/2.0,-9.0/2.0,1.0,(r/Rcc)*(r/Rcc))
		-3.0*std::pow(((Rcc-r)/sigma_sf)*(1.0+(r/Rcc)),-4.0) * Fh(-3.0/2.0,-3.0/2.0,1.0,(r/Rcc)*(r/Rcc)) );
	return phi_ext_out;
}

// e.g., wall potential (carbon cylinder)
double phi_ext(double r){
	//double Dcc;
	//double Rcc = Dcc/2.0;
	double phi_ext_out;
	//phi_ext_out = phi_sf(r) + phi_sf(Rcc-r);
	phi_ext_out = phi_sf(r);
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
	int i,j,k,t;
	double ra;
	double raj;
	double rak;
	double alpha_other_method;
	double alpha_int_i[nstep];
	double alpha_int_j[nstep];
	double alpha_int_k[nhmesh];
	double alpha_int_t[nrmesh];
	double x,y;
	double drad = M_PI/nrmesh;
	for (i=0; i<nstep; i++) {
		alpha_int_i[i] = 0.0;
		for (j=0; j<nstep; j++) {
			alpha_int_j[j] = 0.0;
			for (k=0; k<nhmesh; k++) {
				rak = dh*double(k);
				alpha_int_k[k] = 0.0;
				for (t=0; t<nrmesh; t++) {
					alpha_int_t[t] = 0.0;
					x = r[j]*std::cos(drad*double(t));
					y = r[j]*std::sin(drad*double(t));
					raj = x - r[i];
					ra = raj*raj + y*y + rak*rak;
					ra = std::sqrt(ra);
					alpha_int_t[t] = -phi_att(ra);
					//alpha_int_k[k] += -phi_att(ra);
				}
				//integral_simpson(double *f, int n, double dx)
				alpha_int_k[k] = integral_simpson(alpha_int_t, nrmesh-1, drad);
				//integral_trapezoidal(double *f, int n, double dx)
				//alpha_int_k[k] = integral_trapezoidal(alpha_int_t, nrmesh-1, drad);
			}
			//integral_simpson(double *f, int n, double dx)
			alpha_int_j[j]  = 2.0*r[j]*integral_simpson(alpha_int_k, nhmesh-1, dh)*2.0;
			//alpha_int_j[j]  = 2.0*drad*r[j]*integral_simpson(alpha_int_k, nhmesh-1, dh)*2.0;
			//integral_trapezoidal(double *f, int n, double dx)
			//alpha_int_j[j]  = 2.0*r[j]*integral_trapezoidal(alpha_int_k, nhmesh-1, dh)*2.0;
			//alpha_int_j[j]  = 2.0*drad*r[j]*integral_trapezoidal(alpha_int_k, nhmesh-1, dh)*2.0;
		}
		//integral_simpson(double *f, int n, double dx)
		alpha_int_i[i] = integral_simpson(alpha_int_j, nstep-1, dr) + alpha_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
		//integral_trapezoidal(double *f, int n, double dx)
		//alpha_int_i[i] = integral_trapezoidal(alpha_int_j, nstep-1, dr) + alpha_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	}
	//alpha_other_method = alpha_int_i[0]; // check
	//integral_simpson(double *f, int n, double dx)
	alpha_other_method = integral_simpson(alpha_int_i, nstep-1, dr) / ((Dcc-sigma_ss)/2.0);
	//integral_trapezoidal(double *f, int n, double dx)
	//alpha_other_method = integral_trapezoidal(alpha_int_i, nstep-1, dr) / ((Dcc-sigma_ss)/2.0);
	//
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "average alpha of other method = " << alpha_other_method << " in (carbon) cylinder" << std::endl;
	return alpha_other_method;
}

double phi_att_int(double *r, double *phi_att_int_ij){
	int i,j,k,t;
	double ra;
	double raj;
	double rak;
	double rho_phi_int_k[nhmesh];
	double rho_phi_int_t[nrmesh];
	double x,y;
	double drad = M_PI/nrmesh;
	for (i=0; i<nstep; i++) {
		for (j=0; j<nstep; j++) {
			for (k=0; k<nhmesh; k++) {
				rak = dh*double(k);
				//rho_phi_int_k[k] = 0.0;
				for (t=0; t<nrmesh; t++) {
					x = r[j]*std::cos(drad*double(t));
					y = r[j]*std::sin(drad*double(t));
					raj = x - r[i];
					ra = raj*raj + y*y + rak*rak;
					ra = std::sqrt(ra);
					rho_phi_int_t[t] = phi_att(ra);
					//rho_phi_int_k[k] += phi_att(ra);
				}
				//integral_simpson(double *f, int n, double dx)
				rho_phi_int_k[k] = integral_simpson(rho_phi_int_t, nrmesh-1, drad);
				//integral_trapezoidal(double *f, int n, double dx)
				//rho_phi_int_k[k] = integral_trapezoidal(rho_phi_int_t, nrmesh-1, drad);
			}
			//integral_simpson(double *f, int n, double dx)
			phi_att_int_ij[i*nstep+j] = integral_simpson(rho_phi_int_k, nhmesh-1, dh)*2.0;
			//integral_trapezoidal(double *f, int n, double dx)
			//phi_att_int_ij[i*nstep+j] = integral_trapezoidal(rho_phi_int_k, nhmesh-1, dh)*2.0;
		}
	}
	return 0;
}

// xi include kb1*T*(std::log(rho_b)) type.
// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
double xi(double *rho, double *r, int i, double rho_b, double *rho_sj, double *rho_s0j, double *rho_s1j, double *rho_s2j, double *phi_att_int_ij, double *rho_dfex_int, double *rho_phi_int, double *phi_ext_i){
	int j,k,t;
	double ra;
	double raj;
	double rak;
	double rho_dfex_int_j[nstep];
	double rho_phi_int_j[nstep];
	double rho_dfex_int_k[nhmesh];
	//double rho_phi_int_k[nhmesh];
	double rho_dfex_int_t[nrmesh];
	double x,y;
	double drad = M_PI/nrmesh;
	for (j=0; j<nstep; j++) {
		for (k=0; k<nhmesh; k++) {
			rak = dh*double(k);
			//rho_dfex_int_k[k] = 0.0;
			//rho_phi_int_k[k] = 0.0;
			for (t=0; t<nrmesh; t++) {
				x = r[j]*std::cos(drad*double(t));
				y = r[j]*std::sin(drad*double(t));
				raj = x - r[i];
				ra = raj*raj + y*y + rak*rak;
				ra = std::sqrt(ra);
				rho_dfex_int_t[t] = drhos_per_drho_j(ra, rho_sj[j], rho_s1j[j], rho_s2j[j]);
				//rho_dfex_int_k[k] += drhos_per_drho_j(ra, rho_sj[j], rho_s1j[j], rho_s2j[j]);
				//rho_phi_int_k[k] += phi_att(ra);
			}
			//integral_simpson(double *f, int n, double dx)
			rho_dfex_int_k[k] = integral_simpson(rho_dfex_int_t, nrmesh-1, drad);
			//integral_trapezoidal(double *f, int n, double dx)
			//rho_dfex_int_k[k] = integral_trapezoidal(rho_dfex_int_t, nrmesh-1, drad);
		}
		//integral_simpson(double *f, int n, double dx)
		rho_dfex_int_j[j] = 2.0*r[j]*rho[j]*dfex_per_drhos(rho_sj[j])*integral_simpson(rho_dfex_int_k, nhmesh-1, dh)*2.0;
		//rho_dfex_int_j[j] = 2.0*drad*r[j]*rho[j]*dfex_per_drhos(rho_sj[j])*integral_simpson(rho_dfex_int_k, nhmesh-1, dh)*2.0;
		//rho_phi_int_j[j]  = 2.0*drad*r[j]*rho[j]*integral_simpson(rho_phi_int_k, nhmesh-1, dh)*2.0;
		//integral_trapezoidal(double *f, int n, double dx)
		//rho_dfex_int_j[j] = 2.0*r[j]*rho[j]*dfex_per_drhos(rho_sj[j])*integral_trapezoidal(rho_dfex_int_k, nhmesh-1, dh)*2.0;
		//rho_dfex_int_j[j] = rho[j]*dfex_per_drhos(rho_sj[j])*integral_trapezoidal(rho_dfex_int_k, nhmesh-1, dh)*2.0;
		//rho_phi_int_j[j]  = rho[j]*integral_trapezoidal(rho_phi_int_k, nhmesh-1, dh)*2.0;
		// other method
		rho_phi_int_j[j]  = 2.0*r[j]*rho[j]*phi_att_int_ij[i*nstep+j];
		//rho_phi_int_j[j]  = 2.0*drad*r[j]*rho[j]*phi_att_int_ij[i*nstep+j];
	}
	//integral_simpson(double *f, int n, double dx)
	rho_dfex_int[i] = integral_simpson(rho_dfex_int_j, nstep-1, dr) + rho_dfex_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	rho_phi_int[i]  = integral_simpson(rho_phi_int_j, nstep-1, dr) + rho_phi_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	//integral_trapezoidal(double *f, int n, double dx)
	//rho_dfex_int[i] = integral_trapezoidal(rho_dfex_int_j, nstep-1, dr) + rho_dfex_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	//rho_phi_int[i]  = integral_trapezoidal(rho_phi_int_j, nstep-1, dr) + rho_phi_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	//
	double xi_out;
	//xi_out = kb1*T*std::log(rho_b) + mu_ex(rho_b) - rho_b*alpha - phi_ext(r[i]) - f_ex(rho_sj[i]) - rho_dfex_int - rho_phi_int; // old ver.1.1.1
	xi_out = ( - rho_b*alpha - rho_dfex_int[i] - f_ex(rho_sj[i]) ) + ( mu_ex(rho_b) - rho_phi_int[i] ) + ( kb1*T*std::log(rho_b) - phi_ext_i[i] );
	// debug
	//std::cout << "xi, (kb1*T)*log(rho_b), mu_ex(rho_b), -rho_b*alpha, -phi_ext(r[i]), -f_ex(rho_s(rho,r[i],r)), -rho_dfex_int, -rho_phi_int" << std::endl;
	//std::cout << xi_out << ", " << kb1*T*std::log(rho_b) << ", " << mu_ex(rho_b) << ", " << -rho_b*alpha << ", " << -phi_ext_i[i] << ", " << -f_ex(rho_sj[i]) << ", " << -rho_dfex_int[i] << ", " << -rho_phi_int[i] << std::endl;
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

// grand potential
double omega(double *rho, double *r, double *rho_dfex_int, double *rho_phi_int){
	double omega_out;
	double omega1, omega2, omega3;
	int i;
	double rhor[nstep];
	double rho_x_rho_dfex_int[nstep];
	double rho_x_rho_phi_int[nstep];
	double tpidr = 2.0*M_PI*dr;
	for (i=0; i<nstep; i++){
		rhor[i] = rho[i] * tpidr*r[i];
		rho_x_rho_dfex_int[i] = rho[i] * rho_dfex_int[i] * tpidr*r[i];
		rho_x_rho_phi_int[i]  = rho[i] * rho_phi_int[i] * tpidr*r[i];
	}
	//integral_simpson(double *f, int n, double dx)
	omega1 = -(kb1*T) * integral_simpson(rhor, nstep-1, dr) + rhor[0]/(tpidr*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	omega2 = -integral_simpson(rho_x_rho_dfex_int, nstep-1, dr) + rho_x_rho_dfex_int[0]/(tpidr*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	omega3 = -0.5 * integral_simpson(rho_x_rho_phi_int, nstep-1, dr) + rho_x_rho_phi_int[0]/(tpidr*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	//integral_trapezoidal(double *f, int n, double dx)
	//omega1 = -(kb1*T) * integral_trapezoidal(rho, nstep-1, dr) + rhor[0]/(tpidr*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	//omega2 = -integral_trapezoidal(rho_x_rho_dfex_int, nstep-1, dr) + rho_x_rho_dfex_int[0]/(tpidr*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	//omega3 = -0.5 * integral_trapezoidal(rho_x_rho_phi_int, nstep-1, dr) + rho_x_rho_phi_int[0]/(tpidr*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	omega_out = (omega1 + omega2 + omega3) / epsilon_ff;
	return omega_out;
}

int main(){
	int i,j,k;
	double diff;
	double v_gamma;
	double press_b, press_b0, pp0;
	double rho_b;
	double v_mmol_per_cm3;
	double v_cm3STP_per_g;
	double grand_potential;
	//
	//double Dref=(Dcc-sigma_ss);
	//
	read_parameters();
	double r[nstep];
	double rho[nstep], rho_new[nstep];
	//
	for (i=0; i<nstep; i++){
		r[i] = dr*(0.5+double(i)); // dr = (Dcc-sigma_ss)/double(nstep+1);
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
	ofsppov_vs << "# w = (2.0*Rcc-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_vs << "# P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "w = (2.0*Rcc-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	std::cout << "P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
	double rho_sj[nstep];
	double rho_s0j[nstep];
	double rho_s1j[nstep];
	double rho_s2j[nstep];
	double rho_dfex_int[nstep];
	double rho_phi_int[nstep];
	double phi_ext_i[nstep];
	for (i=0; i<nstep; i++){
		phi_ext_i[i] = phi_ext(r[i]);
		//std::cout << "phi_ext_i[" << i << "] = " << phi_ext_i[i] << std::endl;
	}
	double *phi_att_int_ij = (double *)malloc(sizeof(double)*((nstep+1)*nstep));
	if (phi_att_int_ij == NULL) {
		printf("Memory cannot be allocated.");
		std::exit(1);
	} else {
		printf("Memory has been allocated. The address is %p\n", phi_att_int_ij);
	}
	phi_att_int(r, phi_att_int_ij); // calculate integral phi_att at r[i]
	double diff0;
	double mixing;
	//
	double rho_r[nstep];
	for (k=0; k<100; k++){
		rho_b = rho_b0 * std::exp(-(20.0-2.0*double(k+1.0)/10.0));
		//rho_b = rho_b0 * std::exp(-(20.0-2.0*double(99.0-k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		for (j=0; j<cycle_max; j++){
			rho_s(rho, r, rho_sj, rho_s0j, rho_s1j, rho_s2j);
			for (i=0; i<nstep; i++){
				//rho_new[i] = rho_b*std::exp(xi(rho,r[i],rho_b,r)/(kb1*T)); // this equation occure inf.
				rho_new[i] = std::exp(xi(rho,r,i,rho_b, rho_sj, rho_s0j, rho_s1j, rho_s2j, phi_att_int_ij, rho_dfex_int, rho_phi_int, phi_ext_i)/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				//std::cout << "num of cycle i, r[i], rho_new[i], rho[i]" << std::endl;
				//std::cout << i << ", " << r[i] << ", "<< rho_new[i] << ", " << rho[i] << std::endl;
				//std::cout << i << ", " << rho[i] << ", " << rho_sj[i] << ", " << rho_s0j[i] << ", " << rho_s1j[i] << ", " << rho_s2j[i] << std::endl;
				//
				// overflow about std::exp(730)
				// to avoid overflow
				if (rho_new[i] > 1e9){
					rho_new[i] = 1e9;
				}
				// to avoid -inf or int
				if (rho_new[i] < 1e-18 && rho[i] < 1e-18){
					rho_new[i] = 1e-18;
					rho[i] = 1e-18;
				}
			}
			diff = 0.0;
			for (i=0; i<nstep; i++){
				diff0 = std::abs((rho_new[i]-rho[i])/rho[i]);
				diff = diff + 2.0*diff0;
				mixing = wmixing + wmixing/(0.5+diff0);
				//std::cout << i << ", " << mixing << std::endl;
				rho[i] = mixing*rho_new[i] + (1.0-mixing)*rho[i];
			}
			if ( (diff/nstep*100.0) < 5.0 && j >= 100) {
				break;
			}
		}
		//for (i=0; i<nstep; i++){
		//	std::cout << "--------------------------------------------------" << std::endl;
		//	std::cout << "cycle=" << j << ", r[" << i << "]" << r[i] << " , xi " << xi(rho,r,i,rho_b, rho_sj, rho_s0j, rho_s1j, rho_s2j, phi_att_int_ij, rho_dfex_int, rho_phi_int, phi_ext_i)/(kb1*T) << ", rho[" << i << "]=" << rho[i] << std::endl;
		//}
		
		for (i=0; i<nstep; i++){
			rho_r[i] = rho[i]*2.0*M_PI*r[i];
		}
		//integral_simpson(double *f, int n, double dx)
		v_gamma = integral_simpson(rho_r, nstep-1, dr) + rho[0]*M_PI*(dr/2.0)*(dr/2.0);
		//integral_trapezoidal(double *f, int n, double dx)
		//v_gamma = integral_trapezoidal(rho_r, nstep-1, dr) + rho[0]*M_PI*(dr/2.0)*(dr/2.0);
		v_gamma = v_gamma/(M_PI*(Dcc-sigma_ss)/2.0*(Dcc-sigma_ss)/2.0) - rho_b;
		//v_gamma = 2.0*v_gamma/(Dcc-sigma_ss) - rho_b*(Dref*Dref)/(4.0*(Dcc-sigma_ss));
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10.0; // [mmol/cm3]
		//v_cm3STP_per_g = v_mmol_per_cm3 / 22.414 / (rho_ss*12.0107*(1e7*1e7*1e7)/(6.02214076*1e23)); // [cm3(STP)/g], 2.226 [g/cm3]
		v_cm3STP_per_g = v_mmol_per_cm3 / 22.414 / (rho_ss*12.0107*10.0/6.02214076); // [cm3(STP)/g], 2.226 [g/cm3]
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
		grand_potential = omega(rho, r, rho_dfex_int, rho_phi_int);
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_vs << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_g << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_g << ", " << grand_potential << std::endl;
	}
	// reverse
	std::ofstream ofsppov_ls("./PP0_vs_Vgamma_data_ls.txt");
	ofsppov_ls << "# w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_ls << "# P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	//std::cout << "w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	//std::cout << "P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
	for (k=0; k<100; k++){
		//rho_b = rho_b0 * std::exp(-(20.0-2.0*double(k+1.0)/10.0));
		rho_b = rho_b0 * std::exp(-(20.0-2.0*double(99.0-k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		for (j=0; j<cycle_max; j++){
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dr*nstep/2. 
			rho_s(rho, r, rho_sj, rho_s0j, rho_s1j, rho_s2j);
			for (i=0; i<nstep; i++){
				//rho_new[i] = rho_b*std::exp(xi(rho,r[i],rho_b,r)/(kb1*T)); // this equation occure inf.
				rho_new[i] = std::exp(xi(rho,r,i,rho_b, rho_sj, rho_s0j, rho_s1j, rho_s2j, phi_att_int_ij, rho_dfex_int, rho_phi_int, phi_ext_i)/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				//std::cout << "num of cycle i, r[i], rho_new[i], rho[i]" << std::endl;
				//std::cout << i << ", " << r[i] << ", "<< rho_new[i] << ", " << rho[i] << std::endl;
				//
				// overflow about std::exp(730)
				// to avoid overflow
				if (rho_new[i] > 1e9){
					rho_new[i] = 1e9;
				}
				// to avoid -inf or int
				if (rho_new[i] < 1e-18 && rho[i] < 1e-18){
					rho_new[i] = 1e-18;
					rho[i] = 1e-18;
				}
			}
			diff = 0.0;
			for (i=0; i<nstep; i++){
				diff0 = std::abs((rho_new[i]-rho[i])/rho[i]);
				diff = diff + 2.0*diff0;
				mixing = wmixing + wmixing/(0.5+diff0);
				//std::cout << i << ", " << mixing << std::endl;
				rho[i] = mixing*rho_new[i] + (1.0-mixing)*rho[i];
			}
			if ( (diff/nstep*100.0) < 5.0 && j >= 100) {
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
		for (i=0; i<nstep; i++){
			rho_r[i] = rho[i]*2.0*M_PI*r[i];
		}
		//integral_simpson(double *f, int n, double dx)
		v_gamma = integral_simpson(rho_r, nstep-1, dr) + rho[0]*M_PI*(dr/2.0)*(dr/2.0);
		//integral_trapezoidal(double *f, int n, double dx)
		//v_gamma = integral_trapezoidal(rho_r, nstep-1, dr) + rho[0]*M_PI*(dr/2.0)*(dr/2.0);
		v_gamma = v_gamma/(M_PI*(Dcc-sigma_ss)/2.0*(Dcc-sigma_ss)/2.0) - rho_b;
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10.0; // [mmol/cm3]
		//v_cm3STP_per_g = v_mmol_per_cm3 / 22.414 / (rho_ss*12.0107*(1e7*1e7*1e7)/(6.02214076*1e23)); // [cm3(STP)/g], 2.226 [g/cm3]
		v_cm3STP_per_g = v_mmol_per_cm3 / 22.414 / (rho_ss*12.0107*10.0/6.02214076); // [cm3(STP)/g], 2.226 [g/cm3]
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
		grand_potential = omega(rho, r, rho_dfex_int, rho_phi_int);
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_ls << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_g << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_g << ", " << grand_potential << std::endl;
	}
	return 0;
}
