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
// Units are fundamentaly [K] and [nm] in this routine.

// There are many imperfections, so I hope someone can make it better with a CC0 license. 
// It seems that this code is the first in the world at present (2021/7/5) to be released on CC0 even in NLDFT. 

// compiling: c++ nldft_fmt_cylinder.cpp -O2
// usage: ./a.out

// debag mode
// compiling: c++ nldft_fmt_cylinder.cpp -g -Wall -O0
// run: gdb ./a.out
//      (gdb) run

// ---------- ----------- ------------ ------------
// Adsorbent 
//float H = 1.00; //distace of slit [nm]
//float sigma_ss = 0.34; // [nm]
//int nstep = 100;
//float w_pw = (H-sigma_ss); // pore width [nm]
//float dr = w_pw/float(nstep);
float H;
float Rcc;
float Dcc;
float sigma_ss;
//#define nstep=1001;
//constexpr int nstep = 1001;
int nstep;
float w_pw;
float dr;
// ---------- ----------- ------------ ------------
// assume rho is same value in x-y plane.
// cylinder and normalization, because of cut off (rc).
//int nrmesh = 20; //rho_si and xi function
int nhmesh;
//int ndmesh = d_hs*nhmesh/rc
// nrmesh is on theta range.
int nrmesh;
//int ndmesh = d_hs*nrmesh/rc
int ndmesh;
//float drc = rc/float(nrmesh-1);
float drc;
//float dd = 2.0*d_hs/float(ndmesh-1);
//float dd;
float dh;
// ---------- ----------- ------------ ------------
// iteration of rho
//int cycle_max = 50;
int cycle_max;
//float wmixing = 0.005;
float wmixing;
// ---------- ----------- ------------ ------------
//Carbon dioxide 253.9  [K](epsilon), 0.3454 [nm](sigma), 0.3495 [nm](d_hs)
//Argon          118.05 [K](epsilon), 0.3305 [nm](sigma), 0.3390 [nm](d_hs)
//Nitrogen        94.45 [K](epsilon), 0.3575 [nm](sigma), 0.3575 [nm](d_hs)
//extern float epsilon_ff = 94.45;
//float sigma_ff = 0.3575;
//extern float d_hs = 0.3575; // Maxwell_construction()
//float rc = 1.28; // [nm],cut off, (12.8 [A])
float epsilon_ff;
float sigma_ff;
float d_hs;
float rc;
// ---------- ----------- ------------ ------------
//float rm = std::pow(2.0,1.0/6.0)*sigma_ff; //minimum position of LJ
//float rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
float rm;
// ---------- ----------- ------------ ------------
// Carbon dioxide/Carbon slit 81.5  [K](epsilon), 0.3430 [nm](sigma)
// Nitrogen/Carbon slit       53.72 [K](epsilon), 0.3508 [nm](sigma)
//float epsilon_sf = 53.72; // [K] 
//float sigma_sf = 0.3508; // [nm]
float epsilon_sf;
float sigma_sf;
// ---------- ----------- ------------ ------------
// slit pore (graphite)
//float delta = 0.335; // [nm]
//float rho_ss = 114.0; // [nm^-3], [molecules/nm3]?, 0.114 [A^-3]
float delta;
float rho_ss;
// ---------- ----------- ------------ ------------
//float m = 14.0067*2.0/(6.02214076e23)/1000; // N2 = 4.65173e-26 [kg]
//float m = 4.65173e-26; //[kg] (N2) (e.g., Ar = 6.63e-26 [kg])
float m;
float kb1 = 1.0;
float kb = 1.38e-23; //[J/K] (8.61733262e-5 [eV/K])
//extern float T = 77.347; //[K]
float T;
float h = 6.63e-34; //[Js] (4.135667696e-15 [eVs])
// thermal de Broglie wavelength
//extern float lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm], Maxwell_construction()
float lam;
// Ref: https://www1.doshisha.ac.jp/~bukka/lecture/statistic/pdftext/std-07.pdf
// ---------- ----------- ------------ ------------
// alpha = integal phi_att * -1.0
//extern float alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
//	( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
float alpha;
// ---------- ----------- ------------ ------------
// rho_b0 is related with P0
float rho_b0;
// ---------- ----------- ------------ ------------
// P0
float p0;
// ---------- ----------- ------------ ------------

float integral_trapezoidal(float *f, int n, float dx){
	float sum;
	sum = 0.0;
	int i;
	for(i=1; i<n; i++){
		sum += (f[i-1]+f[i])/2.0*dx;
	}
	return sum;
}

float integral_simpson(float *f, int n, float dx){
	//if( (n+1)%2 == 1 ){
	//	std::cout << "Error, plase change number of data to even ( = array[odd] )" << std::endl;
	//}
	float sum;
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
float d_bh_calc_v1(float epsilon, float sigma){
	//float epsilon = 94.45;
	//float sigma = 0.3575;
	//Lstoskie et al.,
	float xi1 = 0.3837;
	float xi2 = 1.035;
	float xi3 = 0.4249;
	float xi4 = 1.0;
	float d_bh_out;
	d_bh_out = (xi1*kb1*T/epsilon+xi2)/(xi3*kb1*T/epsilon+xi4)*sigma;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "d = d_hs = " << d_bh_out << " [nm] at " << T << " [K] from Barker-Henderson (BH) theory" << std::endl;
	return d_bh_out;
}

//Barker-Henderson (BH) perturbation theory
float d_bh_calc_v2(float epsilon, float sigma){
	//float epsilon = 94.45;
	//float sigma = 0.3575;
	//Lstoskie et al.,
	float d_bh_out;
	float Ts = kb1*T/epsilon;
	d_bh_out = (1.0+0.2977*Ts)/(1.0+0.331637*Ts+0.00104771*Ts*Ts)*sigma;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "d = d_hs = " << d_bh_out << " [nm] at " << T << " [K] from Barker-Henderson (BH) perturbation theory" << std::endl;
	return d_bh_out;
}

void read_parameters(void){
	std::ifstream ifs("parameters.txt");
	std::string str;
	float num[20];
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
		nstep = int(((Dcc-sigma_ss)/2.0)/0.005 + 0.5);
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
	//if ( d_hs == 0.0 ) { d_hs = d_bh_calc_v1(epsilon_ff, sigma_ff); }
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
		nhmesh = int(rc/0.02 + 0.5);
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
	rho_ss = num[13]; // [nm^-3], [mulecules/nm3]
	// ---------- ----------- ------------ ------------
	m = num[14]; //[g/mol], H2=2.01568, Ar=39.948, N2=28.0134, CO2=44.01, O2=31.998
	m = m/(6.02214076e23)/1000; //[kg]
	// ---------- ----------- ------------ ------------
	T = num[15]; // [K]
	if ( d_hs == 0.0 ) { d_hs = d_bh_calc_v1(epsilon_ff, sigma_ff); }
	// ---------- ----------- ------------ ------------
	rho_b0 = num[16];
	// ---------- ----------- ------------ ------------
	nrmesh = num[17];
	//nrmesh = 180/9;
	if ( nrmesh == 0 ) {
		nrmesh = int(((Dcc-sigma_ss)/2.0)/0.01 + 0.5);
		if ( nrmesh%2 == 1 ){
			nrmesh = nrmesh + 1;
		}
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "autoset nrmesh = " << nrmesh << std::endl;
	}
	// ---------- ----------- ------------ ------------
	p0 = num[18];
	// ---------- ----------- ------------ ------------
	
	w_pw = (Dcc-sigma_ss); // pore width [nm]
	dr = ((Dcc-sigma_ss)/2.0)/float(nstep+1);
	rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
	
	// ---------- ----------- ------------ ------------
	
	dh = 2.0*d_hs/float(nhmesh-1);
	
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

float phi_att(float r){
	float e;
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
float wi(float r, int i){
	float wi_out;
	float rpdhs;
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

// Tarazona theory, for SDA
//float rho_si(float *rho, float r1, float *r, int i){
//	int j,k;
//	float ra;
//	float raj;
//	float rak;
//	//float ndmesh = 2*d_hs*nrmesh/rc;
//	//float dd = 2.0*d_hs/float(ndmesh-1);
//	//dd = drc;
//	float tpidd = 2.0*M_PI*dd;
//	float rho_si_out;
//	float rho_si_int_j[nstep];
//	float rho_si_int_k[nrmesh];
//	rho_si_int_k[0] = 0.0;
//	for (j=0; j<nstep; j++) {
//		raj = (r1-r[j]);
//		for (k=1; k<ndmesh; k++) {
//			rak = dd*float(k);
//			//ra = std::pow((r1-r[j]),2.0) + std::pow((float(k)*dd),2.0);
//			//ra = (r1-r[j])*(r1-r[j]) + (float(k)*dd)*(float(k)*dd);
//			ra = raj*raj + rak*rak;
//			//ra = std::pow(ra,0.5);
//			ra = std::sqrt(ra);
//			//std::cout << ra << std::endl;
//			//
//			//rho_si_int_k[k] = rho[j]*wi(ra,i)*(2.0*M_PI*(float(k)*dd)); // old ver.1.1.0
//			//rho_si_int_k[k] = wi(ra,i)*(2.0*M_PI*(float(k)*dd));
//			rho_si_int_k[k] = wi(ra,i)*(tpidd*float(k));
//		}
//		//integral_simpson(float *f, int n, float dx)
//		//rho_si_int_j[j] = integral_simpson(rho_si_int_k, ndmesh, dd); // old ver.1.1.0
//		rho_si_int_j[j] = rho[j]*integral_simpson(rho_si_int_k, ndmesh, dd);
//	}
//	//integral_simpson(float *f, int n, float dx)
//	rho_si_out = integral_simpson(rho_si_int_j, nstep, dr);
//	//
//	return rho_si_out;
//}

// smoothed density approximation (SDA)
//float rho_s(float *rho, float r1, float *r){
//	float rho_den1, rho_den2, rho_s_out;
//	//rho_den1 = std::pow((1.0 - rho_si(rho,r1,r,1)),2.0);
//	rho_den1 = (1.0 - rho_si(rho,r1,r,1));
//	rho_den1 = rho_den1 * rho_den1;
//	//rho_den2 = std::pow((rho_den1 - 4.0*rho_si(rho,r1,r,0)*rho_si(rho,r1,r,2)),0.5);
//	rho_den2 = std::sqrt(rho_den1 - 4.0*rho_si(rho,r1,r,0)*rho_si(rho,r1,r,2));
//	rho_s_out = 2.0*rho_si(rho,r1,r,0)/(1.0 - rho_si(rho,r1,r,1)+rho_den2);
//	return rho_s_out;
//}

// smoothed density approximation (SDA), modified version
//float rho_s(float *rho, float *r, float *rho_sj, float *rho_s0j, float *rho_s1j, float *rho_s2j){
//	int j;
//	float rho_den1j, rho_den2j;
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

// Steele 10-4-3 potential
//float phi_sf(float z){
//	float phi_sf_out;
//	float sigma_sf2 = sigma_sf*sigma_sf;
//	float sfpz = (sigma_sf/z);
//	float sfpz2 = sfpz*sfpz;
//	float dez = (0.61*delta+z);
//	//phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*std::pow(sigma_sf,2.0)*delta*
//	//			( (2.0/5.0)*std::pow((sigma_sf/z),10.0)-std::pow((sigma_sf/z),4.0)-std::pow(sigma_sf,4.0)/
//	//			(3.0*delta*std::pow((0.61*delta+z),3.0)) );
//	phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*(sigma_sf2)*delta*
//				( (2.0/5.0)*std::pow(sfpz2,5.0)-(sfpz2*sfpz2)-(sigma_sf2*sigma_sf2)/
//				(3.0*delta*(dez*dez*dez)) );
//	return phi_sf_out;
//}

// e.g., wall potential (Carbon slit)
//float phi_ext(float z){
//	float phi_ext_out;
//	phi_ext_out = phi_sf(z) + phi_sf(H-z);
//	//std::cout << phi_ext_out << std::endl;
//	return phi_ext_out;
//}

// if either a or b is a nonpositive integer
//float F(float a, float b, float c, float z){
//	int i;
//	float F_out;
//	float bi = 1.0;
//	float ci = 1.0;
//	float ni = 1.0; // factorial
//	float bc = 1.0;
//	F_out = 1.0;
//	for (i=1; i<-a; i++){
//		bc = bc*float((-a+1-i)/i); // binomial coefficients
//		bi = bi*(b+float(i-1));
//		ci = ci*(c+float(i-1));
//		F_out = F_out + std::pow(-1.0,i)*bc*(bi/ci)*std::pow(z,i);
//	}
//	return F_out;
//}

// the hypergeometric function
float Fh(float a, float b, float c, float z){
	int i;
	float Fh_tmp;
	float Fh_out;
	float ai = 1.0;
	float bi = 1.0;
	float ci = 1.0;
	float ni = 1.0; // factorial
	Fh_out = 1.0;
	for (i=1; i<5000; i++){
		ai = ai*(a+float(i-1));
		bi = bi*(b+float(i-1));
		ci = ci*(c+float(i-1));
		ni = ni*float(i);
		Fh_tmp = (ai*bi/ci)*(std::pow(z,i)/ni);
		//std::cout << "F[" << i << "] = " << Fh_tmp << std::endl;
		Fh_out = Fh_out + Fh_tmp;
		if ( std::abs(Fh_tmp) <= 1e-12 ){ break; }
	}
	return Fh_out;
}

// e.g., cylindrial layer
float phi_sf(float r){
	//float Dcc;
	//float Rcc = Dcc/2.0;
	float phi_ext_out;
	phi_ext_out = (M_PI*M_PI)*rho_ss*epsilon_sf*(sigma_sf*sigma_sf) * (
		(63.0/32.0)*std::pow(((Rcc-r)/sigma_sf)*(1.0+(r/Rcc)),-10.0) * Fh(-9.0/2.0,-9.0/2.0,1.0,(r/Rcc)*(r/Rcc))
		-3.0*std::pow(((Rcc-r)/sigma_sf)*(1.0+(r/Rcc)),-4.0) * Fh(-3.0/2.0,-3.0/2.0,1.0,(r/Rcc)*(r/Rcc)) );
	return phi_ext_out;
}

// e.g., wall potential (carbon cylinder)
float phi_ext(float r){
	//float Dcc;
	//float Rcc = Dcc/2.0;
	float phi_ext_out;
	//phi_ext_out = phi_sf(r) + phi_sf(Rcc-r);
	phi_ext_out = phi_sf(r);
	//std::cout << phi_ext_out << std::endl;
	return phi_ext_out;
}

// from Carnahan-Starling (CS) equation of state for NLDFT
//float mu_ex(float rho_b){
//	float y, mu_ex_out;
//	//y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
//	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
//	float den1y = (1.0-y);
//	//mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/std::pow((1.0-y),3.0);
//	//mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/((1.0-y)*(1.0-y)*(1.0-y));
//	mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/(den1y*den1y*den1y);
//	return mu_ex_out;
//}

// The excess hard sphere chemical potential (mu_ex) in the bulk fulid.
// mu_ex is calculated by the PY equation.
float mu_ex(float rho_b){
	float y, mu_ex_out;
	//eta = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	//eta = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
	//float den1e = (1.0-eta);
	float den1y = (1.0-y);
	//mu_ex_out = kb1*T*(-std::log(1-eta) + eta*(14.0 - 13.0*eta + 5.0*eta*eta)/(2.0*std::pow((1.0-eta),3.0)));
	//mu_ex_out = kb1*T*(-std::log(den1e) + eta*(14.0 - 13.0*eta + 5.0*eta*eta)/(2.0*(den1e*den1e*den1e)));
	mu_ex_out = kb1*T*(-std::log(den1y) + y*(14.0 - 13.0*y + 5.0*y*y)/(2.0*den1y*den1y*den1y));
	return mu_ex_out;
}

float mu_b(float rho_b){
	float mu_id, mu_hs, mu_b_out;
	//mu_id = kb1*T*std::log(std::pow(lam,3.0)*rho_b);
	mu_id = kb1*T*std::log((lam*lam*lam)*rho_b);
	mu_hs = mu_id + mu_ex(rho_b);
	mu_b_out = mu_hs - rho_b*alpha;
	return mu_b_out;
}

// For SDA
//float f_ex(float rho_s){
//	float eta, f_ex_out;
//	//eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
//	eta = M_PI*rho_s*(d_hs*d_hs*d_hs)/6.0;
//	float den1e = (1.0-eta);
//	//f_ex_out = kb1*T*eta*(4.0-3.0*eta)/std::pow((1.0-eta),2.0);
//	//f_ex_out = kb1*T*eta*(4.0-3.0*eta)/((1.0-eta)*(1.0-eta));
//	f_ex_out = kb1*T*eta*(4.0-3.0*eta)/(den1e*den1e);
//	return f_ex_out;
//}

// d(f_ex)/d(rho_s), for SDA
//float dfex_per_drhos(float rho_s){
//	float dfex_per_drhos_out;
//	float eta;
//	//eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
//	eta = M_PI*rho_s*(d_hs*d_hs*d_hs)/6.0;
//	float den1e = (1.0-eta);
//	//dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/std::pow((1.0-eta),3.0)*M_PI*std::pow(d_hs,3.0)/6.0;
//	//dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/((1.0-eta)*(1.0-eta)*(1.0-eta))*M_PI*(d_hs*d_hs*d_hs)/6.0;
//	dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/(den1e*den1e*den1e)*(M_PI*(d_hs*d_hs*d_hs)/6.0);
//	return dfex_per_drhos_out;
//}

// d(rho_s)/d(rho), for SDA
//float drhos_per_drho(float *rho, float r1, float r2, float *r, float ra){
//	float w, drhos_per_drho_out;
//	// Percus-Yevick approximation, Tarazona theory
//	w = wi(ra,0) + wi(ra,1)*rho_s(rho,r1,r) + wi(ra,2)*std::pow(rho_s(rho,r1,r),2.0);
//	drhos_per_drho_out = w/(1.0-rho_si(rho,r2,r,1)-2.0*rho_si(rho,r2,r,2)*rho_s(rho,r2,r));
//	return drhos_per_drho_out;
//}

// d(rho_s)/d(rho), modified version, for SDA
//float drhos_per_drho_j(float ra, float rho_sj, float rho_s1j, float rho_s2j){
//	float w, drhos_per_drho_out;
//	// Percus-Yevick approximation, Tarazona theory
//	//w = wi(ra,0) + wi(ra,1)*rho_sj + wi(ra,2)*std::pow(rho_sj,2.0);
//	w = wi(ra,0) + wi(ra,1)*rho_sj + wi(ra,2)*(rho_sj*rho_sj);
//	drhos_per_drho_out = w/(1.0-rho_s1j-2.0*rho_s2j*rho_sj);
//	return drhos_per_drho_out;
//}

float ni(float *rho, float *r, int i, float *n0_j, float *n1_j, float *n2_j, float *n3_j, float *nv1_j, float *nv2_j,
		  float *n0, float *n1, float *n2, float *n3, float *nv1, float *nv2){
	float Rif;
	Rif = d_hs/2.0; // [nm], Rif is the hard-sphere radius of fluid
	//float Ris;
	//Ris = 0.2217/2.0; // [nm] Ris is the hard-sphere radius of solid (for QSDFT)
	// 2.217e-1 [m], The hard sphere diameter of carbon atoms
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
	int j;
	int t;
	//
	float nraj;
	float nxf, nxf2;
	float ntr, ntr2;
	float new_nrad, old_nrad, dnrad;
	float nrho, old_nrho;
	float ny, ny2, old_ny, dny;
	float nypw;
	//
	float praj;
	float pxf, pxf2;
	float ptr, ptr2;
	float new_prad, old_prad, dprad;
	float prho, old_prho;
	float py, py2, old_py, dpy;
	float pypw;
	//
	for (j=0; j<nstep; j++) {
		//std::cout << j << " " << nr[j] << std::endl;
		n0_j[j] = 0.0;
		n1_j[j] = 0.0;
		n2_j[j] = 0.0;
		n3_j[j] = 0.0;
		nv1_j[j] = 0.0;
		nv2_j[j] = 0.0;
		//
		// negative
		//
		nraj = (-r[j]-r[i]);
		nxf2 = (Rif*Rif-nraj*nraj);
		// nxf2 >= 0.0 is in Rif range.
		if ( nxf2 >= 0.0 ){
			old_nrad = 0.0;
			old_nrho = rho[j];
			old_ny = 0.0;
			//
			for (t=j+1; t<nstep; t++){
				ny2 = -r[t]*-r[t] - -r[j]*-r[j];
				ny = std::sqrt(ny2);
				dny = ny - old_ny;
				//
				ntr2 = nxf2 + -r[j]*-r[j];
				ntr = std::sqrt(ntr2);
				if ( ntr < (Dcc-sigma_ss)/2.0 ){
					//
					nxf = std::sqrt(nxf2);
					if ( old_ny <= nxf ){
						if ( ny <= nxf ){
							new_nrad = std::asin(ny/nxf); // radian
						} else {
							new_nrad = M_PI;
						}
						dnrad = new_nrad - old_nrad;
						//
						nrho = 4.0*nxf*(old_nrho*dnrad/2.0 + rho[t]*dnrad/2.0);
						//
						//n0_j[j] += 4.0*(M_PI/180.0*xf)*(old_rho*dtheta/2.0 + rho[t]*dtheta/2.0)/(4.0*M_PI*Rif*Rif);
						n0_j[j] += nrho/(4.0*M_PI*Rif*Rif);
						//
						//n1_j[j] += 4.0*(M_PI/180.0*xf)*(old_rho*dtheta/2.0 + rho[t]*dtheta/2.0)/(4.0*M_PI*Rif);
						n1_j[j] += nrho/(4.0*M_PI*Rif);
						//
						//n2_j[j] += 4.0*(M_PI/180.0*xf)*(old_rho*dtheta/2.0 + rho[t]*dtheta/2.0);
						n2_j[j] += nrho;
						//
						n3_j[j] += 4.0*nxf*(old_nrho*(dny/2.0)*std::cos(old_nrad) + rho[t]*(dny/2.0)*std::cos(new_nrad));
						//
						//nv1_j[j] += 4.0*(M_PI/180.0*xf)*(old_rho*dtheta/2.0 + rho[t]*dtheta/2.0)/(4.0*M_PI*Rif)*(raj/Rif);
						nv1_j[j] += nrho/(4.0*M_PI*Rif)*(nraj/Rif); // total is 0.0 on x-y plane.
						//
						//nv1_j[j] += 4.0*(M_PI/180.0*xf)*(old_rho*dtheta/2.0 + rho[t]*dtheta/2.0)*(raj/Rif);
						nv2_j[j] += nrho*(nraj/Rif); // total is 0.0 on x-y plane.
					}
				} else {
					//
					nypw = std::sqrt((Dcc-sigma_ss)/2.0*(Dcc-sigma_ss)/2.0 - -r[j]*-r[j]);
					if ( old_ny < nypw ){
						if ( ny < nypw ){
							new_nrad = std::asin(ny/nypw); // radian
						} else {
							new_nrad = M_PI;
						}
						n3_j[j] += 4.0*nypw*(old_nrho*(dny/2.0)*std::cos(old_nrad) + rho[t]*(dny/2.0)*std::cos(new_nrad));
					}
				}
				old_ny = ny;
				old_nrad = new_nrad;
				old_nrho = rho[t];
			}
		}
		//
		// positive
		//
		praj = (r[j]-r[i]);
		pxf2 = (Rif*Rif-praj*praj);
		// pxf2 >= 0.0 is in Rif range.
		if ( pxf2 >= 0.0  ){
			old_prad = 0.0;
			old_prho = rho[j];
			old_py = 0.0;
			//
			for (t=j+1; t<nstep; t++){
				py2 = r[t]*r[t] - r[j]*r[j];
				py = std::sqrt(py2);
				dpy = py - old_py;
				//
				ptr2 = pxf2 + r[j]*r[j];
				ptr = std::sqrt(ptr2);
				if ( ptr < (Dcc-sigma_ss)/2.0 ){
					//
					pxf = std::sqrt(pxf2);
					if ( old_py <= pxf ){
						if ( py <= pxf ){
							new_prad = std::asin(py/pxf); // radian
						} else {
							new_prad = M_PI;
						}
						dprad = new_prad - old_prad;
						//
						prho = 4.0*pxf*(old_prho*dprad/2.0 + rho[t]*dprad/2.0);
						//
						//n0_j[j] += 4.0*(M_PI/180.0*xf)*(old_rho*dtheta/2.0 + rho[t]*dtheta/2.0)/(4.0*M_PI*Rif*Rif);
						n0_j[j] += prho/(4.0*M_PI*Rif*Rif);
						//
						//n1_j[j] += 4.0*(M_PI/180.0*xf)*(old_rho*dtheta/2.0 + rho[t]*dtheta/2.0)/(4.0*M_PI*Rif);
						n1_j[j] += prho/(4.0*M_PI*Rif);
						//
						//n2_j[j] += 4.0*(M_PI/180.0*xf)*(old_rho*dtheta/2.0 + rho[t]*dtheta/2.0);
						n2_j[j] += prho;
						//
						n3_j[j] += 4.0*pxf*(old_prho*(dpy/2.0)*std::cos(old_prad) + rho[t]*(dpy/2.0)*std::cos(new_prad));
						//
						//nv1_j[j] += 4.0*(M_PI/180.0*xf)*(old_rho*dtheta/2.0 + rho[t]*dtheta/2.0)/(4.0*M_PI*Rif)*(raj/Rif);
						nv1_j[j] += prho/(4.0*M_PI*Rif)*(praj/Rif); // total is 0.0 on x-y plane.
						//
						//nv1_j[j] += 4.0*(M_PI/180.0*xf)*(old_rho*dtheta/2.0 + rho[t]*dtheta/2.0)*(raj/Rif);
						nv2_j[j] += prho*(praj/Rif); // total is 0.0 on x-y plane.
					}
				} else {
					//
					pypw = std::sqrt((Dcc-sigma_ss)/2.0*(Dcc-sigma_ss)/2.0 - r[j]*r[j]);
					if ( old_py < pypw ){
						if ( py < pypw ){
							new_prad = std::asin(py/pypw); // radian
						} else {
							new_prad = M_PI;
						}
						n3_j[j] += 4.0*pypw*(old_prho*(dpy/2.0)*std::cos(old_prad) + rho[t]*(dpy/2.0)*std::cos(new_prad));
					}
				}
				old_py = py;
				old_prad = new_prad;
				old_prho = rho[t];
			}
		}
		//
		//n0_j[j] = (rho[j])/(4.0*M_PI*Rif*Rif)*(2.0*M_PI*x);
		//n0_j[j] = (rho[j])/(2.0*Rif*Rif)*x;
		//n0_j[j] = (rho[j])/(2.0*Rif*Rif)*xf + (rho_ssq(r[j])+rho_ssq(H-r[j]))/(2.0*Ris*Ris)*xs; // For QSDFT
		//n0_j[j] = (rho[j])/(2.0*Rif*Rif)*xf;
		//
		//n1_j[j] = (rho[j])/(4.0*M_PI*Rif)*(2.0*M_PI*x);
		//n1_j[j] = (rho[j])/(2.0*Rif)*x;
		//n1_j[j] = (rho[j])/(2.0*Rif)*xf + (rho_ssq(r[j])+rho_ssq(H-r[j]))/(2.0*Ris)*xs; // For QSDFT
		//n1_j[j] = (rho[j])/(2.0*Rif)*xf;
		//
		//n2_j[j] = (rho[j])*(2.0*M_PI*xf) + (rho_ssq(r[j])+rho_ssq(H-r[j]))*(2.0*M_PI*xs); // For QSDFT
		//n2_j[j] = (rho[j])*(2.0*M_PI*xf);
		//
		//n3_j[j] = (rho[j])*(M_PI*xf*xf) + (rho_ssq(r[j])+rho_ssq(H-r[j]))*(M_PI*xs*xs); // For QSDFT
		//n3_j[j] = (rho[j])*(M_PI*xf*xf);
		//
		//nv1_j[j] = (rho[j])/(4.0*M_PI*Rif)*(raj/Rif)*(2.0*M_PI*x);
		//nv1_j[j] = (rho[j])/(2.0*Rif)*(raj/Rif)*xf + (rho_ssq(r[j])+rho_ssq(H-r[j]))/(2.0*Ris)*(raj/Ris)*xs; // For QSDFT
		//
		//nv2 = integral rho(z)*(z/Ri)*(2.0*M_PI*x);
		//nv2_j[j] = (rho[j])*(raj/Rif)*(2.0*M_PI*xf) + (rho_ssq(r[j])+rho_ssq(H-r[j]))*(raj/Ris)*(2.0*M_PI*xs); // For QSDFT
		
		//
		//std::cout << i << ", " << j << ", " << r[i] << ", " << r[j] << ", " << raj << ", " << xf << std::endl;
		//std::cout << "i, j, rho[j], n0_j[j], n1_j[j], n2_j[j], n3_j[j], nv1_j[j], nv2_j[j]" << std::endl;
		//std::cout << i << ", " << j << ", " << rho[j] << ", " << n0_j[j] << ", " << n1_j[j] << ", " << n2_j[j] << ", " << n3_j[j] << ", " << nv1_j[j] << ", " << nv2_j[j] << std::endl;
	}
    //integral_trapezoidal(float *f, int n, float dx)
	n0[i] = integral_trapezoidal(n0_j, nstep-1, dr) + 2.0*n0_j[0]*(dr/2.0);
	n1[i] = integral_trapezoidal(n1_j, nstep-1, dr) + 2.0*n1_j[0]*(dr/2.0);
	n2[i] = integral_trapezoidal(n2_j, nstep-1, dr) + 2.0*n2_j[0]*(dr/2.0);
	n3[i] = integral_trapezoidal(n3_j, nstep-1, dr) + 2.0*n3_j[0]*(dr/2.0);
	nv1[i] = integral_trapezoidal(nv1_j, nstep-1, dr) + 2.0*nv1_j[0]*(dr/2.0);
	nv2[i] = integral_trapezoidal(nv2_j, nstep-1, dr) + 2.0*nv2_j[0]*(dr/2.0);
	//
	//integral_simpson(float *f, int n, float dx)
	//n0[i] = integral_simpson(n0_j, nstep-1, dr) + 2.0*n0_j[0]*(dr/2.0);
	//n1[i] = integral_simpson(n1_j, nstep-1, dr) + 2.0*n1_j[0]*(dr/2.0);
	//n2[i] = integral_simpson(n2_j, nstep-1, dr) + 2.0*n2_j[0]*(dr/2.0);
	//n3[i] = integral_simpson(n3_j, nstep-1, dr) + 2.0*n3_j[0]*(dr/2.0);
	//nv1[i] = integral_simpson(nv1_j, nstep-1, dr) + 2.0*nv1_j[0]*(dr/2.0);
	//nv2[i] = integral_simpson(nv2_j, nstep-1, dr) + 2.0*nv2_j[0]*(dr/2.0);
	//
	//std::cout << "i, r[i], j, r[j], nraj, nxf, n0[i], n1[i], n2[i], n3[i], nv1[i], nv2[i]" << std::endl;
	//std::cout << i << ", " << r[i] << ", " << j-1 << ", " << r[j-1] << ", " << nraj << ", " << nxf << ", " << n0[i] << ", " << n1[i] << ", " << n2[i] << ", " << n3[i] << ", " << nv1[i] << ", " << nv2[i] << ", " << std::endl;
	return 0;
}

// c(1)(r) = (-k*T)*dF_HS(rho)/drho(r) = - integral f_ex dr' = dfex
// f_ex = sigma dphi/dn * w(r-r')
// phi = phi1 + phi2 + phi3 (RSLT2, the second modification)
// phi1 = (-n0*ln(1-n3)
// phi2 = (n1*n2-nv1*nv2)/(1-n3)
// phi3 = n2*n2*n2/(24*pi*(1-n3)*(1-n3)) * (1-3*|nv2/n2|*|nv2/n2|+2*|nv2/n2|*|nv2/n2|*|nv2/n2|)
// Only the z-axis components of nv1 and nv2 (with positive and negative) are remained due to symmetry.
float dfex(float *r, int i, float *n0, float *n1, float *n2, float *n3, float *nv1, float *nv2){
	float dfex_out;
	float sxi;
	float sign;
	float Rif;
	Rif = d_hs/2.0; // [nm] Rif is the hard-sphere radius of fluid
	//float Ris;
	//Ris = 0.2217/2.0; // [nm] Ris is the hard-sphere radius of solid (for QSDFT)
	// 2.217e-1 [m], The hard sphere diameter of carbon atoms
	//
	// Memo
	// df(x)/dx = [d/dx1,...,d/dxn]t * [f1(x),...,fn(x)]
	// df(r)/dr = [d/dx,d/dy,d/dz]t * [fx(r),fy(r),fz(r)]
	// dr/dr = [d/dx,d/dy,d/dz]t * [x,y,z] = [dx/dx,dy/dy,dz/dz] = [1,1,1] = [I] = 1
	// r > 0: d|r|/dr = dr/dr = [I] = 1
	// r < 0: d|r|/dr = d(-r)/dr = -[I] = -1
	//
	int j;
	int t;
	//
	float nraj;
	float nx, nx2;
	float ntr, ntr2;
	float new_nrad, old_nrad, dnrad;
	//float nrho, old_nrho;
	float ny, ny2, old_ny, dny;
	float nypw;
	float old_ndpn0, old_ndpn1, old_ndpn2, old_ndpn3;
	float old_ndpnv1, old_ndpnv2;
	float new_ndpn0, new_ndpn1, new_ndpn2, new_ndpn3;
	float new_ndpnv1, new_ndpnv2;
	//
	float praj;
	float px, px2;
	float ptr, ptr2;
	float new_prad, old_prad, dprad;
	//float prho, old_prho;
	float py, py2, old_py, dpy;
	float pypw;
	float old_pdpn0, old_pdpn1, old_pdpn2, old_pdpn3;
	float old_pdpnv1, old_pdpnv2;
	float new_pdpn0, new_pdpn1, new_pdpn2, new_pdpn3;
	float new_pdpnv1, new_pdpnv2;
	//
	float dphi_per_n0, dphi_per_n0_j[nstep];
	float dphi_per_n1, dphi_per_n1_j[nstep];
	float dphi_per_n2, dphi_per_n2_j[nstep];
	float dphi_per_n3, dphi_per_n3_j[nstep];
	float dphi_per_nv1, dphi_per_nv1_j[nstep];
	float dphi_per_nv2, dphi_per_nv2_j[nstep];
	//
	for (j=0; j<nstep; j++) {
		dphi_per_n0_j[j] = 0.0;
		dphi_per_n1_j[j] = 0.0;
		dphi_per_n2_j[j] = 0.0;
		dphi_per_n3_j[j] = 0.0;
		dphi_per_nv1_j[j] = 0.0;
		dphi_per_nv2_j[j] = 0.0;
		//
		// negative
		//
		if ( n2[j] > 0.0 ){
			sxi = std::abs(-nv2[j]/n2[j]);
		} else {
			sxi = 0.0;
		}
		//
		if ( -nv2[j]/n2[j] >= 0.0 ){
			// nv2/n2 > 0  ->  nv2/n2 = sxi
			sign = 1.0;
		} else if ( -nv2[j]/n2[j] < 0.0 ) {
			// nv2/n2 < 0  ->  -nv2/n2 = -sxi
			sign = -1.0;
		}
		//
		nraj = (-r[j]-r[i]);
		nx2 = (Rif*Rif-nraj*nraj);
		// nx2 >= 0.0 is in Rif range.
		if ( nx2 >= 0.0 ){
			old_nrad = 0.0;
			old_ndpn0 = -std::log(1.0-n3[j]);
			//std::cout << old_dpn0 << std::endl;
			old_ndpn1 = n2[j]/(1.0-n3[j]);
			//std::cout << old_dpn1 << std::endl;
			old_ndpn2 = ( n1[j]/(1.0-n3[j])
						+ 3.0*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
						+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
						* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(nv2[j]/(n2[j]*n2[j])*sign)
					   );
			//std::cout << old_dpn2 << std::endl;
			old_ndpn3 = ( n0[j]/(1.0-n3[j])
						+ (n1[j]*n2[j] - nv1[j]*nv2[j])/((1.0-n3[j])*(1.0-n3[j])) 
						+ 2.0*n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
					   );
			//std::cout << old_dpn3 << std::endl;
			old_ndpnv1 = nv2[j]/(1.0-n3[j]);
			//std::cout << old_dpnv1 << std::endl;
			old_ndpnv2 = ( nv1[j]/(1.0-n3[j])
						+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
						* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(1.0/n2[j])
					   );
			//std::cout << old_dpnv2 << std::endl;
			//
			for (t=j+1; t<nstep; t++){
				ny2 = -r[t]*-r[t] - -r[j]*-r[j];
				ny = std::sqrt(ny2);
				dny = ny - old_ny;
				//
				ntr2 = nx2 + -r[j]*-r[j];
				ntr = std::sqrt(ntr2);
				if ( ntr < (Dcc-sigma_ss)/2.0 ) {
					//
					nx = std::sqrt(nx2);
					if ( old_ny <= nx ){
						if ( ny <= nx ){
							new_nrad = std::asin(ny/nx); // radian
						} else {
							new_nrad = M_PI;
						}
						dnrad = new_nrad - old_nrad;
						//
						// dphi/dn0
						new_ndpn0 = -std::log(1.0-n3[t]);
						dphi_per_n0_j[j] += 4.0*nx*(old_ndpn0*dnrad/2.0 + new_ndpn0*dnrad/2.0)/(4.0*M_PI*Rif*Rif); // PHYSICAL REVIEW E, VOLUME 64, 011602
						//
						// dphi/dn1
						new_ndpn1 = n2[t]/(1.0-n3[t]);
						dphi_per_n1_j[j] += 4.0*nx*(old_ndpn1*dnrad/2.0 + new_ndpn1*dnrad/2.0)/(4.0*M_PI*Rif); // PHYSICAL REVIEW E, VOLUME 64, 011602
						//
						// dphi/dn2, q=2 case, RSLT version // PHYSICAL REVIEW E 64 011602
						new_ndpn2 = ( n1[t]/(1.0-n3[t])
							+ 3.0*n2[t]*n2[t]/(24.0*M_PI*(1.0-n3[t])*(1.0-n3[t]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
							+ n2[t]*n2[t]*n2[t]/(24.0*M_PI*(1.0-n3[t])*(1.0-n3[t])*(1.0-n3[t]))
							* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(nv2[t]/(n2[t]*n2[t])*sign)
						);
						dphi_per_n2_j[j] += 4.0*nx*(old_ndpn2*dnrad/2.0 + new_ndpn2*dnrad/2.0);
						//
						// dphi/dn3, q=2 case, RSLT version // PHYSICAL REVIEW E, VOLUME 64, 011602
						new_ndpn3 = ( n0[t]/(1.0-n3[t])
							+ (n1[t]*n2[t] - nv1[t]*nv2[t])/((1.0-n3[t])*(1.0-n3[t])) 
							+ 2.0*n2[t]*n2[t]*n2[t]/(24.0*M_PI*(1.0-n3[t])*(1.0-n3[t])*(1.0-n3[t]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
						);
						dphi_per_n3_j[j] += 4.0*nx*(old_ndpn3*(dny/2.0)*std::cos(old_nrad) + new_ndpn3*(dny/2.0)*std::cos(new_nrad));
						//
						// dphi/dnv1, total is 0.0 on x-y plane.
						new_ndpnv1 = nv2[t]/(1.0-n3[t]);
						dphi_per_nv1_j[j] += 4.0*nx*(old_ndpnv1*dnrad/2.0 + new_ndpnv1*dnrad/2.0)/(4.0*M_PI*Rif)*(nraj/Rif);
						//
						// dphi/dnv2, total is 0.0 on x-y plane.
						new_ndpnv2 = ( nv1[t]/(1.0-n3[t])
							+ n2[t]*n2[t]*n2[t]/(24.0*M_PI*(1.0-n3[t])*(1.0-n3[t])*(1.0-n3[t]))
							* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(1.0/n2[t])
						);
						dphi_per_nv2_j[j] += 4.0*nx*(old_ndpnv2*dnrad/2.0 + new_ndpnv2*dnrad/2.0)*(nraj/Rif);
					}
				} else {
					//
					nypw = std::sqrt((Dcc-sigma_ss)/2.0*(Dcc-sigma_ss)/2.0 - -r[j]*-r[j]);
					if ( old_ny < nypw ){
						if ( ny < nypw ){
							new_nrad = std::asin(ny/nypw); // radian
						} else {
							new_nrad = M_PI;
						}
						// dphi/dn3, q=2 case, RSLT version // PHYSICAL REVIEW E, VOLUME 64, 011602
						new_ndpn3 = ( n0[t]/(1.0-n3[t])
							+ (n1[t]*n2[t] - nv1[t]*nv2[t])/((1.0-n3[t])*(1.0-n3[t])) 
							+ 2.0*n2[t]*n2[t]*n2[t]/(24.0*M_PI*(1.0-n3[t])*(1.0-n3[t])*(1.0-n3[t]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
						);
						dphi_per_n3_j[j] += 4.0*nypw*(old_ndpn3*(dny/2.0)*std::cos(old_nrad) + new_ndpn3*(dny/2.0)*std::cos(new_nrad));
					}
				}
				old_ny = ny;
				old_nrad = new_nrad;
				old_ndpn0 = new_ndpn0;
				old_ndpn1 = new_ndpn1;
				old_ndpn2 = new_ndpn2;
				old_ndpn3 = new_ndpn3;
				old_ndpnv1 = new_ndpnv1;
				old_ndpnv2 = new_ndpnv2;
			}
		}
		//
		// positive
		//
		if ( n2[j] > 0.0 ){
			sxi = std::abs(nv2[j]/n2[j]);
		} else {
			sxi = 0.0;
		}
		//
		if ( nv2[j]/n2[j] >= 0.0 ){
			// nv2/n2 > 0  ->  nv2/n2 = sxi
			sign = 1.0;
		} else if ( nv2[j]/n2[j] < 0.0 ) {
			// nv2/n2 < 0  ->  -nv2/n2 = -sxi
			sign = -1.0;
		}
		//
		praj = (r[j]-r[i]);
		px2 = (Rif*Rif-praj*praj);
		// px2 >= 0.0 is in Rif range.
		if ( px2 >= 0.0 ){
			old_prad = 0.0;
			old_pdpn0 = -std::log(1.0-n3[j]);
			//std::cout << old_dpn0 << std::endl;
			old_pdpn1 = n2[j]/(1.0-n3[j]);
			//std::cout << old_dpn1 << std::endl;
			old_pdpn2 = ( n1[j]/(1.0-n3[j])
						+ 3.0*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
						+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
						* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(-nv2[j]/(n2[j]*n2[j])*sign)
					   );
			//std::cout << old_dpn2 << std::endl;
			old_pdpn3 = ( n0[j]/(1.0-n3[j])
						+ (n1[j]*n2[j] - nv1[j]*nv2[j])/((1.0-n3[j])*(1.0-n3[j])) 
						+ 2.0*n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
					   );
			//std::cout << old_dpn3 << std::endl;
			old_pdpnv1 = -nv2[j]/(1.0-n3[j]);
			//std::cout << old_dpnv1 << std::endl;
			old_pdpnv2 = ( -nv1[j]/(1.0-n3[j])
						+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
						* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(1.0/n2[j])
					   );
			//std::cout << old_dpnv2 << std::endl;
			//
			for (t=j+1; t<nstep; t++){
				py2 = r[t]*r[t] - r[j]*r[j];
				py = std::sqrt(py2);
				dpy = py - old_py;
				//
				ptr2 = px2 + r[j]*r[j];
				ptr = std::sqrt(ptr2);
				if ( ptr < (Dcc-sigma_ss)/2.0) {
					//
					px = std::sqrt(px2);
					if ( old_py <= px ){
						if ( py <= px ){
							new_prad = std::asin(py/px); // radian
						} else {
							new_prad = M_PI;
						}
						dprad = new_prad - old_prad;
						//
						// dphi/dn0
						new_pdpn0 = -std::log(1.0-n3[t]);
						dphi_per_n0_j[j] += 4.0*px*(old_pdpn0*dprad/2.0 + new_pdpn0*dprad/2.0)/(4.0*M_PI*Rif*Rif); // PHYSICAL REVIEW E, VOLUME 64, 011602
						//
						// dphi/dn1
						new_pdpn1 = n2[t]/(1.0-n3[t]);
						dphi_per_n1_j[j] += 4.0*px*(old_pdpn1*dprad/2.0 + new_pdpn1*dprad/2.0)/(4.0*M_PI*Rif); // PHYSICAL REVIEW E, VOLUME 64, 011602
						//
						// dphi/dn2, q=2 case, RSLT version // PHYSICAL REVIEW E 64 011602
						new_pdpn2 = ( n1[t]/(1.0-n3[t])
							+ 3.0*n2[t]*n2[t]/(24.0*M_PI*(1.0-n3[t])*(1.0-n3[t]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
							+ n2[t]*n2[t]*n2[t]/(24.0*M_PI*(1.0-n3[t])*(1.0-n3[t])*(1.0-n3[t]))
							* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(-nv2[t]/(n2[t]*n2[t])*sign)
						);
						dphi_per_n2_j[j] += 4.0*px*(old_pdpn2*dprad/2.0 + new_pdpn2*dprad/2.0);
						//
						// dphi/dn3, q=2 case, RSLT version // PHYSICAL REVIEW E, VOLUME 64, 011602
						new_pdpn3 = ( n0[t]/(1.0-n3[t])
							+ (n1[t]*n2[t] - nv1[t]*nv2[t])/((1.0-n3[t])*(1.0-n3[t])) 
							+ 2.0*n2[t]*n2[t]*n2[t]/(24.0*M_PI*(1.0-n3[t])*(1.0-n3[t])*(1.0-n3[t]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
						);
						dphi_per_n3_j[j] += 4.0*px*(old_pdpn3*(dpy/2.0)*std::cos(old_prad) + new_pdpn3*(dpy/2.0)*std::cos(new_prad));
						//
						// dphi/dnv1, total is 0.0 on x-y plane.
						new_pdpnv1 = -nv2[t]/(1.0-n3[t]);
						dphi_per_nv1_j[j] += 4.0*px*(old_pdpnv1*dprad/2.0 + new_pdpnv1*dprad/2.0)/(4.0*M_PI*Rif)*(praj/Rif);
						//
						// dphi/dnv2, total is 0.0 on x-y plane.
						new_pdpnv2 = ( -nv1[t]/(1.0-n3[t])
							+ n2[t]*n2[t]*n2[t]/(24.0*M_PI*(1.0-n3[t])*(1.0-n3[t])*(1.0-n3[t]))
							* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(1.0/n2[t])
						);
						dphi_per_nv2_j[j] += 4.0*px*(old_pdpnv2*dprad/2.0 + new_pdpnv2*dprad/2.0)*(praj/Rif);
					}
				} else {
					//
					pypw = std::sqrt((Dcc-sigma_ss)/2.0*(Dcc-sigma_ss)/2.0 - r[j]*r[j]);
					if ( old_py < pypw ){
						if ( py < pypw ){
							new_prad = std::asin(py/pypw); // radian
						} else {
							new_prad = M_PI;
						}
						// dphi/dn3, q=2 case, RSLT version // PHYSICAL REVIEW E, VOLUME 64, 011602
						new_pdpn3 = ( n0[t]/(1.0-n3[t])
							+ (n1[t]*n2[t] - nv1[t]*nv2[t])/((1.0-n3[t])*(1.0-n3[t])) 
							+ 2.0*n2[t]*n2[t]*n2[t]/(24.0*M_PI*(1.0-n3[t])*(1.0-n3[t])*(1.0-n3[t]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)
						);
						dphi_per_n3_j[j] += 4.0*pypw*(old_pdpn3*(dpy/2.0)*std::cos(old_prad) + new_pdpn3*(dpy/2.0)*std::cos(new_prad));
					}
				}
				old_py = py;
				old_prad = new_prad;
				old_pdpn0 = new_pdpn0;
				old_pdpn1 = new_pdpn1;
				old_pdpn2 = new_pdpn2;
				old_pdpn3 = new_pdpn3;
				old_pdpnv1 = new_pdpnv1;
				old_pdpnv2 = new_pdpnv2;
			}
		}
		//
		//sxi = std::abs(nv2[j]/n2[j]);
		//std::cout << j << ", sxi = " << sxi << std::endl;
		//
		// dphi/dn0
		//dphi_per_n0[j] = -std::log(1.0-n3[j])/(4.0*M_PI*Rif*Rif)*(2.0*M_PI*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
		//dphi_per_n0_j[j] = -std::log(1.0-n3[j])/(2.0*Rif*Rif)*x; // PHYSICAL REVIEW E, VOLUME 64, 011602
		//
		// dphi/dn1
		//dphi_per_n1_j[j] = ( n2[j]/(1.0-n3[j]) )/(4.0*M_PI*Rif)*(2.0*M_PI*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
		//dphi_per_n1_j[j] = ( n2[j]/(1.0-n3[j]) )/(2.0*Rif)*x; // PHYSICAL REVIEW E, VOLUME 64, 011602
		//
		// dphi/dn2
		//if ( nv2[j]/n2[j] >= 0.0 ){
		//	// nv2/n2 > 0  ->  nv2/n2 = sxi
		//	sign = 1.0;
		//}else if ( nv2[j]/n2[j] < 0.0 ) {
		//	// nv2/n2 < 0  ->  -nv2/n2 = -sxi
		//	sign = -1.0;
		//}
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
		//dphi_per_n2_j[j] = ( n1[j]/(1.0-n3[j])
		//	+ 3.0*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j]))*(1.0-3.0*sxi*sxi+2.0*sxi*sxi*sxi*sign)
		//	+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//		* (1.0-6.0*sxi*sign+6.0*sxi*sxi)*(-nv2[j]/(n2[j]*n2[j])*sign)
		//)*(2.0*M_PI*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
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
		// dphi/dn3, RSLT2 version
		//dphi_per_n3_j[j] = ( n0[j]/(1.0-n3[j])
		//	+ (n1[j]*n2[j] - nv1[j]*nv2[j])/((1.0-n3[j])*(1.0-n3[j])) 
		//	+ 2.0*n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))*(1.0-3.0*sxi*sxi+2.0*sxi*sxi*sxi)
		//)*(M_PI*x*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
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
		// dphi/dnv1, total is 0.0 on x-y plane.
		//
		// dphi/dnv1
		//dphi_per_nv1_j[j] = ( -nv2[j]/(1.0-n3[j]) )/(4.0*M_PI*Ri)*(raj/Ri)*(2.0*M_PI*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
		//dphi_per_nv1_j[j] = ( -nv2[j]/(1.0-n3[j]) )/(2.0*Rif)*(raj/Rif)*x; // PHYSICAL REVIEW E, VOLUME 64, 011602
		//dphi_per_nv1_j[j] = ( -nnv2[t]/(1.0-nn3[t]) )/(2.0*Rif)*(raj/Rif)*x; // PHYSICAL REVIEW E, VOLUME 64, 011602
		//
		// dphi/dnv2, total is 0.0 on x-y plane.
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
		// dphi/dnv2, q=2 case, RSLT version // PHYSICAL REVIEW E, VOLUME 64, 011602
		//dphi_per_nv2_j[j] = ( -nnv1[t]/(1.0-nn3[t])
		//	+ nn2[t]*nn2[t]*nn2[t]/(24.0*M_PI*(1.0-nn3[t])*(1.0-nn3[t])*(1.0-nn3[t]))
		//	* 2.0*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(1.0/nn2[t])
		//)*(raj/Rif)*(2.0*M_PI*x);
		//
		// dphi/dnv2, q=3 case, RSLT version // PHYSICAL REVIEW E, VOLUME 64, 011602
		//dphi_per_nv2_j[j] = ( -nv1[j]/(1.0-n3[j])
		//	+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//	* 3.0*(1.0-sxi*sxi)*(1.0-sxi*sxi)*(-2.0*sxi*sign)*(1.0/n2[j])
		//)*(raj/Rif)*(2.0*M_PI*x);
		//
		// dphi/dnv2, RSLT2 version
		//dphi_per_nv2_j[j] = ( -nv1[j]/(1.0-n3[j])
		//	+ n2[j]*n2[j]*n2[j]/(24.0*M_PI*(1.0-n3[j])*(1.0-n3[j])*(1.0-n3[j]))
		//		* (1.0-6.0*sxi*sign+6.0*sxi*sxi)*(sign*1.0/n2[j])
		//)*(raj/Rif)*(2.0*M_PI*x); // PHYSICAL REVIEW E, VOLUME 64, 011602
		//
		// dphi/dnv2, The modified fundamental-measure theory (MFMT) // Langmuir 2008, 24, 12431-12439
		//dphi_per_nv2_j[j] = ( -nv1[j]/(1.0-n3[j])
		//	-2.0*n2[j]*nv2[j]*std::log(1.0-n3[j])/(12.0*M_PI*n3[j]*n3[j])
		//	-2.0*n2[j]*nv2[j]/(12.0*M_PI*n3[j]*(1.0-n3[j])*(1.0-n3[j]))
		//)*(raj/Rif)*(2.0*M_PI*x);
	}
	//
    //integral_trapezoidal(float *f, int n, float dx)
	dphi_per_n0  = integral_trapezoidal(dphi_per_n0_j, nstep-1, dr) + 2.0*dphi_per_n0_j[0]*(dr/2.0);
	dphi_per_n1  = integral_trapezoidal(dphi_per_n1_j, nstep-1, dr) + 2.0*dphi_per_n1_j[0]*(dr/2.0);
	dphi_per_n2  = integral_trapezoidal(dphi_per_n2_j, nstep-1, dr) + 2.0*dphi_per_n2_j[0]*(dr/2.0);
	dphi_per_n3  = integral_trapezoidal(dphi_per_n3_j, nstep-1, dr) + 2.0*dphi_per_n3_j[0]*(dr/2.0);
	dphi_per_nv1 = integral_trapezoidal(dphi_per_nv1_j, nstep-1, dr) + 2.0*dphi_per_nv1_j[0]*(dr/2.0);
	dphi_per_nv2 = integral_trapezoidal(dphi_per_nv2_j, nstep-1, dr) + 2.0*dphi_per_nv2_j[0]*(dr/2.0);
	//
	//integral_simpson(float *f, int n, float dx)
	//dphi_per_n0 = integral_simpson(dphi_per_n0_j, nstep-1, dr) + 2.0*dphi_per_n0_j[0]*(dr/2.0);
	//dphi_per_n1 = integral_simpson(dphi_per_n1_j, nstep-1, dr) + 2.0*dphi_per_n1_j[0]*(dr/2.0);
	//dphi_per_n2 = integral_simpson(dphi_per_n2_j, nstep-1, dr) + 2.0*dphi_per_n2_j[0]*(dr/2.0);
	//dphi_per_n3 = integral_simpson(dphi_per_n3_j, nstep-1, dr) + 2.0*dphi_per_n3_j[0]*(dr/2.0);
	//dphi_per_nv1 = integral_simpson(dphi_per_nv1_j, nstep-1, dr) + 2.0*dphi_per_nv1_j[0]*(dr/2.0);
	//dphi_per_nv2 = integral_simpson(dphi_per_nv2_j, nstep-1, dr) + 2.0*dphi_per_nv2_j[0]*(dr/2.0);
	//
	//std::cout << "i, dphi_per_n0, dphi_per_n1, dphi_per_n2, dphi_per_n3, dphi_per_nv1, dphi_per_nv2" << std::endl;
	//std::cout << i << ", " << dphi_per_n0 << "," << dphi_per_n1 << "," << dphi_per_n2 << "," << dphi_per_n3 << "," << dphi_per_nv1 << "," << dphi_per_nv2 << "," << std::endl;
	//
	dfex_out = -(dphi_per_n0 + dphi_per_n1 + dphi_per_n2 + dphi_per_n3 + dphi_per_nv1 + dphi_per_nv2);
	//std::cout << dfex_out << std::endl;
	return dfex_out;
}

float calc_alpha(float *r){
	int i,j,k,t;
	float ra;
	float raj;
	float rak;
	float alpha_other_method;
	float alpha_int_i[nstep];
	float alpha_int_j[nstep];
	float alpha_int_k[nhmesh];
	float alpha_int_t[nrmesh];
	float x,y;
	float drad = M_PI/nrmesh;
	for (i=0; i<nstep; i++) {
		alpha_int_i[i] = 0.0;
		for (j=0; j<nstep; j++) {
			alpha_int_j[j] = 0.0;
			for (k=0; k<nhmesh; k++) {
				rak = dh*float(k);
				alpha_int_k[k] = 0.0;
				for (t=0; t<nrmesh; t++) {
					alpha_int_t[t] = 0.0;
					x = r[j]*std::cos(drad*float(t));
					y = r[j]*std::sin(drad*float(t));
					raj = x - r[i];
					ra = raj*raj + y*y + rak*rak;
					ra = std::sqrt(ra);
					alpha_int_t[t] = -phi_att(ra);
					//alpha_int_k[k] += -phi_att(ra);
				}
				//integral_simpson(float *f, int n, float dx)
				alpha_int_k[k] = integral_simpson(alpha_int_t, nrmesh-1, drad);
				//integral_trapezoidal(float *f, int n, float dx)
				//alpha_int_k[k] = integral_trapezoidal(alpha_int_t, nrmesh-1, drad);
			}
			//integral_simpson(float *f, int n, float dx)
			alpha_int_j[j]  = 2.0*r[j]*integral_simpson(alpha_int_k, nhmesh-1, dh)*2.0;
			//alpha_int_j[j]  = 2.0*drad*r[j]*integral_simpson(alpha_int_k, nhmesh-1, dh)*2.0;
			//integral_trapezoidal(float *f, int n, float dx)
			//alpha_int_j[j]  = 2.0*r[j]*integral_trapezoidal(alpha_int_k, nhmesh-1, dh)*2.0;
			//alpha_int_j[j]  = 2.0*drad*r[j]*integral_trapezoidal(alpha_int_k, nhmesh-1, dh)*2.0;
		}
		//integral_simpson(float *f, int n, float dx)
		alpha_int_i[i] = integral_simpson(alpha_int_j, nstep-1, dr) + alpha_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
		//integral_trapezoidal(float *f, int n, float dx)
		//alpha_int_i[i] = integral_trapezoidal(alpha_int_j, nstep-1, dr) + alpha_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	}
	//alpha_other_method = alpha_int_i[0]; // check
	//integral_simpson(float *f, int n, float dx)
	alpha_other_method = integral_simpson(alpha_int_i, nstep-1, dr) / ((Dcc-sigma_ss)/2.0);
	//integral_trapezoidal(float *f, int n, float dx)
	//alpha_other_method = integral_trapezoidal(alpha_int_i, nstep-1, dr) / ((Dcc-sigma_ss)/2.0);
	//
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "average alpha of other method = " << alpha_other_method << " in (carbon) cylinder" << std::endl;
	return alpha_other_method;
}

float phi_att_int(float *r, float *phi_att_int_ij){
	int i,j,k,t;
	float ra;
	float raj;
	float rak;
	float rho_phi_int_k[nhmesh];
	float rho_phi_int_t[nrmesh];
	float x,y;
	float drad = M_PI/nrmesh;
	for (i=0; i<nstep; i++) {
		for (j=0; j<nstep; j++) {
			for (k=0; k<nhmesh; k++) {
				rak = dh*float(k);
				//rho_phi_int_k[k] = 0.0;
				for (t=0; t<nrmesh; t++) {
					x = r[j]*std::cos(drad*float(t));
					y = r[j]*std::sin(drad*float(t));
					raj = x - r[i];
					ra = raj*raj + y*y + rak*rak;
					ra = std::sqrt(ra);
					rho_phi_int_t[t] = phi_att(ra);
					//rho_phi_int_k[k] += phi_att(ra);
				}
				//integral_simpson(float *f, int n, float dx)
				rho_phi_int_k[k] = integral_simpson(rho_phi_int_t, nrmesh-1, drad);
				//integral_trapezoidal(float *f, int n, float dx)
				//rho_phi_int_k[k] = integral_trapezoidal(rho_phi_int_t, nrmesh-1, drad);
			}
			//integral_simpson(float *f, int n, float dx)
			phi_att_int_ij[i*nstep+j] = integral_simpson(rho_phi_int_k, nhmesh-1, dh)*2.0;
			//integral_trapezoidal(float *f, int n, float dx)
			//phi_att_int_ij[i*nstep+j] = integral_trapezoidal(rho_phi_int_k, nhmesh-1, dh)*2.0;
		}
	}
	return 0;
}

// xi include kb1*T*(std::log(rho_b)) type.
// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
float xi(float *rho, float *r, int i, float rho_b, float *phi_att_int_ij, float *rho_phi_int, float *phi_ext_i){
	int j;
	//float ra;
	//float raj;
	//float rak;
	float rho_phi_int_j[nstep];
	//float x,y;
	//float drad = M_PI/nrmesh;
	for (j=0; j<nstep; j++) {
		rho_phi_int_j[j]  = 2.0*r[j]*rho[j]*phi_att_int_ij[i*nstep+j];
	}
	//integral_simpson(float *f, int n, float dx)
	//rho_dfex_int[i] = integral_simpson(rho_dfex_int_j, nstep, dr);
	rho_phi_int[i]  = integral_simpson(rho_phi_int_j, nstep-1, dr) + rho_phi_int_j[0]/(2.0*M_PI*r[0])*M_PI*(dr/2.0)*(dr/2.0);
	//
	float xi_out;
	//xi_out = kb1*T*std::log(rho_b) + mu_ex(rho_b) - rho_b*alpha - phi_ext(r[i]) - f_ex(rho_sj[i]) - rho_dfex_int - rho_phi_int; // old ver.1.1.1
	//xi_out = ( - rho_b*alpha - rho_dfex_int[i] - f_ex(rho_sj[i]) ) + ( mu_ex(rho_b) - rho_phi_int[i] ) + ( kb1*T*std::log(rho_b) - phi_ext(r[i]) );
	xi_out = ( - rho_b*alpha ) + ( mu_ex(rho_b) - rho_phi_int[i] ) + ( kb1*T*std::log(rho_b) - phi_ext_i[i] );
	// debug
	//std::cout << "i, xi_out, -rho_b*alpha, mu_ex(rho_b), -rho_phi_int[i], kb1*T*std::log(rho_b), -phi_ext_i[i]" << std::endl;
	//std::cout << i << ", " << xi_out << ", " << -rho_b*alpha << ", " << mu_ex(rho_b) << ", " << -rho_phi_int[i] << ", " << kb1*T*std::log(rho_b) << ", " << - phi_ext_i[i] << std::endl;
	//std::cout << "xi, (kb1*T)*log(rho_b), mu_ex(rho_b), -rho_b*alpha, -phi_ext(r[i]), -f_ex(rho_s(rho,r[i],r)), -rho_dfex_int, -rho_phi_int" << std::endl;
	//std::cout << xi_out << ", " << kb1*T*std::log(rho_b) << ", " << mu_ex(rho_b) << ", " << -rho_b*alpha << ", " << -phi_ext(r[i]) << ", " << -f_ex(rho_sj[i]) << ", " << -rho_dfex_int << ", " << -rho_phi_int << std::endl;
	//if ( std::isnan(rho_sj[i]) || std::isnan(f_ex(rho_sj[i])) || std::isnan(rho_dfex_int) || std::isnan(rho_phi_int) ){
	//	std::cout << i << ", " << rho[i] << ", " << rho_sj[i] << ", " << f_ex(rho_sj[i]) << ", " << rho_dfex_int << ", " << rho_phi_int << std::endl;
	//	std::exit(1);
	//}
	return xi_out;
}

// For SDA, from Carnahan-Starling (CS) equation
//float press_hs(float rho_b){
//	float y, press_hs_out;
//	//y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
//	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
//	float den1y = (1.0-y);
//	//press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/std::pow((1.0-y),3.0);
//	//press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/((1.0-y)*(1.0-y)*(1.0-y));
//	press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/(den1y*den1y*den1y);
//	return press_hs_out;
//}

// For FMT, from Percus-Yevick (PY) equation
// http://www.sklogwiki.org/SklogWiki/index.php/Exact_solution_of_the_Percus_Yevick_integral_equation_for_hard_spheres
float press_hs(float rho_b){
	float eta, press_hs_out;
	//eta = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	eta = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
	float den1e = (1.0-eta);
	//press_hs_out = rho_b*kb1*T* (1.0 + eta + eta*eta)/std::pow((1.0-eta),3.0);
	//press_hs_out = rho_b*kb1*T* (1.0 + eta + eta*eta)/((1.0-eta)*(1.0-eta)*(1.0-eta));
	press_hs_out = rho_b*kb1*T* (1.0 + eta + eta*eta)/(den1e*den1e*den1e);
	return press_hs_out;
}

float Maxwell_construction(void){
	int i,j;
	int iter_max_drhob0 = 500000;
	int iter_max_dmue = 50000;
	float drhob0 = 0.00005;
	float dmue = 0.0005;
	float threshold_diff = 0.05;
	float threshold_find = 0.05;
	//
	float mu_b_per_epsilon_ff[iter_max_drhob0];
	float mu_e_per_epsilon_ff;
	float diff,diffp;
	int flag;
	float rho_b0_out;
	float rho_b0_gas, rho_b0_metastable, rho_b0_liquid;
	float press_b0;
	//
	// rho_b vs. mu_b/epsilon_ff
	std::ofstream ofs("./Maxwell_construction_data.txt");
	ofs << "# Chemical_potential(mu_b/epsilon_ff), Density(rho_b*d_hs^3)" << std::endl;
	for (i=0; i<iter_max_drhob0; i++){
		rho_b0_out = drhob0*float(i+1.0);
		mu_b_per_epsilon_ff[i] = mu_b(rho_b0_out)/epsilon_ff;
		//ofs << mu_b_per_epsilon_ff[i] << ", " << rho_b0_out*std::pow(d_hs,3.0) << std::endl;
		ofs << mu_b_per_epsilon_ff[i] << ", " << rho_b0_out*(d_hs*d_hs*d_hs) << std::endl;
		//std::cout << "rho_b0 = "<< rho_b0_out << ", mu_b/epsilon_ff = " << mu_b_per_epsilon_ff[i] << std::endl;
	}
	// Maxwell equal area rule
	for (j=0; j<iter_max_dmue; j++){
		mu_e_per_epsilon_ff = dmue*float(j+1.0) - 12.0;
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
		rho_b0_out = drhob0*float(j+1.0);
		if (std::abs(diff) <= threshold_diff) {
			//std::cout << "mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << ", diff = " << diff << std::endl;
			break;
		}
	}
	// find rho_b0
	flag = 0;
	for (i=0; i<iter_max_drhob0; i++){
		rho_b0_out = drhob0*float(i+1.0);
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
	std::cout << "Bulk pressure, P0 = " << press_b0*kb*1e27 << " [Pa] = " << press_b0*kb*1e27/101325.0 << " [atm], (rho_b0 = " << rho_b0_out << ")" <<std::endl;
	std::cout << std::endl;
	std::cout << "gas phase   : rho_b0_gas        = " << rho_b0_gas        << ", rho_b0_gas*d_hs^3        = " << rho_b0_gas*std::pow(d_hs,3.0)        << std::endl;
	std::cout << "metastable  : rho_b0_metastable = " << rho_b0_metastable << ", rho_b0_metastable*d_hs^3 = " << rho_b0_metastable*std::pow(d_hs,3.0) << std::endl;
	std::cout << "liquid phase: rho_b0_liquid     = " << rho_b0_liquid     << ", rho_b0_liquid*d_hs^3     = " << rho_b0_liquid*std::pow(d_hs,3.0)     << std::endl;
	return rho_b0_out;
}

// grand potential for SDA
//float omega(float *rho, float *r, float *rho_dfex_int, float *rho_phi_int){
//	float omega_out;
//	float omega1, omega2, omega3;
//	int i;
//	int omega_nstep = (nstep-2)/2;
//	float rho_x_rho_dfex_int[omega_nstep];
//	float rho_x_rho_phi_int[omega_nstep];
//	for (i=0; i<=omega_nstep; i++){
//		rho_x_rho_dfex_int[i] = rho[i] * rho_dfex_int[i];
//		rho_x_rho_phi_int[i]  = rho[i] * rho_phi_int[i];
//	}
//	omega1 = -(kb1*T) * integral_simpson(rho, omega_nstep, dr);
//	omega2 = -integral_simpson(rho_x_rho_dfex_int, omega_nstep, dr);
//	omega3 = -0.5 * integral_simpson(rho_x_rho_phi_int, omega_nstep, dr);
//	omega_out = (omega1 + omega2 + omega3) * 2.0 / epsilon_ff;
//	return omega_out;
//}

float fex(int i, float *n0, float *n1, float *n2, float *n3, float *nv1, float *nv2){
	float fex_out;
	float sxi = std::abs(nv2[i]/n2[i]);
	float phi1, phi2, phi3;
	phi1 = -n0[i]*std::log(1.0-n3[i]);
	phi2 = (n1[i]*n2[i] - nv1[i]*nv2[i])/(1.0-n3[i]);
	//
	//phi3 = ((1.0/3.0)*n2[i]*n2[i]*n2[i] - n2[i]*(nv2[i]*nv2[i]))/(8.0*M_PI*(1.0-n3[i])*(1.0-n3[i]));
	//
	// RSLT1, q=2
	phi3 = n2[i]*n2[i]*n2[i]/(24.0*M_PI*(1.0-n3[i])*(1.0-n3[i]))*(1.0-sxi*sxi)*(1.0-sxi*sxi);
	//
	// RSLT1, q=3
	//phi3 = n2[i]*n2[i]*n2[i]/(24.0*M_PI*(1.0-n3[i])*(1.0-n3[i]))*(1.0-sxi*sxi)*(1.0-sxi*sxi)*(1.0-sxi*sxi);
	//
	// RSLT2 version // PHYSICAL REVIEW E 64 011602
	//if ( nv2[i]/n2[i] < 0.0 ){
	//	sxi = sxi*-1.0;
	//}
	//phi3 = n2[i]*n2[i]*n2[i]/(24.0*M_PI*(1.0-n3[i])*(1.0-n3[i]))*(1.0-3.0*sxi*sxi+2.0*sxi*sxi*sxi);
	//
	fex_out = phi1 + phi2 + phi3;
	return fex_out;
}

// grand potential for FMT
float omega(float *rho, float *r, float *fex_i, float *rho_phi_int, float rho_b){
	float omega_out;
	omega_out = 1.0;
	float omega1, omega2, omega3, omega4;
	int i;
	float fid[nstep];
	float rho_x_rho_phi_int[nstep];
	float rho_x_phi_ext_mu[nstep];
	float mu = (kb1*T)*std::log(rho_b*lam*lam*lam) + mu_ex(rho_b) - rho_b*alpha;
	float tpi = 2.0*M_PI;
	for (i=0; i<nstep; i++){
		fid[i] = tpi*r[i]*rho[i] * (std::log(rho[i]*lam*lam*lam)-1.0);
		rho_x_rho_phi_int[i]  = tpi*r[i]*rho[i] * rho_phi_int[i];
		rho_x_phi_ext_mu[i] = tpi*r[i]*rho[i] * (phi_ext(r[i]) - mu);
		//std::cout << i << ", " << rho[i] << ", " << fid[i] << std::endl;
	}
	float spr2 = M_PI*(dr/2.0)*(dr/2.0)/(tpi*r[0]);
	omega1 = (kb1*T) * (integral_simpson(fid, nstep-1, dr) + fid[0]*spr2); // Fid
	omega2 = (kb1*T) * (integral_simpson(fex_i, nstep-1, dr) + fex_i[0]*spr2); // Fex
	omega3 = 0.5 * (integral_simpson(rho_x_rho_phi_int, nstep-1, dr) + rho_x_rho_phi_int[0]*spr2);
	omega4 = integral_simpson(rho_x_phi_ext_mu, nstep-1, dr) + rho_x_phi_ext_mu[0]*spr2;
	//
	omega_out = (omega1 + omega2 + omega3 + omega4) / epsilon_ff;
	//std::cout << omega1 << ", " << omega2 << ", " << omega3 << ", " << omega4 << std::endl;
	return omega_out;
}

int main(){
	int i,j,k;
	float v_gamma;
	float press_b, press_b0, pp0;
	float rho_b;
	float v_mmol_per_cm3;
	float v_cm3STP_per_g;
	float grand_potential;
	//
	read_parameters();
	float r[nstep];
	float rho[nstep], rho_new[nstep];
	//
	for (i=0; i<nstep; i++){
		r[i] = dr*(0.5+float(i)); // dr = (Dcc-sigma_ss)/float(nstep+1);
		//std::cout << i << ", " << r[i] << std::endl;
	}
	
	// show alpha
	//calc_alpha(r);
	// alpha = calc_alpha(r);
	
	// set rho_b0
	float y, a, b, c;
	float flag_P; flag_P = 0.0;
	float rho_b1; rho_b1 = 0.0;
	if ( rho_b0 == 0.0 ) {
		rho_b0 = Maxwell_construction();
	} else if ( rho_b0 < 0.0 ) {
		// rho_b0 < 0.0
		flag_P = -1.0;
		y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
		a = -0.5*alpha;
		b = kb1*T*(1.0 + y + y*y - y*y*y)/((1.0-y)*(1.0-y)*(1.0-y));
		c = -1.0*p0/(kb*1e27);
		rho_b1 = (-b+std::pow((b*b-4.0*a*c),0.5))/(2.0*a);
		if ( rho_b0==-10.0 ) {
			// change [Pa] to [atm]
			flag_P=-10.0;
		} else if ( rho_b0<=-100.0 ) {
			// change High pressure range to 1-100 [atm]
			flag_P=rho_b0;
			c = -1.0*101325.0*(flag_P*-1.0)/(kb*1e27);
		}
		rho_b0 = (-b+std::pow((b*b-4.0*a*c),0.5))/(2.0*a);
	}
	press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
	pp0 = press_b0*kb*1e27;
	std::cout << "rho_b0 = " << rho_b0 << std::endl;
	std::cout << "Pressure     : " << pp0 << " [Pa]" << std::endl;
	if ( flag_P>-100.0 ) {
		//std::cout << "Ref. Pressure: 101325 [Pa] = 1 [atm]" << std::endl;
		std::cout << "P0           : " << p0 << " [Pa] = " << p0/101325.0 << " [atm]" << std::endl;
	} else if ( flag_P<=-100.0 ) {
		std::cout << "Ref. Pressure: 1.01325e+07 [Pa] = 100 [atm] (10.1325 [MPa])" << std::endl;
	}
	
	//std::cout << rho_b0 << std::endl;
	// initialization
	for (i=0; i<nstep; i++){
		rho[i] = rho_b0/(nstep*dr);
		rho_new[i] = 0.0;
	}
	//
	//float rho_sj[nstep];
	//float rho_s0j[nstep];
	//float rho_s1j[nstep];
	//float rho_s2j[nstep];
	float n0_j[nstep], n0[nstep];   // For FMT
	float n1_j[nstep], n1[nstep];   // For FMT
	float n2_j[nstep], n2[nstep];   // For FMT
	float n3_j[nstep], n3[nstep];   // For FMT
	float nv1_j[nstep], nv1[nstep]; // For FMT
	float nv2_j[nstep], nv2[nstep]; // For FMT
	//float rho_dfex_int[nstep];
	//float dfex_int[nstep];
	float rho_phi_int[nstep];
	float phi_ext_i[nstep];
	for (i=0; i<nstep; i++){
		phi_ext_i[i] = phi_ext(r[i]);
		//std::cout << "phi_ext_i[" << i << "] = " << phi_ext_i[i] << std::endl;
	}
	float *phi_att_int_ij = (float *)malloc(sizeof(float)*((nstep+1)*nstep));
	if (phi_att_int_ij == NULL) {
		printf("Memory cannot be allocated.");
		std::exit(1);
	} else {
		printf("Memory has been allocated. The address is %p\n", phi_att_int_ij);
	}
	phi_att_int(r, phi_att_int_ij); // calculate integral phi_att at r[i]
	//
	float c1;
	float fex_i[nstep];  // For grand potential, Omega
	float rho_r[nstep];
	float check_c1xi;
	//
	float diff_old1 = 1.0;
	float diff;
	float diff0;
	float threshold = 0.5/100*nstep;
	//
	float xio;
	float rho_b_k[182]={3.91276e-08,7.56979e-08,1.42189e-07,2.59316e-07,4.59813e-07,
						7.65e-07,1.48e-06,2.78e-06,5.07e-06,8.99e-06,1.55e-05,2.61e-05,4.28e-05,6.87e-05,0.000107744,
						0.000165450,0.000249000,0.000367617,0.000532901,0.000759151,0.001063641,0.001466842,0.001992605,0.002668158,0.003524105,
						0.004594237,0.005915211,0.007526184,0.009468211,0.011783671,0.014515526,0.017706579,0.021398421,0.025631184,0.030442237,
						0.035865395,0.041930789,0.048663553,0.056083947,0.064206711,0.073040395,0.082588289,0.092847105,0.103808026,0.115456447,
						0.127772237,0.140730263,0.154301316,0.168452632,0.183144737,0.198336842,0.213988158,0.230053947,0.246484211,0.263235526,
						0.280259211,0.297506579,0.314930263,0.332486842,0.350127632,0.367810526,0.385494737,0.403138158,0.420705263,0.438159211,
						0.455467105,0.472598684,0.489525000,0.506219737,0.522661842,0.538827631,0.554700000,0.570261842,0.585498684,0.600400000,
						0.614953947,0.629153947,0.642990789,0.656461842,0.669563158,0.682292105,0.694648684,0.706632895,0.718247368,0.729494737,
						0.740377632,0.750900000,0.761068421,0.770886842,0.780363158,0.789501316,0.798311842,0.806798684,0.814971053,0.822838158,
						0.830405263,0.837682895,0.844677632,0.851397368,0.857852632,0.864050000,0.869998684,0.875705263,0.881178947,0.886427632,
						0.891459211,0.896281579,0.900901316,0.905326316,0.909563158,0.913619737,0.917503947,0.921219737,0.924775000,0.928177632,
						0.931430263,0.934542105,0.937517105,0.940361842,0.943080263,0.945678947,0.948161842,0.950534211,0.952801316,0.954965789,
						0.957034211,0.959009211,0.960896053,0.962697368,0.964417105,0.966059211,0.967626316,0.969122368,0.970550000,0.971913158,
						0.973214474,0.974455263,0.975640789,0.976771053,0.977848684,0.978877632,0.979859211,0.980796053,0.981689474,0.982542105,
						0.983353947,0.984128947,0.984868421,0.985573684,0.986247368,0.986888158,0.987500000,0.988082895,0.988639474,0.989169737,
						0.989676316,0.990157895,0.990618421,0.991056579,0.991475000,0.991873684,0.992253947,0.992615789,0.992961842,0.993290789,
						0.993603947,0.993903947,0.994189474,0.994461842,0.994721053,0.994968421,0.995203947,0.995428947,0.995642105,0.995846053,
						0.996040789,0.996226316,0.996403947,0.996572368,0.996732895,0.996885526,0.997031579};
	//
	// P/P0, V[molecules/nm^3], Omega/epsilon_ff[nm^-2]
	string Punit;
	string Punits;
	if(flag_P>=0.0){
		Punit = "P/P0";
		Punits = "PP0";
	} else if(flag_P<=-10.0){
		Punit = "atm";
		Punits= "atm";
	} else{
		Punit = "Pa";
		Punits= "Pa";
	}

	std::ofstream ofsppov_vs("./"+Punits+"_vs_Vgamma_data_vs.txt");
	ofsppov_vs << "# w = (2.0*Rcc-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_vs << "# P[" << Punit << "], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "w = (2.0*Rcc-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	std::cout << "P[" << Punit << "], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
	//
	for (k=0; k<=181; k++){
		//rho_b = rho_b0 * rho_b_k[k];
		if(flag_P<=-100.0){
			rho_b = (rho_b0 - rho_b1) * (rho_b_k[k] - 3.91276e-08) + rho_b1;
		} else {
			rho_b = rho_b0 * rho_b_k[k];
		}
		//
		for (j=0; j<cycle_max; j++){
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dr*nstep/2. 
			//rho_s(rho, r, rho_sj, rho_s0j, rho_s1j, rho_s2j);
			for (i=0; i<nstep; i++){
				ni(rho, r, i, n0_j, n1_j, n2_j, n3_j, nv1_j, nv2_j, n0, n1, n2, n3, nv1, nv2);
			}
			//for (i=0; i<=(nstep-2)/2; i++){
			for (i=0; i<nstep; i++){
				c1 = dfex(r, i, n0, n1, n2, n3, nv1, nv2);
				xio = c1+xi(rho,r,i,rho_b, phi_att_int_ij, rho_phi_int, phi_ext_i)/(kb1*T); // xi include kb1*T*(std::log(rho_b)) type.
				if (-14 < xio && xio < 12){
					rho_new[i] = std::exp(xio); // xi include kb1*T*(std::log(rho_b)) type.
				} else if (xio < -14){
					rho_new[i] = 1e-6;
				} else {
					// overflow about std::exp(730)
				    // to avoid overflow
					rho_new[i] = (2.0*rho_b0/dr + rho[i])*1.2;
				}
			}
			diff_old1 = diff;
			diff = 0.0;
			for (i=0; i<nstep; i++){
				diff0 = std::abs((rho_new[i]-rho[i])/rho[i]);
				diff = diff + diff0;
				rho[i] = wmixing*rho_new[i] + (1.0-wmixing)*rho[i];
			}
			//
			if (diff < threshold && diff_old1 < threshold && j>=500) {
				break;
			}
		}
		//for (i=0; i<nstep; i++){
		//	std::cout << "--------------------------------------------------" << std::endl;
		//	std::cout << "cycle=" << j << ", r[" << i << "]" << r[i] << " , -ext " << - phi_ext(r[i]) << ", rho[" << i << "]=" << rho[i] << std::endl;
		//}
		//
		for (i=0; i<nstep; i++){
			rho_r[i] = rho[i]*2.0*M_PI*r[i];
		}
		//integral_simpson(float *f, int n, float dx)
		v_gamma = integral_simpson(rho_r, nstep-1, dr) + rho[0]*M_PI*(dr/2.0)*(dr/2.0);
		//integral_trapezoidal(float *f, int n, float dx)
		//v_gamma = integral_trapezoidal(rho_r, nstep-1, dr) + rho[0]*M_PI*(dr/2.0)*(dr/2.0);
		v_gamma = v_gamma/(M_PI*(Dcc-sigma_ss)/2.0*(Dcc-sigma_ss)/2.0) - rho_b;
		//v_gamma = 2.0*v_gamma/(Dcc-sigma_ss) - rho_b*(Dref*Dref)/(4.0*(Dcc-sigma_ss));
		if (v_gamma < 0) { v_gamma = 0.0; }
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10.0; // [mmol/cm3]
		//v_cm3STP_per_g = v_mmol_per_cm3 * 22.414 / (rho_ss*12.0107*(1e7*1e7*1e7)/(6.02214076*1e23)); // [cm3(STP)/g], 2.226 [g/cm3]
		v_cm3STP_per_g = v_mmol_per_cm3 * 22.414 / (rho_ss*12.0107*10.0/6.02214076); // [cm3(STP)/g], 2.226 [g/cm3]
		//v_gamma = v_gamma * (0.8064/28.0134/1e21*6.02214e23)/rho_b;
		// N2(77K): 0.8064 g/mL, 0.8064/28.0134 mol/mL, 0.8064/28.0134/1e21 mol/nm3, 0.8064/28.0134/1e21*6.02214e23 molecules/nm3
		//std::cout << "V= " << v_gamma << std::endl;
		// press_hs(rho_b) from Carnahan-Starling (CS) equation of state
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		//
		if(flag_P==0.0){
			pp0 = press_b/press_b0;
		} else if (flag_P<=-10.0){
			pp0 = press_b*kb*1e27/p0;
		} else {
			// kb1=1, kb = 1.38e-23 [J/K], T [K], rho_b [N/nm^3], 1 [atm] = 101325 [Pa]
			pp0 = press_b*kb*1e27;
		}
		//
		// grand ppotential, Omega
		for (i=0; i<nstep; i++){
			fex_i[i] = fex(i, n0, n1, n2, n3, nv1, nv2);
		}
		grand_potential = omega(rho, r, fex_i, rho_phi_int, rho_b);
		//
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_vs << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_g << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_g << ", " << grand_potential << std::endl;
	}
	// reverse
	// P/P0, V[molecules/nm^3], Omega/epsilon_ff[nm^-2]
	std::ofstream ofsppov_ls("./"+Punits+"_vs_Vgamma_data_ls.txt");
	ofsppov_ls << "# w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_ls << "# P[" << Punit << "], V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	//std::cout << "w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	//std::cout << "P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
	for (k=181; k>=0; k--){
		//rho_b = rho_b0 * rho_b_k[k];
		if(flag_P<=-100.0){
			rho_b = (rho_b0 - rho_b1) * (rho_b_k[k] - 3.91276e-08) + rho_b1;
		} else {
			rho_b = rho_b0 * rho_b_k[k];
		}
		//
		for (j=0; j<cycle_max; j++){
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dr*nstep/2. 
			//rho_s(rho, r, rho_sj, rho_s0j, rho_s1j, rho_s2j);
			for (i=0; i<nstep; i++){
				ni(rho, r, i, n0_j, n1_j, n2_j, n3_j, nv1_j, nv2_j, n0, n1, n2, n3, nv1, nv2);
			}
			//for (i=0; i<=(nstep-2)/2; i++){
			for (i=0; i<nstep; i++){
				c1 = dfex(r, i, n0, n1, n2, n3, nv1, nv2);
				xio = c1+xi(rho,r,i,rho_b, phi_att_int_ij, rho_phi_int, phi_ext_i)/(kb1*T); // xi include kb1*T*(std::log(rho_b)) type.
				if (-14 < xio && xio < 12){
					rho_new[i] = std::exp(xio); // xi include kb1*T*(std::log(rho_b)) type.
				} else if (xio < -14){
					rho_new[i] = 1e-6;
				} else {
					// overflow about std::exp(730)
				    // to avoid overflow
					rho_new[i] = (2.0*rho_b0/dr + rho[i])*1.2;
				}
			}
			diff_old1 = diff;
			diff = 0.0;
			for (i=0; i<nstep; i++){
				diff0 = std::abs((rho_new[i]-rho[i])/rho[i]);
				diff = diff + diff0;
				rho[i] = wmixing*rho_new[i] + (1.0-wmixing)*rho[i];
			}
			//
			if (diff < threshold && diff_old1 < threshold && j>=500) {
				break;
			}
		}
		//for (i=0; i<nstep; i++){
		//	std::cout << "--------------------------------------------------" << std::endl;
		//	std::cout << "cycle=" << j << ", r[" << i << "]" << r[i] << " , -ext " << - phi_ext(r[i]) << ", rho[" << i << "]=" << rho[i] << std::endl;
		//}
		//
		for (i=0; i<nstep; i++){
			rho_r[i] = rho[i]*2.0*M_PI*r[i];
		}
		//integral_simpson(float *f, int n, float dx)
		v_gamma = integral_simpson(rho_r, nstep-1, dr) + rho[0]*M_PI*(dr/2.0)*(dr/2.0);
		//integral_trapezoidal(float *f, int n, float dx)
		//v_gamma = integral_trapezoidal(rho_r, nstep-1, dr) + rho[0]*M_PI*(dr/2.0)*(dr/2.0);
		v_gamma = v_gamma/(M_PI*(Dcc-sigma_ss)/2.0*(Dcc-sigma_ss)/2.0) - rho_b;
		if (v_gamma < 0) { v_gamma = 0.0; }
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10.0; // [mmol/cm3]
		//v_cm3STP_per_g = v_mmol_per_cm3 * 22.414 / (rho_ss*12.0107*(1e7*1e7*1e7)/(6.02214076*1e23)); // [cm3(STP)/g], 2.226 [g/cm3]
		v_cm3STP_per_g = v_mmol_per_cm3 * 22.414 / (rho_ss*12.0107*10.0/6.02214076); // [cm3(STP)/g], 2.226 [g/cm3]
		//v_gamma = v_gamma * (0.8064/28.0134/1e21*6.02214e23)/rho_b;
		// N2(77K): 0.8064 g/mL, 0.8064/28.0134 mol/mL, 0.8064/28.0134/1e21 mol/nm3, 0.8064/28.0134/1e21*6.02214e23 molecules/nm3
		//std::cout << "V= " << v_gamma << std::endl;
		// press_hs(rho_b) from Carnahan-Starling (CS) equation of state for SDA
		// press_hs(rho_b) from Percus-Yevick (PY) equation of state for FMT
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		//
		if(flag_P==0.0){
			pp0 = press_b/press_b0;
		} else if (flag_P<=-10.0){
			pp0 = press_b*kb*1e27/p0;
		} else {
			// kb1=1, kb = 1.38e-23 [J/K], T [K], rho_b [N/nm^3], 1 [atm] = 101325 [Pa]
			pp0 = press_b*kb*1e27;
		}
		//
		// grand ppotential, Omega
		for (i=0; i<nstep; i++){
			fex_i[i] = fex(i, n0, n1, n2, n3, nv1, nv2);
		}
		grand_potential = omega(rho, r, fex_i, rho_phi_int, rho_b);
		//
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_ls << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_g << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_g << ", " << grand_potential << std::endl;
	}
	free(phi_att_int_ij);
	return 0;
}
