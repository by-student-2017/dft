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

// compiling: c++ nldft_sda_slit.cpp -O2
// usage: ./a.out

// debag mode
// compiling: c++ nldft_sda_slit.cpp -g -Wall -O0
// run: gdb ./a.out
//      (gdb) run

// ---------- ----------- ------------ ------------
// Adsorbent
float H;    //distace of slit [nm]
float sigma_ss; // sigma_ss = 0.34 [nm]
int nstep;   // number of mesh on z axis
float w_pw; // w_pw = (H-sigma_ss), pore width [nm]
float dr;   // dr = (H-sigma_ss)/float(nstep-1)
// ---------- ----------- ------------ ------------
// assume rho is same value in x-y plane.
// cylinder and normalization, because of cut off (rc).
int nrmesh; //rho_si and xi function
float drc; //float drc = rc/float(nrmesh-1);
// ---------- ----------- ------------ ------------
// iteration of rho
int cycle_max;  //int cycle_max = 50;
float wmixing; //float wmixing = 0.005;
// ---------- ----------- ------------ ------------
//Carbon dioxide 253.9  [K](epsilon), 0.3454 [nm](sigma), 0.3495 [nm](d_hs)
//Argon          118.05 [K](epsilon), 0.3305 [nm](sigma), 0.3390 [nm](d_hs)
//Nitrogen        94.45 [K](epsilon), 0.3575 [nm](sigma), 0.3575 [nm](d_hs), P0=4.97528=101325 Pa
//Hydrogen        34.3  [K](epsilon), 0.3040 [nm](sigma), 0.3040 [nm](d_hs) (J. Jagiello 2007)
//extern float epsilon_ff = 94.45;
//float sigma_ff = 0.3575;
//extern float d_hs = 0.3575; // Maxwell_construction()
//float rc = 1.28; // [nm],cut off, (12.8 [A])
float epsilon_ff;
float sigma_ff;
float d_hs;
float rc;
// ---------- ----------- ------------ ------------
float rm; // rm = std::pow(2.0,1.0/6.0)*sigma_ff = 1.12246205*sigma_ff; //minimum position of LJ
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

//Barker-Henderson (BH) theory
float d_bh_calc(float epsilon, float sigma){
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

void read_parameters(void){
	std::ifstream ifs("parameters.txt");
	std::string str;
	//
	int i,j;
	float num[20];
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
		nstep = int((H-sigma_ss)/0.01 + 0.5);
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
	nrmesh = int(num[9]);
	if ( nrmesh == 0 ) {
		nrmesh = int(rc/0.04 + 0.5);
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
	rho_b0 = num[16];
	// ---------- ----------- ------------ ------------
	
	w_pw = (H-sigma_ss); // pore width [nm]
	dr = (H-sigma_ss)/float(nstep-1);
	drc = rc/float(nrmesh-1); // rho_si(), calc_alpha(), xi()
	rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
	
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

float rho_si_int_k(float *r, float *rho_si_int_ijrj){
	int i;
	int jr;
	int j,k;
	float ra;
	float raj;
	float rak;
	//
	//float rho_si_int_j[nstep];
	float tmp_rho_si_int_k[nrmesh];
	//
	// drc = rc/float(nrmesh-1);
	float tpidrc = 2.0*M_PI*drc;
	//
	for (i=0; i<3; i++) {
		for (jr=0; jr<nstep; jr++) {
			//
			for (j=0; j<nstep; j++) {
				raj = (r[j]-r[jr]);
				tmp_rho_si_int_k[0] = 0.0;
				for (k=1; k<nrmesh; k++) {
					rak = drc*float(k);
					ra =  rak*rak + raj*raj;
					ra = std::sqrt(ra);
					//
					tmp_rho_si_int_k[k] = wi(ra,i)*(tpidrc*float(k));
				}
				//integral_simpson(float *f, int n, float dx)
				rho_si_int_ijrj[i*nstep*nstep+jr*nstep+j] = integral_simpson(tmp_rho_si_int_k, nrmesh-1, drc);
			}
			//
		}
	}
	return 0;
}

// Tarazona theory
float rho_si(float *rho, float *r, int jr, int i, float *rho_si_int_ijrj){
	int j,k;
	float ra;
	float raj;
	float rak;
	//
	float rho_si_int_j[nstep];
	//float rho_si_int_k[nrmesh];
	//
	// drc = rc/float(nrmesh-1);
	float tpidrc = 2.0*M_PI*drc;
	//rho_si_int_k[0] = 0.0;
	//
	for (j=0; j<nstep; j++) {
		//raj = (r[j]-r[jr]);
		//for (k=1; k<nrmesh; k++) {
		//	rak = drc*float(k);
		//	ra =  rak*rak + raj*raj;
		//	ra = std::sqrt(ra);
		//	//
		//	rho_si_int_k[k] = wi(ra,i)*(tpidrc*float(k));
		//}
		//integral_simpson(float *f, int n, float dx)
		//rho_si_int_j[j] = rho[j]*integral_simpson(rho_si_int_k, nrmesh-1, drc);
		rho_si_int_j[j] = rho[j] * rho_si_int_ijrj[i*nstep*nstep+jr*nstep+j];
	}
	float rho_si_out;
	rho_si_out = integral_simpson(rho_si_int_j, nstep-1, dr);
	//
	return rho_si_out;
}

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
float rho_s(float *rho, float *r, float *rho_sj, float *rho_s0j, float *rho_s1j, float *rho_s2j, float *rho_si_int_ijrj){
	int j;
	float rho_den1j, rho_den2j;
	//for (j=0; j<(nstep-2)/2; j++) {
	for (j=0; j<nstep; j++) {
		rho_s0j[j] = rho_si(rho, r, j, 0, rho_si_int_ijrj);
		rho_s1j[j] = rho_si(rho, r, j, 1, rho_si_int_ijrj);
		rho_s2j[j] = rho_si(rho, r, j, 2, rho_si_int_ijrj);
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
		//rho_s0j[(nstep-1)-j] = rho_s0j[j];
		//rho_s1j[(nstep-1)-j] = rho_s1j[j];
		//rho_s2j[(nstep-1)-j] = rho_s2j[j];
		//rho_sj[(nstep-1)-j] = rho_sj[j];
	}
	return 0;
}

// Steele 10-4-3 potential
float phi_sf(float z){
	float phi_sf_out;
	float sigma_sf2 = sigma_sf*sigma_sf;
	float sfpz = (sigma_sf/z);
	float sfpz2 = sfpz*sfpz;
	float dez = (0.61*delta+z);
	//phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*std::pow(sigma_sf,2.0)*delta*
	//			( (2.0/5.0)*std::pow((sigma_sf/z),10.0)-std::pow((sigma_sf/z),4.0)-std::pow(sigma_sf,4.0)/
	//			(3.0*delta*std::pow((0.61*delta+z),3.0)) );
	phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*(sigma_sf2)*delta*
				( (2.0/5.0)*std::pow(sfpz2,5.0)-(sfpz2*sfpz2)-(sigma_sf2*sigma_sf2)/
				(3.0*delta*(dez*dez*dez)) );
	return phi_sf_out;
}

// e.g., wall potential (Carbon slit)
float phi_ext(float z){
	float phi_ext_out;
	phi_ext_out = phi_sf(z) + phi_sf(H-z);
	//std::cout << phi_ext_out << std::endl;
	return phi_ext_out;
}

// from Carnahan-Starling (CS) equation of state
float mu_ex(float rho_b){
	float y, mu_ex_out;
	//y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
	float den1y = (1.0-y);
	//mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/std::pow((1.0-y),3.0);
	//mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/((1.0-y)*(1.0-y)*(1.0-y));
	mu_ex_out = kb1*T*(8.0*y-9.0*y*y+3.0*y*y*y)/(den1y*den1y*den1y);
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

float f_ex(float rho_s){
	float eta, f_ex_out;
	//eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
	eta = M_PI*rho_s*(d_hs*d_hs*d_hs)/6.0;
	float den1e = (1.0-eta);
	//f_ex_out = kb1*T*eta*(4.0-3.0*eta)/std::pow((1.0-eta),2.0);
	//f_ex_out = kb1*T*eta*(4.0-3.0*eta)/((1.0-eta)*(1.0-eta));
	f_ex_out = kb1*T*eta*(4.0-3.0*eta)/(den1e*den1e);
	return f_ex_out;
}

// d(f_ex)/d(rho_s)
float dfex_per_drhos(float rho_s){
	float dfex_per_drhos_out;
	float eta;
	//eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
	eta = M_PI*rho_s*(d_hs*d_hs*d_hs)/6.0;
	float den1e = (1.0-eta);
	//dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/std::pow((1.0-eta),3.0)*M_PI*std::pow(d_hs,3.0)/6.0;
	//dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/((1.0-eta)*(1.0-eta)*(1.0-eta))*M_PI*(d_hs*d_hs*d_hs)/6.0;
	dfex_per_drhos_out = kb1*T*(4.0-2.0*eta)/(den1e*den1e*den1e)*(M_PI*(d_hs*d_hs*d_hs)/6.0);
	return dfex_per_drhos_out;
}

// d(rho_s)/d(rho)
//float drhos_per_drho(float *rho, float r1, float r2, float *r, float ra){
//	float w, drhos_per_drho_out;
//	// Percus-Yevick approximation, Tarazona theory
//	w = wi(ra,0) + wi(ra,1)*rho_s(rho,r1,r) + wi(ra,2)*std::pow(rho_s(rho,r1,r),2.0);
//	drhos_per_drho_out = w/(1.0-rho_si(rho,r2,r,1)-2.0*rho_si(rho,r2,r,2)*rho_s(rho,r2,r));
//	return drhos_per_drho_out;
//}

// d(rho_s)/d(rho), modified version
float drhos_per_drho_j(float ra, float rho_sj, float rho_s1j, float rho_s2j){
	float w, drhos_per_drho_out;
	// Percus-Yevick approximation, Tarazona theory
	//w = wi(ra,0) + wi(ra,1)*rho_sj + wi(ra,2)*std::pow(rho_sj,2.0);
	w = wi(ra,0) + wi(ra,1)*rho_sj + wi(ra,2)*(rho_sj*rho_sj);
	drhos_per_drho_out = w/(1.0-rho_s1j-2.0*rho_s2j*rho_sj);
	return drhos_per_drho_out;
}

float calc_alpha(float *r){
	int i,j,k;
	float ra;
	float raj;
	float rak;
	//
	float alpha_other_method;
	float alpha_int_j[nstep];
	float alpha_int_k[nrmesh];
	//
	//float drc = rc/float(nrmesh-1);
	float tpidrc = 2.0*M_PI*drc;
	alpha_int_k[0] = 0.0;
	//
	for (i=0; i<=(nstep-2)/2; i++){
		for (j=0; j<nstep; j++) {
			raj = (r[i]-r[j]);
			for (k=1; k<nrmesh; k++) {
				rak = drc*float(k);
				ra = raj*raj + rak*rak;
				ra = std::sqrt(ra);
				//
				alpha_int_k[k]  = -phi_att(ra)*(tpidrc*float(k));
			}
			//integral_simpson(float *f, int n, float dx)
			alpha_int_j[j]  = integral_simpson(alpha_int_k, nrmesh-1, drc);
		}
		alpha_other_method  = alpha_other_method + integral_simpson(alpha_int_j, nstep-1, dr);
	}
	alpha_other_method  = alpha_other_method * 2.0 * dr / (H-sigma_ss);
	//std::cout << "--------------------------------------------------" << std::endl;
	//std::cout << "average alpha of other method = " << alpha_other_method << " in (carbon) slit" << std::endl;
	return alpha_other_method;
}

float phi_att_int(float *r, float *phi_att_int_ij){
	int i,j,k;
	float ra;
	float raj;
	float rak;
	//
	float phi_int_k[nrmesh];
	//
	//float drc = rc/float(nrmesh-1);
	float tpidrc = 2.0*M_PI*drc;
	phi_int_k[0] = 0.0;
	//
	for (i=0; i<nstep; i++) {
		for (j=0; j<nstep; j++) {
			raj = (r[i]-r[j]);
			for (k=1; k<nrmesh; k++) {
				rak = drc*float(k);
				ra =  rak*rak + raj*raj;
				ra = std::sqrt(ra);
				//
				//phi_int_k[k]  = phi_att(ra)*(2.0*M_PI*(float(k)*drc));
				phi_int_k[k]  = phi_att(ra)*(tpidrc*float(k));
			}
			phi_att_int_ij[i*nstep+j] = integral_simpson(phi_int_k, nrmesh-1, drc);
		}
	}
	return 0;
}

// xi include kb1*T*(std::log(rho_b)) type.
// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
float xi(float *rho, float *r, int i, float rho_b, float *rho_sj, float *rho_s0j, float *rho_s1j, float *rho_s2j, float *phi_att_int_ij, float *rho_dfex_int, float *rho_phi_int, float *phi_ext_i){
	int j,k;
	float ra;
	float raj;
	float rak;
	//
	float rho_phi_int_j[nstep];
	float rho_dfex_int_j[nstep];
	float rho_dfex_int_k[nrmesh];
	//
	// drc = rc/float(nrmesh-1);
	float tpidrc = 2.0*M_PI*drc;
	rho_dfex_int_k[0] = 0.0;
	//
	for (j=0; j<nstep; j++) {
		raj = (r[i]-r[j]);
		for (k=1; k<nrmesh; k++) {
			rak = drc*float(k);
			ra =  rak*rak + raj*raj;
			ra = std::sqrt(ra);
			//
			// d(f_ex)/d(rho) = d(f_ex)/d(rho_s) * d(rho_s)/d(rho)
			rho_dfex_int_k[k] = drhos_per_drho_j(ra, rho_sj[j], rho_s1j[j], rho_s2j[j])*(tpidrc*float(k));
		}
		//integral_simpson(float *f, int n, float dx)
		rho_dfex_int_j[j] = rho[j]*dfex_per_drhos(rho_sj[j])*integral_simpson(rho_dfex_int_k, nrmesh-1, drc);
		rho_phi_int_j[j]  = rho[j]*phi_att_int_ij[i*nstep+j];
	}
	//integral_simpson(float *f, int n, float dx)
	rho_dfex_int[i] = integral_simpson(rho_dfex_int_j, nstep-1, dr);
	rho_phi_int[i]  = integral_simpson(rho_phi_int_j, nstep-1, dr);
	//
	float xi_out;
	xi_out = ( - rho_b*alpha - rho_dfex_int[i] - f_ex(rho_sj[i]) ) + ( mu_ex(rho_b) - rho_phi_int[i] ) + ( kb1*T*std::log(rho_b) - phi_ext_i[i] );
	//
	return xi_out;
}

float press_hs(float rho_b){
	float y, press_hs_out;
	//y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
	float den1y = (1.0-y);
	//press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/std::pow((1.0-y),3.0);
	//press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/((1.0-y)*(1.0-y)*(1.0-y));
	press_hs_out = rho_b*kb1*T* (1.0 + y + y*y - y*y*y)/(den1y*den1y*den1y);
	return press_hs_out;
}

float Maxwell_construction(void){
	int i,j;
	int iter_max_drhob0 = 500000;
	int iter_max_dmue = 50000;
	float drhob0 = 0.00005;
	float dmue = 0.0005;
	float threshold_diff = 0.03;
	float threshold_find = 0.03;
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
	std::cout << "Bulk pressure, P0 = " << press_b0 << " (rho_b0 = " << rho_b0_out << ")" <<std::endl;
	std::cout << std::endl;
	std::cout << "gas phase   : rho_b0_gas        = " << rho_b0_gas        << ", rho_b0_gas*d_hs^3        = " << rho_b0_gas*std::pow(d_hs,3.0)        << std::endl;
	std::cout << "metastable  : rho_b0_metastable = " << rho_b0_metastable << ", rho_b0_metastable*d_hs^3 = " << rho_b0_metastable*std::pow(d_hs,3.0) << std::endl;
	std::cout << "liquid phase: rho_b0_liquid     = " << rho_b0_liquid     << ", rho_b0_liquid*d_hs^3     = " << rho_b0_liquid*std::pow(d_hs,3.0)     << std::endl;
	return rho_b0_out;
}

// grand potential
float omega(float *rho, float *r, float *rho_dfex_int, float *rho_phi_int){
	float omega_out;
	float omega1, omega2, omega3;
	int i;
	int omega_nstep = (nstep-2)/2;
	float rho_x_rho_dfex_int[omega_nstep+1];
	float rho_x_rho_phi_int[omega_nstep+1];
	for (i=0; i<=omega_nstep; i++){
		rho_x_rho_dfex_int[i] = rho[i] * rho_dfex_int[i];
		rho_x_rho_phi_int[i]  = rho[i] * rho_phi_int[i];
	}
	omega1 = -(kb1*T) * integral_simpson(rho, omega_nstep, dr);
	omega2 = -integral_simpson(rho_x_rho_dfex_int, omega_nstep, dr);
	omega3 = -0.5 * integral_simpson(rho_x_rho_phi_int, omega_nstep, dr);
	omega_out = (omega1 + omega2 + omega3) * 2.0 / epsilon_ff;
	return omega_out;
}

int main(){
	int i,j,k;
	float v_gamma;
	float press_b, press_b0, pp0;
	float rho_b;
	float v_mmol_per_cm3;
	float v_cm3STP_per_cm3;
	float grand_potential;
	//
	read_parameters();
	float r[nstep];
	float rho[nstep], rho_new[nstep];
	//
	for (i=0; i<nstep; i++){
		r[i] = sigma_ss/2.0 + dr*float(i); // dr = (H-sigma_ss)/float(nstep+1);
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
	
	std::cout << "--------------------------------------------------" << std::endl;
	float rho_sj[nstep];
	float rho_s0j[nstep];
	float rho_s1j[nstep];
	float rho_s2j[nstep];
	//
	float rho_dfex_int[nstep];
	float rho_phi_int[nstep];
	float phi_ext_i[nstep];
	for (i=0; i<nstep; i++){
		phi_ext_i[i] = phi_ext(r[i]);
	}
	std::cout << "phi_ext_i calculation was finished" << std::endl;
	//float phi_att_int_ij[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	float *phi_att_int_ij = (float *)malloc(sizeof(float)*((nstep+1)*nstep));
	if (phi_att_int_ij == NULL) {
		printf("Memory cannot be allocated.");
		std::exit(1);
	} else {
		printf("Memory has been allocated. The address is %p\n", phi_att_int_ij);
	}
	phi_att_int(r, phi_att_int_ij); // calculate integral phi_att at r[i]
	std::cout << "phi_att_int calculation was finished" << std::endl;
	//
	float *rho_si_int_ijrj = (float *)malloc(sizeof(float)*(2*nstep*nstep+nstep*nstep+nstep));
	rho_si_int_k(r, rho_si_int_ijrj);
	std::cout << "rho_si_int_k calculation was finished" << std::endl;
	//
	//float diff_old2 = 1.0;
	float diff_old1 = 1.0;
	float diff;
	float diff0;
	float mixing;
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
	std::ofstream ofsppov_vs("./PP0_vs_Vgamma_data_vs.txt");
	ofsppov_vs << "# w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_vs << "# P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	std::cout << "P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega/epsilon_ff[1/nm2]" << std::endl;
	for (k=0; k<=181; k++){
		rho_b = rho_b0 * rho_b_k[k];
		// Hill Equation
		//rho_b = rho_b0 * (0.0 + (1.0 - -0.0))*
		//	(std::pow(float(k),4.2323)/(std::pow(float(k),4.2323)+std::pow(62.997,4.2323)));
	//for (k=0; k<100; k++){
		//rho_b = rho_b0 * std::exp(-(20.0-2.0*float(k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		for (j=0; j<cycle_max; j++){
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dr*nstep/2. 
			rho_s(rho, r, rho_sj, rho_s0j, rho_s1j, rho_s2j, rho_si_int_ijrj);
			for (i=0; i<=(nstep-2)/2; i++){
				//rho_new[i] = rho_b*std::exp(xi(rho,r[i],rho_b,r)/(kb1*T)); // this equation occure inf.
				//rho_new[i] = std::exp(xi(rho,r,i,rho_b, rho_sj, rho_s0j, rho_s1j, rho_s2j, phi_att_int_ij, rho_dfex_int, rho_phi_int, phi_ext_i)/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				xio = xi(rho,r,i,rho_b, rho_sj, rho_s0j, rho_s1j, rho_s2j, phi_att_int_ij, rho_dfex_int, rho_phi_int, phi_ext_i)/(kb1*T);
				if (-14 < xio && xio < 12){
					rho_new[i] = std::exp(xio); // xi include kb1*T*(std::log(rho_b)) type.
				} else if (xio < -14){
					rho_new[i] = 1e-6;
				} else {
					// overflow about std::exp(730)
				    // to avoid overflow
					rho_new[i] = rho[i] / 10.0;
				}
			}
			//diff_old2 = diff_old1;
			diff_old1 = diff;
			diff = 0.0;
			for (i=0; i<=(nstep-2)/2; i++){
				diff0 = std::abs((rho_new[i]-rho[i])/rho[i]);
				diff = diff + 2.0*diff0;
				//mixing = wmixing + wmixing/(0.5+diff0);
				//std::cout << i << ", " << mixing << std::endl;
				rho[i] = wmixing*rho_new[i] + (1.0-wmixing)*rho[i];
				rho[(nstep-1)-i] = rho[i]; // The rest is filled with mirror symmetry. 
			}
			//if (diff/nstep < 0.005 && diff_old/nstep < 0.005 && j >= 20) {
			//float threshold = 0.5/100*nstep;
			//if (diff < threshold && diff_old1 < threshold && diff_old2 < threshold) {
			if (diff < threshold && diff_old1 < threshold) {
				break;
			}
		}
		//
		v_gamma = integral_simpson(rho, nstep-1, dr);
		v_gamma = v_gamma/(H-sigma_ss) - rho_b;
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10.0; // [mmol/cm3]
		v_cm3STP_per_cm3 = v_mmol_per_cm3 * 22.414;
		if (v_gamma < 0) { v_gamma = 0.0; }
		//v_gamma = v_gamma * (0.8064/28.0134/1e21*6.02214e23)/rho_b;
		// N2(77K): 0.8064 g/mL, 0.8064/28.0134 mol/mL, 0.8064/28.0134/1e21 mol/nm3, 0.8064/28.0134/1e21*6.02214e23 molecules/nm3
		//std::cout << "V= " << v_gamma << std::endl;
		//
		// press_hs(rho_b) from Carnahan-Starling (CS) equation of state
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		//
		pp0 = press_b/press_b0;
		grand_potential = omega(rho, r, rho_dfex_int, rho_phi_int);
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_vs << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
	}
	// reverse
	std::ofstream ofsppov_ls("./PP0_vs_Vgamma_data_ls.txt");
	ofsppov_ls << "# w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_ls << "# P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	//std::cout << "w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	//std::cout << "P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega/epsilon_ff[1/nm2]" << std::endl;
	for (k=181; k>=0; k--){
		rho_b = rho_b0 * rho_b_k[k];
		// Hill Equation
		//rho_b = rho_b0 * (0.0 + (1.0 - -0.0))*
		//	(std::pow(float(k),4.2323)/(std::pow(float(k),4.2323)+std::pow(62.997,4.2323)));
	//for (k=0; k<100; k++){
		//rho_b = rho_b0 * std::exp(-(20.0-2.0*float(99.0-k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		for (j=0; j<cycle_max; j++){
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dr*nstep/2. 
			rho_s(rho, r, rho_sj, rho_s0j, rho_s1j, rho_s2j, rho_si_int_ijrj);
			for (i=0; i<=(nstep-2)/2; i++){
				//rho_new[i] = rho_b*std::exp(xi(rho,r[i],rho_b,r)/(kb1*T)); // this equation occure inf.
				//rho_new[i] = std::exp(xi(rho,r,i,rho_b, rho_sj, rho_s0j, rho_s1j, rho_s2j, phi_att_int_ij, rho_dfex_int, rho_phi_int, phi_ext_i)/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				xio = xi(rho,r,i,rho_b, rho_sj, rho_s0j, rho_s1j, rho_s2j, phi_att_int_ij, rho_dfex_int, rho_phi_int, phi_ext_i)/(kb1*T);
				if (-14 < xio && xio < 12){
					rho_new[i] = std::exp(xio); // xi include kb1*T*(std::log(rho_b)) type.
				} else if (xio < -14){
					rho_new[i] = 1e-6;
				} else {
					// overflow about std::exp(730)
				    // to avoid overflow
					rho_new[i] = rho[i] / 10.0;
				}
			}
			//diff_old2 = diff_old1;
			diff_old1 = diff;
			diff = 0.0;
			for (i=0; i<=(nstep-2)/2; i++){
				diff0 = std::abs((rho_new[i]-rho[i])/rho[i]);
				diff = diff + 2.0*diff0;
				//mixing = wmixing + wmixing/(0.5+diff0);
				//std::cout << i << ", " << mixing << std::endl;
				rho[i] = wmixing*rho_new[i] + (1.0-wmixing)*rho[i];
				rho[(nstep-1)-i] = rho[i]; // The rest is filled with mirror symmetry. 
			}
			//if (diff/nstep < 0.005 && diff_old/nstep < 0.005 && j >= 20) {
			//float threshold = 0.5/100*nstep;
			//if (diff < threshold && diff_old1 < threshold && diff_old2 < threshold) {
			if (diff < threshold && diff_old1 < threshold) {
				break;
			}
		}
		//
		v_gamma = integral_simpson(rho, nstep-1, dr);
		v_gamma = v_gamma/(H-sigma_ss) - rho_b;
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10.0; // [mmol/cm3]
		v_cm3STP_per_cm3 = v_mmol_per_cm3 * 22.414;
		if (v_gamma < 0) { v_gamma = 0.0; }
		//v_gamma = v_gamma * (0.8064/28.0134/1e21*6.02214e23)/rho_b;
		// N2(77K): 0.8064 g/mL, 0.8064/28.0134 mol/mL, 0.8064/28.0134/1e21 mol/nm3, 0.8064/28.0134/1e21*6.02214e23 molecules/nm3
		//std::cout << "V= " << v_gamma << std::endl;
		//
		// press_hs(rho_b) from Carnahan-Starling (CS) equation of state
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		//
		pp0 = press_b/press_b0;
		grand_potential = omega(rho, r, rho_dfex_int, rho_phi_int);
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_ls << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
	}
	free(phi_att_int_ij);
	free(rho_si_int_ijrj);
	return 0;
}
