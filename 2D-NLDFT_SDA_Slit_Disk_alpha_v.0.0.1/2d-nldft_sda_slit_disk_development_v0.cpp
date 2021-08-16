#include <fstream>   // for file in and out
#include <iostream>  // for cout
#include <cmath>     // for log, exp
#include <sstream>   // for read parameters

//#include "maxwell_construction.h"

using namespace std;

//non-local smoothed density approximation：SDA
//The two dimentional non-local density functional theory（2D-NLDFT)
//Reference: https://www.j-ad.org/adsorption_news/30_1.pdf

// Note
// Units are fundamentaly [K] and [nm] in this routine.

// There are many imperfections, so I hope someone can make it better with a CC0 license. 
// It seems that this code is the first in the world at present (2021/7/5) to be released on CC0 even in NLDFT. 

// compiling: c++ nldft.cpp -O2
// usage: ./a.out

// debag mode
// compiling: c++ nldft.cpp -g -Wall -O0
// run: gdb ./a.out
//      (gdb) run

// ---------- ----------- ------------ ------------
// Adsorbent 
double H; // distace of slit. [nm]
double sigma_ss; // [nm]
int nzstep; // number of step between disks
int nxstep; // number of step
double w_pw; // pore width, [nm]
double dz; //
double dx; //
double D;  // diameter of disk
// ---------- ----------- ------------ ------------
// assume rho is same value in x-y plane.
// cylinder and normalization, because of cut off (rc).
int nrmesh; // number of mesh, radius of disk
int ntmesh; // theta
// ---------- ----------- ------------ ------------
// iteration of rho
int cycle_max;
double wmixing; // mixing parameter
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
double rc;   // fluid-fluid
// ---------- ----------- ------------ ------------
double rm; // rm = std::pow(2.0,1.0/6.0)*sigma_ff = 1.12246205*sigma_ff, minimum position of LJ
double rmsf; // rmsf = 1.12246205*sigma_sf; // 2^(1/6)=1.12246205
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
double h0;
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

double rho_ssq(double z){
	double rho_ssq_out;
	//double h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	//double rho_ss = 114.0; // [molecules/nm3] (the density of bulk carbon)
	//double delta = 0.13; // [nm] (the roughness parameter) (the half-width of the density ramp)
	if ( 0.0 <= z && z < h0 ){
		rho_ssq_out = rho_ss;
	} else {
		rho_ssq_out = 0.0;
	}
	return rho_ssq_out;
}

void read_parameters(void){
	std::ifstream ifs("parameters.txt");
	std::string str;
	//
	int i,j;
	double num[25];
	j = 0;
	//
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
	nzstep = int(num[2]);
	if ( nzstep == 0 ) {
		nzstep = int((H-sigma_ss)/0.02 + 0.5);
		if ( nzstep%2 == 1 ){
			nzstep = nzstep + 1;
		}
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "autoset nzstep = " << nzstep << std::endl;
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
	ntmesh = int(num[9]);
	if ( ntmesh == 0 ) {
		ntmesh = 20;
		std::cout << "autoset ntmesh = " << ntmesh << std::endl;
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
	h0 = num[17];
	// ---------- ----------- ------------ ------------
	double DH;
	DH = num[18]; // D/H ratio
	D = DH * H;
	// ---------- ----------- ------------ ------------
	nxstep = int(num[19]);
	if ( nxstep == 0 ) {
		nxstep = int((D/2.0)/0.02 + 0.5);
		if ( nxstep%2 == 1 ){
			nxstep = nxstep + 1;
		}
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "autoset nxstep = " << nxstep << std::endl;
	}
	// ---------- ----------- ------------ ------------
	
	w_pw = (H-sigma_ss); // pore width [nm]
	dz = (H-sigma_ss)/double(nzstep-1);
	dx = (D/2.0)/double(nxstep); // x-y plane
	rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
	//rmsf = 1.12246205*sigma_sf; // 2^(1/6)=1.12246205
	
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

double integral_trapezoidal(double *f, int n, double dq){
	double sum;
	sum = 0.0;
	int i;
	for(i=1; i<n; i++){
		sum += (f[i-1]+f[i])/2.0*dq;
	}
	return sum;
}

double integral_simpson(double *f, int n, double dq){
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
	return (dq/3.0)*sum;
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
// In this case, WCA type is different from result of steele potential.
double phi_att_sf(double r){
	double e;
	// Lennard-Jones（LJ) potential
	e = std::pow((sigma_sf/r),6.0);
	e = 4.0*epsilon_sf*( e*e - e );
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
double rho_si(double *rho, double x0, double z0, double *x, double *z, int i){
	int ix; // x axis for rho, x = y
	int iz; // z axis for rho
	int t;  // theta
	double raiz; // z axis
	double ra;
	//
	double rho_si_out;
	double rho_si_int_iz[nzstep];
	double rho_si_int_ix[nxstep];
	double rho_si_int_t[ntmesh];
	//
	double xt,yt;
	double drad = M_PI/nrmesh; // radian
	double spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	//
	for (iz=0; iz<nzstep; iz++) {
		raiz = (z[iz]-z0);
		for (ix=0; ix<nxstep; ix++) {
			for (t=0; t<ntmesh; t++) {
				xt = x[ix]*std::cos(drad*double(t));
				yt = x[ix]*std::sin(drad*double(t));
				ra= (xt-x0)*(xt-x0) + yt*yt + raiz*raiz; // x, y, z
				ra = std::sqrt(ra);
				rho_si_int_t[t]  = wi(ra,i);
			}
			//integral_simpson(double *f, int n, double dx)
			rho_si_int_ix[ix] = 2.0*x[ix]*rho[ix*nzstep+iz] * integral_simpson(rho_si_int_t, ntmesh-1, drad);
		}
		//integral_simpson(double *f, int n, double dx)
		//rho_si_int_iz[iz] = integral_simpson(rho_si_int_ix, nxstep-1, dx) + rho_si_int_ix[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
		rho_si_int_iz[iz] = integral_simpson(rho_si_int_ix, nxstep-1, dx) + rho_si_int_ix[0]*spr2;
	}
	//integral_simpson(double *f, int n, double dx)
	rho_si_out = integral_simpson(rho_si_int_iz, nzstep-1, dz);
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
double rho_s(double *rho, double *x, double *z, double *rho_s_ixiz, double *rho_s0_ixiz, double *rho_s1_ixiz, double *rho_s2_ixiz){
	int ix; // x axis for rho, x = y
	int iz; // z axis for rho
	int t;  // theta
	//
	double rho_den1;
	double rho_den2;
	for (ix=0; ix<nxstep; ix++) {
		for (iz=0; iz<nzstep; iz++) {
			rho_s0_ixiz[ix*nzstep+iz] = rho_si(rho, x[ix], z[iz], x, z, 0);
			rho_s1_ixiz[ix*nzstep+iz] = rho_si(rho, x[ix], z[iz], x, z, 1);
			rho_s2_ixiz[ix*nzstep+iz] = rho_si(rho, x[ix], z[iz], x, z, 2);
			//rho_den1 = std::pow((1.0 - rho_s1_ixiz[ix*nzstep+iz]),2.0);
			rho_den1 = (1.0 - rho_s1_ixiz[ix*nzstep+iz]);
			rho_den1 = rho_den1 * rho_den1;
			//rho_den2 = std::pow((rho_den1 - 4.0*rho_s0_ixiz[ix*nzstep+iz]*rho_s2_ixiz[ix*nzstep+iz]),0.5);
			//rho_den2 = std::sqrt(rho_den1 - 4.0*rho_s0_ixiz[ix*nzstep+iz]*rho_s2_ixiz[ix*nzstep+iz]);
			rho_den2 = rho_den1 - 4.0*rho_s0_ixiz[ix*nzstep+iz]*rho_s2_ixiz[ix*nzstep+iz];
			// to avoide nan
			if ( rho_den2 > 0 ) {
				rho_den2 = std::sqrt(rho_den2);
			} else {
				rho_den2 = 0.0;
			}
			rho_s_ixiz[ix*nzstep+iz] = 2.0*rho_s0_ixiz[ix*nzstep+iz]/(1.0 - rho_s1_ixiz[ix*nzstep+iz]+rho_den2);
			//std::cout << iz << ", " << rho_ixiz[ix*nzstep+iz] << ", " << rho_s_ixiz[ix*nzstep+iz] << ", " << rho_s0_ixiz[ix*nzstep+iz] << ", " << rho_s1_ixiz[ix*nzstep+iz] << ", " << rho_s2_ixiz[ix*nzstep+iz] << std::endl;
			//std::cout << rho_den1 << ", " << rho_den2 << std::endl;
		}
		//std::cout << ix << ", " << rho_ixiz[ix*nzstep+iz] << ", " << rho_s_ixiz[ix*nzstep+iz] << ", " << rho_s0_ixiz[ix*nzstep+iz] << ", " << rho_s1_ixiz[ix*nzstep+iz] << ", " << rho_s2_ixiz[ix*nzstep+iz] << std::endl;
	}
	return 0;
}

// Steele 10-4-3 potential for NLDFT
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

double calc_alpha(double *x, double *z){
	int ix; // x axis for rho
	int iz; // z axis for rho
	int jz; // z axis
	int jx; // x axis
	int t;  // theta
	double ra;  // distance
	double rajz; // z axis
	//
	double phi_att_ff_int_t[ntmesh]; // theta
	double phi_att_ff_int_jx[nxstep]; // x axis
	double phi_att_ff_int_jz[nzstep]; // z axis
	//
	double xt,yt;
	double drad = M_PI/nrmesh; // radian
	double spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	//
	double alpha_int_ix[nxstep];
	double alpha_int_iz[nzstep];
	double alpha_other_method;
	//
	for (ix=0; ix<nxstep; ix++) {
		for (iz=0; iz<nzstep; iz++) {
			//
			// j, z
			for (jz=0; jz<nzstep; jz++) {
				rajz = (z[jz]-z[iz]);
				// k, x
				for (jx=0; jx<nxstep; jx++) {
					for (t=0; t<ntmesh; t++) {
						xt = x[jx]*std::cos(drad*double(t));
						yt = x[jx]*std::sin(drad*double(t));
						ra = (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz*rajz; // x, y, z
						ra = std::sqrt(ra);
						phi_att_ff_int_t[t]  = -phi_att_ff(ra);
					}
					//integral_simpson(double *f, int n, double dx)
					phi_att_ff_int_jx[jx] = 2.0*x[jx]*integral_simpson(phi_att_ff_int_t, ntmesh-1, drad);
				}
				//integral_simpson(double *f, int n, double dx)
				//phi_att_ff_int_jz[jz] = integral_simpson(phi_att_ff_int_jx, nxstep-1, dx) + phi_att_ff_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
				phi_att_ff_int_jz[jz] = integral_simpson(phi_att_ff_int_jx, nxstep-1, dx) + phi_att_ff_int_jx[0]*spr2;
			}
			//integral_simpson(double *f, int n, double dx)
			alpha_int_iz[iz] = integral_simpson(phi_att_ff_int_jz, nzstep-1, dz);
		}
		//integral_simpson(double *f, int n, double dx)
		alpha_int_ix[ix] = integral_simpson(alpha_int_iz, nzstep-1, dz);
	}
	//integral_simpson(double *f, int n, double dx)
	alpha_other_method = integral_simpson(alpha_int_ix, nxstep-1, dx) / ((H-sigma_ss)*D);
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "average alpha of other method = " << alpha_other_method << " " << std::endl;
	return alpha_other_method;
}

double phi_att_ff_int(double *x, double *z, double *phi_att_ff_int_ixizjxjz){
	int ix; // x axis for rho
	int iz; // z axis for rho
	int jz;  // z axis
	int jx;  // x axis
	int t;  // theta
	double ra;  // distance
	double rajz; // z axis
	//
	double phi_att_ff_int_t[ntmesh];  // theta
	double phi_att_ff_int_jx[nxstep]; // x axis
	double phi_att_ff_int_jz[nzstep]; // z axis
	//
	double xt,yt;
	double drad = M_PI/nrmesh; // radian
	double spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	//
	for (ix=0; ix<nxstep; ix++) {
		for (iz=0; iz<nzstep; iz++) {
			//
			// j, z
			for (jz=0; jz<nzstep; jz++) {
				rajz = (z[jz]-z[iz]);
				// k, x
				for (jx=0; jx<nxstep; jx++) {
					for (t=0; t<ntmesh; t++) {
						xt = x[jx]*std::cos(drad*double(t));
						yt = x[jx]*std::sin(drad*double(t));
						ra = (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz*rajz; // x, y, z
						ra = std::sqrt(ra);
						phi_att_ff_int_t[t]  = phi_att_ff(ra);
					}
					//phi_att_ff_int_ixizjxjz[nxstep][nzstep][nxstep][nzstep]
					phi_att_ff_int_ixizjxjz[ix*nzstep*nxstep*nzstep+iz*nxstep*nzstep+jx*nzstep+jz] = 2.0*x[jx]*integral_simpson(phi_att_ff_int_t, ntmesh-1, drad);
				}
				//phi_att_ff_int_jz[jz] = integral_simpson(phi_att_ff_int_jx, nxstep-1, dx);
			}
			//phi_att_ff_int_ixiz[ix*nzstep+iz] = integral_simpson(phi_att_ff_int_jz, nzstep-1, dz);
		}
	}
	return 0;
}

// solid-fluid
double phi_att_sf_int(double *x, double *z, double *rhos_phi_sf_int_ixiz){
	int ix; // x axis for rho
	int iz; // z axis for rho
	int jz;  // wall area
	int jx;  // radius on x-y plane
	int t;  // theta
	double ra;  // distance
	double rajz; // z axis
	//
	int sfzmesh = 20000; // number of step in wall area
	//double h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	double dsf = (h0+2.0*delta)/(sfzmesh-1);
	double wwidth = (h0+2.0*delta);
	//
	double phi_sf_int_jz[sfzmesh];
	double phi_sf_int_jx[nrmesh];
	double phi_sf_int_t[ntmesh];
	//
	double xt,yt;
	double drad = M_PI/nrmesh; // radian
	double spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	//
	for (ix=0; ix<nxstep; ix++) {
		for (iz=0; iz<nzstep; iz++) {
			// under side, z <= 0
			for (jz=0; jz<sfzmesh; jz++) {
				rajz = ((double(jz)*dsf-wwidth)-z[iz]);
				for (jx=0; jx<nxstep; jx++) {
					for (t=0; t<ntmesh; t++) {
						xt = x[jx]*std::cos(drad*double(t));
						yt = x[jx]*std::sin(drad*double(t));
						ra = (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz*rajz; // x, y, z
						ra = std::sqrt(ra);
						phi_sf_int_t[t]  = phi_att_sf(ra);
					}
					phi_sf_int_jx[jx] = 2.0*x[jx]*rho_ssq(double(jz)*dsf)*integral_simpson(phi_sf_int_t, ntmesh-1, drad);
				}
				//phi_sf_int_jz[jz] = integral_simpson(phi_sf_int_jx, nxstep-1, dx) + phi_sf_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
				phi_sf_int_jz[jz] = integral_simpson(phi_sf_int_jx, nxstep-1, dx) + phi_sf_int_jx[0]*spr2;
			}
			// top side, z >= H, z go to a positive value.
			for (jz=0; jz<sfzmesh; jz++) {
				rajz = ((H+double(jz)*dsf)-z[iz]);
				for (jx=1; jx<nxstep; jx++) {
					for (t=0; t<ntmesh; t++) {
						xt = x[jx]*std::cos(drad*double(t));
						yt = x[jx]*std::sin(drad*double(t));
						ra = (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz*rajz; // x, y, z
						ra = std::sqrt(ra);
						phi_sf_int_t[t]  = phi_att_sf(ra);
					}
					phi_sf_int_jx[jx] = 2.0*x[jx]*rho_ssq(double(jz)*dsf)*integral_simpson(phi_sf_int_t, ntmesh-1, drad);
				}
				//phi_sf_int_jz[jz] = phi_sf_int_jz[jz] + integral_simpson(phi_sf_int_jx, nxstep-1, dx) + phi_sf_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
				phi_sf_int_jz[jz] = phi_sf_int_jz[jz] + integral_simpson(phi_sf_int_jx, nxstep-1, dx) + phi_sf_int_jx[0]*spr2;
			}
			//
			//rhos_phi_sf_int_ixiz[ix*nzstep+iz] = integral_simpson(phi_sf_int_jz, sfzmesh-1, dsf);
			//std::cout << "x[ix]=" << x[ix] << ", z[iz]=" << z[iz] << ", rhos_phi_sf_int=" << rhos_phi_sf_int_ixiz[ix*nzstep+iz] << std::endl;
		}
	}
	return 0;
}

// xi include kb1*T*(std::log(rho_b)) type.
// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
double xi(double *rho, double *x, double *z, int x0, int z0, double rho_b, double *rho_s_ixiz, double *rho_s0_ixiz, double *rho_s1_ixiz, double *rho_s2_ixiz, double *phi_att_ff_int_ixizjxjz, double *rho_dfex_int_ixiz, double *rho_phi_int_ixiz){
	int ix; // x axis for rho
	int iz; // z axis for rho
	int jz;  // wall area
	int jx;  // radius on x-y plane
	int t;  // theta
	double ra;  // distance
	double rajz; // z axis
	//
	double rho_dfex_int_t[ntmesh];
	double rho_dfex_int_jx[nxstep];
	double rho_dfex_int_jz[nzstep];
	double rho_phi_int_jx[nxstep];
	double rho_phi_int_jz[nzstep];
	//
	double xt,yt;
	double drad = M_PI/nrmesh; // radian
	double spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	// j, z
	for (jz=0; jz<nzstep; jz++) {
		rajz = (z[jz]-z[z0]);
		// k, x
		for (jx=0; jx<nxstep; jx++) {
			for (t=0; t<ntmesh; t++) {
				xt = x[jx]*std::cos(drad*double(t));
				yt = x[jx]*std::sin(drad*double(t));
				ra = (xt-x[x0])*(xt-x[x0]) + yt*yt + rajz*rajz; // x, y, z
				ra = std::sqrt(ra);
				rho_dfex_int_t[t] = drhos_per_drho_j(ra, rho_s_ixiz[ix*nzstep+iz], rho_s1_ixiz[ix*nzstep+iz], rho_s2_ixiz[ix*nzstep+iz]);
			}
			rho_dfex_int_jx[jx] = 2.0*x[jx]*rho[jx*nzstep+jz]*dfex_per_drhos(rho_s_ixiz[ix*nzstep+iz])*integral_simpson(rho_dfex_int_t, ntmesh-1, drad);
			rho_phi_int_jx[jx]  = rho[jx*nzstep+jz]*phi_att_ff_int_ixizjxjz[ix*nzstep*nxstep*nzstep+iz*nxstep*nzstep+jx*nzstep+jz];
		}
		//rho_dfex_int_jz[jz] = integral_simpson(rho_dfex_int_jx, nxstep-1, dx) + rho_dfex_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
		//rho_phi_int_jz[jz]  = integral_simpson(rho_phi_int_jx, nxstep-1, dx) + rho_phi_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
		rho_dfex_int_jz[jz] = integral_simpson(rho_dfex_int_jx, nxstep-1, dx) + rho_dfex_int_jx[0]*spr2;
		rho_phi_int_jz[jz]  = integral_simpson(rho_phi_int_jx, nxstep-1, dx) + rho_phi_int_jx[0]*spr2;
	}
	//integral_simpson(double *f, int n, double dx)
	rho_dfex_int_ixiz[x0*nzstep+z0] = integral_simpson(rho_dfex_int_jz, nzstep-1, dz);
	rho_phi_int_ixiz[x0*nzstep+z0]  = integral_simpson(rho_phi_int_jz, nzstep-1, dz);
	//
	double xi_out;
	xi_out = ( - rho_b*alpha - rho_dfex_int_ixiz[x0*nzstep+z0] - f_ex(rho_s_ixiz[x0*nzstep+z0]) ) + ( mu_ex(rho_b) - rho_phi_int_ixiz[x0*nzstep+z0] ) + ( kb1*T*std::log(rho_b) );
	// debug
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
double omega(double *rho, double *x, double *z, double *rho_dfex_int_ixiz, double *rho_phi_int_ixiz){
	int ix;
	int iz;
	int omega_nstep = (nzstep-2)/2;
	// x
	double rho_int_ix[nxstep];
	double rho_x_rho_dfex_int_ix[nxstep];
	double rho_x_rho_phi_int_ix[nxstep];
	// z
	double rho_int_iz[omega_nstep+1];
	double rho_x_rho_dfex_int_iz[omega_nstep+1];
	double rho_x_rho_phi_int_iz[omega_nstep+1];
	//
	double omega1, omega2, omega3;
	double omega_out;
	//
	double tpi = 2.0*M_PI;
	double spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	//
	for (iz=0; iz<=omega_nstep; iz++){
		for (ix=0; ix<nxstep; ix++){
			rho_int_ix[ix] = tpi * x[ix] * rho[ix*nzstep+iz];
			rho_x_rho_dfex_int_ix[ix] = tpi * x[ix] * rho[ix*nzstep+iz] * rho_dfex_int_ixiz[ix*nzstep+iz];
			rho_x_rho_phi_int_ix[ix]  = tpi * x[ix] * rho[ix*nzstep+iz] * rho_phi_int_ixiz[ix*nzstep+iz];
		}
		rho_int_iz[iz] = integral_simpson(rho_int_ix, nxstep, dx) + rho_int_ix[0]*spr2;
		rho_x_rho_dfex_int_iz[iz] = integral_simpson(rho_x_rho_dfex_int_ix, nxstep, dx) + rho_x_rho_dfex_int_ix[0]*spr2;
		rho_x_rho_phi_int_iz[iz] = integral_simpson(rho_x_rho_phi_int_ix, nxstep, dx) + rho_x_rho_phi_int_ix[0]*spr2;
	}
	omega1 = -(kb1*T) * integral_simpson(rho_int_iz, omega_nstep, dz);
	omega2 = -integral_simpson(rho_x_rho_dfex_int_iz, omega_nstep, dz);
	omega3 = -0.5 * integral_simpson(rho_x_rho_phi_int_iz, omega_nstep, dz);
	omega_out = (omega1 + omega2 + omega3) * 2.0 / epsilon_ff;
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
	read_parameters();
	double z[nzstep];
	double x[nxstep]; // x = y
	//double rho[nxstep][nzstep]; // [(nzstep+1)*nxstep], a[x][z]= a[x*n+z] for a[][n]
	double *rho     = (double *)malloc(sizeof(double)*((nzstep+1)*nxstep));
	double *rho_new = (double *)malloc(sizeof(double)*((nzstep+1)*nxstep));
	//
	for (i=0; i<nxstep; i++){
		x[i] = dx/2.0 + dx*double(i); // dx = (D/2.0)/double(nxstep);
	}
	for (i=0; i<nzstep; i++){
		z[i] = sigma_ss/2.0 + dz*double(i); // dz = (H-sigma_ss)/double(nzstep+1);
	}
	
	// show alpha
	//calc_alpha(x,z);
	// alpha = calc_alpha(x,z);
	
	// set rho_b0
	if ( rho_b0 != 0.0 ){
		std::cout << "rho_b0 = " << rho_b0 << std::endl;
	} else {
		rho_b0 = Maxwell_construction();
	}
	
	//std::cout << rho_b0 << std::endl;
	// initialization
	int ix, iz;
	for (ix=0; ix<nxstep; ix++){
		for (iz=0; iz<nzstep; iz++){
			rho[iz*nxstep+ix] = rho_b0/(nzstep*dz)/(nxstep*dx);
			rho_new[ix*nzstep+iz] = 0.0;
		}
	}
	// P/P0, V[molecules/nm^3], Omega/epsilon_ff[nm^-2]
	std::ofstream ofsppov_vs("./PP0_vs_Vgamma_data_vs.txt");
	ofsppov_vs << "# w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_vs << "# P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	std::cout << "P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
	//double rho[nxstep][nzstep]; // [(nzstep+1)*nxstep], a[x][z]= a[x*n+z] for a[][n]
	double *rho_s_ixiz  = (double *)malloc(sizeof(double)*((nzstep+1)*nxstep));
	double *rho_s0_ixiz = (double *)malloc(sizeof(double)*((nzstep+1)*nxstep));
	double *rho_s1_ixiz = (double *)malloc(sizeof(double)*((nzstep+1)*nxstep));
	double *rho_s2_ixiz = (double *)malloc(sizeof(double)*((nzstep+1)*nxstep));
	//double phi_att_ff_int_ij[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//Memo: a[i][j][k]= a[i*n*o+j*n+k] for a[][n][o], a[i][j][k][l]= a[i*n*o*p+j*o*p+k*p+l] for a[][n][o][p]
	double *phi_att_ff_int_ixizjxjz = (double *)malloc(sizeof(double)*(nxstep*nzstep*nxstep*nzstep+nzstep*nxstep*nzstep+nxstep*nzstep+nzstep));
	phi_att_ff_int(x, z, phi_att_ff_int_ixizjxjz); // calculate integral phi_att_ff at r[i]
	//double rho[nxstep][nzstep]; // [(nzstep+1)*nxstep], a[x][z]= a[x*n+z] for a[][n]
	double *rhos_phi_sf_int_ixiz  = (double *)malloc(sizeof(double)*((nzstep+1)*nxstep));
	phi_att_sf_int(x, z, rhos_phi_sf_int_ixiz); // calculate integral phi_att_sf at r[i] -> rhos * phi_att_sf
	//
	//double rho[nxstep][nzstep]; // [(nzstep+1)*nxstep], a[x][z]= a[x*n+z] for a[][n]
	double *rho_dfex_int_ixiz  = (double *)malloc(sizeof(double)*((nzstep+1)*nxstep));
	double *rho_phi_int_ixiz   = (double *)malloc(sizeof(double)*((nzstep+1)*nxstep));
	//
	double diff0;
	double mixing;
	//
	double rho_int_ix[nxstep];
	double rho_int_iz[((nzstep-2)/2)];
	for (k=0; k<100; k++){
		rho_b = rho_b0 * std::exp(-(20.0-2.0*double(k+1.0)/10.0));
		//rho_b = rho_b0 * std::exp(-(20.0-2.0*double(99.0-k+1.0)/10.0));
		//std::cout << "--------------------------------------------------" << std::endl;
		//std::cout << "rho_b = " << rho_b << std::endl;
		for (j=0; j<cycle_max; j++){
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dz*nzstep/2. 
			rho_s(rho, x, z, rho_s_ixiz, rho_s0_ixiz, rho_s1_ixiz, rho_s2_ixiz);
			for (ix=0; ix<nxstep; ix++){
				for (iz=0; iz<=(nzstep-2)/2; iz++){
					rho_new[i] = std::exp(xi(rho, x, z, ix, iz, rho_b, rho_s_ixiz, rho_s0_ixiz, rho_s1_ixiz, rho_s2_ixiz, phi_att_ff_int_ixizjxjz, rho_dfex_int_ixiz, rho_phi_int_ixiz)/(kb1*T)-rhos_phi_sf_int_ixiz[ix*nzstep+iz]/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
					//
					// overflow about std::exp(730)
					// to avoid overflow
					if (rho_new[ix*nzstep+iz] > 1e9){
						rho_new[ix*nzstep+iz] = 1e9;
					}
					// to avoid -inf or int
					if (rho_new[ix*nzstep+iz] < 1e-18 && rho[ix*nzstep+iz] < 1e-18){
						rho_new[ix*nzstep+iz] = 1e-18;
						rho[ix*nzstep+iz] = 1e-18;
					}
				}
				diff = 0.0;
				for (ix=0; ix<nxstep; ix++){
					for (iz=0; iz<=(nzstep-2)/2; iz++){
						diff0 = std::abs((rho_new[ix*nzstep+iz]-rho[ix*nzstep+iz])/rho[ix*nzstep+iz]);
						diff = diff + 2.0*diff0;
						mixing = wmixing + wmixing/(0.5+diff0);
						//std::cout << i << ", " << mixing << std::endl;
						rho[ix*nzstep+iz] = mixing*rho_new[ix*nzstep+iz] + (1.0-mixing)*rho[ix*nzstep+iz];
						rho[ix*nzstep+((nzstep-1)-iz)] = rho[ix*nzstep+iz]; // The rest is filled with mirror symmetry. 
					}
				}
				if ( (diff/(nxstep*nzstep)*100.0) < 5.0 && j >= 100) {
					break;
				}
				//std::cout << "--------------------------------------------------" << std::endl;
				//std::cout << "cycle=" << j << ", diff=" << diff << ", rho[nzstep/2]=" << rho[nzstep/2] << std::endl;
			}
		}
		//for (i=0; i<nzstep; i++){
		//	std::cout << "--------------------------------------------------" << std::endl;
		//	std::cout << "cycle=" << j << ", r[" << i << "]" << r[i] << " , -ext " << - phi_ext(r[i]) << ", rho[" << i << "]=" << rho[i] << std::endl;
		//}
		//
		double spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
		for (iz=0; iz<=(nzstep-2)/2; iz++){
			for (ix=0; ix<nxstep; ix++){
				rho_int_ix[ix] = 2.0*M_PI * x[ix] * rho[ix*nzstep+iz];
			}
			rho_int_iz[iz] = integral_simpson(rho_int_ix, nxstep, dx) + rho_int_ix[0]*spr2;
		}
		v_gamma = 2.0 * integral_simpson(rho_int_iz, ((nzstep-2)/2), dz);
		v_gamma = v_gamma/((H-sigma_ss)*M_PI*(D/2.0)*(D/2.0)) - rho_b;
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10.0; // [mmol/cm3]
		//v_cm3STP_per_g = v_mmol_per_cm3 / 22.414 / (rho_ss*12.0107*(1e7*1e7*1e7)/(6.02214076*1e23)); // [cm3(STP)/g], 2.226 [g/cm3]
		v_cm3STP_per_g = v_mmol_per_cm3 / 22.414 / (rho_ss*12.0107*10.0/6.02214076); // [cm3(STP)/g], 2.226 [g/cm3]
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
		grand_potential = omega(rho, x, z, rho_dfex_int_ixiz, rho_phi_int_ixiz);
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_vs << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_g << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_g << ", " << grand_potential << std::endl;
	}
	return 0;
}
