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
float H; // distace of slit. [nm]
float sigma_ss; // [nm]
int nzstep; // number of step between disks
int nxstep; // number of step
float w_pw; // pore width, [nm]
float dz; //
float dx; //
float D;  // diameter of disk
// ---------- ----------- ------------ ------------
// assume rho is same value in x-y plane.
// cylinder and normalization, because of cut off (rc).
int ntmesh; // theta
// ---------- ----------- ------------ ------------
// iteration of rho
int cycle_max;
float wmixing; // mixing parameter
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
float rc;   // fluid-fluid
// ---------- ----------- ------------ ------------
float rm; // rm = std::pow(2.0,1.0/6.0)*sigma_ff = 1.12246205*sigma_ff, minimum position of LJ
// ---------- ----------- ------------ ------------
// Carbon dioxide/Carbon slit 81.5  [K](epsilon), 0.3430 [nm](sigma)
// Nitrogen/Carbon slit       53.72 [K](epsilon), 0.3508 [nm](sigma)
//float epsilon_sf = 53.72; // [K] 
//float sigma_sf = 0.3508; // [nm]
float epsilon_sf;
float sigma_sf;
// C: 52.84 K, 0.343 nm
// N2: 104.2 K, 0.3632 nm
// H2: 36.7 K, 0.2958 nm
// CH4-N2: 95.2 K, 0.3745 nm
// CO2-N2: 101.6 K, 0.3636 nm
// ---------- ----------- ------------ ------------
// slit pore (graphite)
//float delta = 0.335; // [nm]
//float rho_ss = 114.0; // [nm^-3], [molecules/nm3]?, 0.114 [A^-3]
float delta;
float rho_ss;
float h0;
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

float rho_ssq(float z){
	float rho_ssq_out;
	//float h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	//float rho_ss = 114.0; // [molecules/nm3] (the density of bulk carbon)
	//float delta = 0.13; // [nm] (the roughness parameter) (the half-width of the density ramp)
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
	float num[25];
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
		nzstep = int((H-sigma_ss)/0.01 + 0.5);
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
	// move blow (nxtep)
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
	float DH;
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
	// ntmesh
	if ( ntmesh == 0 ) {
		ntmesh = int(nxstep/1.0);
		std::cout << "autoset ntmesh = " << ntmesh << std::endl;
	}
	// ---------- ----------- ------------ ------------
	
	w_pw = (H-sigma_ss); // pore width [nm]
	dz = (H-sigma_ss)/(nzstep-1);
	dx = (D/2.0)/nxstep; // x-y plane
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

float integral_trapezoidal(float *f, int n, float dq){
	float sum;
	sum = 0.0;
	int i;
	for(i=1; i<n; i++){
		sum += (f[i-1]+f[i])/2.0*dq;
	}
	return sum;
}

float integral_simpson(float *f, int n, float dq){
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
	return (dq/3.0)*sum;
}

// The attractive potentials of fluid-fluid interactions.
float phi_att_ff(float r){
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

// The attractive potentials of solid-fluid interactions.
// This case use normal Lennard-Jones（LJ) potential, because the WCA type is different from result of steele potential.
float phi_att_sf(float r){
	float e;
	// Lennard-Jones（LJ) potential
	e = std::pow((sigma_sf/r),6.0);
	e = 4.0*epsilon_sf*( e*e - e );
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

float rho_si_int_t(float *rho_si_int_t_iixizjxjz, float *x, float *z){
	int i;  // type of wi
	int ix; // x axis for rho, x = y, old x0
	int iz; // z axis for rho, old y0
	int jx; // x axis for rho, x = y
	int jz; // z axis for rho
	int t;  // theta
	float rajz; // z axis
	float ra;
	//
	float tmp_rho_si_int_t[ntmesh];
	//
	float xt,yt;
	float drad = M_PI/(ntmesh-1); // radian
	//float spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	//
	for (i=0; i<3; i++) {
		for (ix=0; ix<nxstep; ix++) {
			for (iz; iz<(nzstep-2)/2; iz++) {
				//
				for (jz=0; jz<nzstep; jz++) {
					rajz = (z[jz]-z[iz]);
					for (jx=0; jx<nxstep; jx++) {
						for (t=0; t<ntmesh; t++) {
							xt = x[jx]*std::cos(drad*float(t));
							yt = x[jx]*std::sin(drad*float(t));
							ra= (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz*rajz; // x, y, z
							ra = std::sqrt(ra);
							tmp_rho_si_int_t[t]  = wi(ra,i);
						}
						//integral_simpson(float *f, int n, float dx)
						rho_si_int_t_iixizjxjz[i*nxstep*nzstep*nxstep*nzstep+ix*nzstep*nxstep*nzstep+iz*nxstep*nzstep+jx*nzstep+jz] = 
							2.0*x[jx]*integral_simpson(tmp_rho_si_int_t, ntmesh-1, drad);
						rho_si_int_t_iixizjxjz[i*nxstep*nzstep*nxstep*nzstep+ix*nzstep*nxstep*nzstep+((nzstep-1)-iz)*nxstep*nzstep+jx*nzstep+jz] = 
							rho_si_int_t_iixizjxjz[i*nxstep*nzstep*nxstep*nzstep+ix*nzstep*nxstep*nzstep+iz*nxstep*nzstep+jx*nzstep+jz];
					}
				}
				//
			}
		}
	}
	return 0;
}

// Tarazona theory
float rho_si(float *rho, float *x, float *z, int ix, int iz, int i, float *rho_si_int_t_iixizjxjz){
	int jx; // x axis for rho, x = y
	int jz; // z axis for rho
	int t;  // theta
	float rajz; // z axis
	float ra;
	//
	float rho_si_int_jz[nzstep];
	float rho_si_int_jx[nxstep];
	//float rho_si_int_t[ntmesh];
	//
	float xt,yt;
	//float drad = M_PI/(ntmesh-1); // radian
	float spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	//
	for (jz=0; jz<nzstep; jz++) {
		//rajz = (z[jz]-z[iz]);
		for (jx=0; jx<nxstep; jx++) {
			//for (t=0; t<ntmesh; t++) {
			//	xt = x[jx]*std::cos(drad*float(t));
			//	yt = x[jx]*std::sin(drad*float(t));
			//	ra= (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz*rajz; // x, y, z
			//	ra = std::sqrt(ra);
			//	rho_si_int_t[t]  = wi(ra,i);
			//}
			//integral_simpson(float *f, int n, float dx)
			//rho_si_int_jx[jx] = 2.0*x[jx]*rho[jx*nzstep+jz] * integral_simpson(rho_si_int_t, ntmesh-1, drad);
			rho_si_int_jx[jx] = rho[jx*nzstep+jz] * 
				rho_si_int_t_iixizjxjz[i*nxstep*nzstep*nxstep*nzstep+ix*nzstep*nxstep*nzstep+iz*nxstep*nzstep+jx*nzstep+jz];
		}
		//integral_simpson(float *f, int n, float dx)
		//rho_si_int_jz[jz] = integral_simpson(rho_si_int_jx, nxstep-1, dx) + rho_si_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
		rho_si_int_jz[jz] = integral_simpson(rho_si_int_jx, nxstep-1, dx) + rho_si_int_jx[0]*spr2;
	}
	float rho_si_out;
	//integral_simpson(float *f, int n, float dx)
	rho_si_out = integral_simpson(rho_si_int_jz, nzstep-1, dz);
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
float rho_s(float *rho, float *x, float *z, float *rho_s_ixiz, float *rho_s0_ixiz, float *rho_s1_ixiz, float *rho_s2_ixiz, float *rho_si_int_t_iixizjxjz){
	int ix; // x axis for rho, x = y
	int iz; // z axis for rho
	//
	float rho_den1;
	float rho_den2;
	for (ix=0; ix<nxstep; ix++) {
		for (iz=0; iz<nzstep; iz++) {
			rho_s0_ixiz[ix*nzstep+iz] = rho_si(rho, x, z, ix, iz, 0, rho_si_int_t_iixizjxjz);
			rho_s1_ixiz[ix*nzstep+iz] = rho_si(rho, x, z, ix, iz, 1, rho_si_int_t_iixizjxjz);
			rho_s2_ixiz[ix*nzstep+iz] = rho_si(rho, x, z, ix, iz, 2, rho_si_int_t_iixizjxjz);
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

float calc_alpha(float *x, float *z){
	int ix; // x axis for rho
	int iz; // z axis for rho
	int jz; // z axis
	int jx; // x axis
	int t;  // theta
	float ra;  // distance
	float rajz; // z axis
	//
	float phi_att_ff_int_t[ntmesh]; // theta
	float phi_att_ff_int_jx[nxstep]; // x axis
	float phi_att_ff_int_jz[nzstep]; // z axis
	//
	float xt,yt;
	float drad = M_PI/(ntmesh-1); // radian
	float spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	//
	float alpha_int_ix[nxstep];
	float alpha_int_iz[nzstep];
	float alpha_other_method;
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
						xt = x[jx]*std::cos(drad*float(t));
						yt = x[jx]*std::sin(drad*float(t));
						ra = (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz*rajz; // x, y, z
						ra = std::sqrt(ra);
						phi_att_ff_int_t[t]  = -phi_att_ff(ra);
					}
					//integral_simpson(float *f, int n, float dx)
					phi_att_ff_int_jx[jx] = 2.0*x[jx]*integral_simpson(phi_att_ff_int_t, ntmesh-1, drad);
				}
				//integral_simpson(float *f, int n, float dx)
				//phi_att_ff_int_jz[jz] = integral_simpson(phi_att_ff_int_jx, nxstep-1, dx) + phi_att_ff_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
				phi_att_ff_int_jz[jz] = integral_simpson(phi_att_ff_int_jx, nxstep-1, dx) + phi_att_ff_int_jx[0]*spr2;
			}
			//integral_simpson(float *f, int n, float dx)
			alpha_int_iz[iz] = integral_simpson(phi_att_ff_int_jz, nzstep-1, dz);
		}
		//integral_simpson(float *f, int n, float dx)
		alpha_int_ix[ix] = integral_simpson(alpha_int_iz, nzstep-1, dz);
	}
	//integral_simpson(float *f, int n, float dx)
	alpha_other_method = integral_simpson(alpha_int_ix, nxstep-1, dx) / ((H-sigma_ss)*D);
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "average alpha of other method = " << alpha_other_method << " " << std::endl;
	return alpha_other_method;
}

float phi_att_ff_int(float *x, float *z, float *phi_att_ff_int_ixizjxjz){
	int ix; // x axis for rho
	int iz; // z axis for rho
	int jz;  // z axis
	int jx;  // x axis
	int t;  // theta
	float ra;  // distance
	float rajz; // z axis
	//
	float phi_att_ff_int_t[ntmesh];  // theta
	//float phi_att_ff_int_jx[nxstep]; // x axis
	//float phi_att_ff_int_jz[nzstep]; // z axis
	//
	float xt,yt;
	float drad = M_PI/(ntmesh-1); // radian
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
						xt = x[jx]*std::cos(drad*float(t));
						yt = x[jx]*std::sin(drad*float(t));
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

// solid-fluid, old version
//float phi_att_sf_int(float *x, float *z, float *rhos_phi_sf_int_ixiz){
//	int ix; // x axis for rho
//	int iz; // z axis for rho
//	int jz; // wall area
//	int jx; // radius on x-y plane
//	int t;  // theta
//	float ra;  // distance
//	float rajz; // z axis
//	//
//	int sfzmesh; // number of step in wall area
//	//float h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
//	float dsf;
//	float wwidth;
//	//
//	delta = 0.0;
//	if ( delta == 0.0 ) {
//		sfzmesh = 50; // number of step in wall area
//		dsf = (h0+2.0*delta)/(sfzmesh-1);
//	} else {
//		sfzmesh = 2; // number of step in wall area
//		h0 = delta*sfzmesh;
//		dsf = delta/(sfzmesh-1);
//	}
//	//
//	float phi_sf_int_jz[sfzmesh];
//	float phi_sf_int_jx[nxstep];
//	float phi_sf_int_t[ntmesh];
//	//
//	float xt,yt;
//	float drad = M_PI/(ntmesh-1); // radian
//	float spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
//	//
//	//std::cout << "--------------------------------------------------" << std::endl;
//	//std::cout << "phi_att_sf_int calculation: start" << std::endl;
//	//
//	for (ix=0; ix<nxstep; ix++) {
//		for (iz=0; iz<nzstep; iz++) {
//			rhos_phi_sf_int_ixiz[ix*nzstep+iz] = 0.0;
//			// under side, z <= 0
//			for (jz=0; jz<sfzmesh; jz++) {
//				rajz = (-float(jz)*dsf-z[iz]);
//				for (jx=0; jx<nxstep; jx++) {
//					for (t=0; t<ntmesh; t++) {
//						xt = x[jx]*std::cos(drad*float(t));
//						yt = x[jx]*std::sin(drad*float(t));
//						ra = (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz*rajz; // x, y, z
//						ra = std::sqrt(ra);
//						phi_sf_int_t[t]  = phi_att_sf(ra);
//					}
//					phi_sf_int_jx[jx] = 2.0*x[jx]*rho_ssq(float(jz)*dsf)*delta*integral_simpson(phi_sf_int_t, ntmesh-1, drad);
//					//phi_sf_int_jx[jx] = 2.0*x[jx]*rho_ssq(float(jz)*dsf)*delta*integral_trapezoidal(phi_sf_int_t, ntmesh-1, drad);
//				}
//				//phi_sf_int_jz[jz] = integral_simpson(phi_sf_int_jx, nxstep-1, dx) + phi_sf_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
//				phi_sf_int_jz[jz] = integral_simpson(phi_sf_int_jx, nxstep-1, dx) + phi_sf_int_jx[0]*spr2;
//				//phi_sf_int_jz[jz] = integral_trapezoidal(phi_sf_int_jx, nxstep-1, dx) + phi_sf_int_jx[0]*spr2;
//			}
//			// top side, z >= H, z go to a positive value.
//			for (jz=0; jz<sfzmesh; jz++) {
//				rajz = ((H+float(jz)*dsf)-z[iz]);
//				for (jx=0; jx<nxstep; jx++) {
//					for (t=0; t<ntmesh; t++) {
//						xt = x[jx]*std::cos(drad*float(t));
//						yt = x[jx]*std::sin(drad*float(t));
//						ra = (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz*rajz; // x, y, z
//						ra = std::sqrt(ra);
//						phi_sf_int_t[t]  = phi_att_sf(ra);
//					}
//					phi_sf_int_jx[jx] = 2.0*x[jx]*rho_ssq(float(jz)*dsf)*delta*integral_simpson(phi_sf_int_t, ntmesh-1, drad);
//					//phi_sf_int_jx[jx] = 2.0*x[jx]*rho_ssq(float(jz)*dsf)*delta*integral_trapezoidal(phi_sf_int_t, ntmesh-1, drad);
//				}
//				//phi_sf_int_jz[jz] = phi_sf_int_jz[jz] + integral_simpson(phi_sf_int_jx, nxstep-1, dx) + phi_sf_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
//				phi_sf_int_jz[jz] = phi_sf_int_jz[jz] + integral_simpson(phi_sf_int_jx, nxstep-1, dx) + phi_sf_int_jx[0]*spr2;
//				//phi_sf_int_jz[jz] = phi_sf_int_jz[jz] + integral_trapezoidal(phi_sf_int_jx, nxstep-1, dx) + phi_sf_int_jx[0]*spr2;
//				rhos_phi_sf_int_ixiz[ix*nzstep+iz] = rhos_phi_sf_int_ixiz[ix*nzstep+iz] + phi_sf_int_jz[jz];
//			}
//			//rhos_phi_sf_int_ixiz[ix*nzstep+iz] = integral_simpson(phi_sf_int_jz, sfzmesh-1, dsf);
//			//rhos_phi_sf_int_ixiz[ix*nzstep+iz] = integral_trapezoidal(phi_sf_int_jz, sfzmesh-1, dsf);
//			std::cout << "x[ix]=" << x[ix] << ", z[iz]=" << z[iz] << ", rhos_phi_sf_int=" << rhos_phi_sf_int_ixiz[ix*nzstep+iz] << std::endl;
//		}
//	}
//	//std::cout << "--------------------------------------------------" << std::endl;
//	//std::cout << "phi_att_sf_int calculation: end" << std::endl;
//	return 0;
//}

// solid-fluid, modified version
float phi_att_sf_int(float *x, float *z, float *rhos_phi_sf_int_ixiz){
	int ix; // x axis for rho
	int iz; // z axis for rho
	int jz; // wall area
	int jx; // radius on x-y plane
	int t;  // theta
	//
	float ra_under;  // distance
	float rajz_under; // z axis
	//
	float ra_top;  // distance
	float rajz_top; // z axis
	//
	int sfzmesh; // number of step in wall area
	//float h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	float dsfz;
	//
	if ( delta == 0.0 ) {
		sfzmesh = 100; // number of step in wall area
		dsfz = (h0+2.0*delta)/(sfzmesh-1);
	} else {
		sfzmesh = 2; // number of step in wall area
		h0 = delta*(sfzmesh-1);
		dsfz = h0/(sfzmesh-1);
	}
	int sfnxstep = int(D/0.02);
	int sfntmesh = 2*int(sfnxstep/7.0);
	//
	float rhos_phi_sf_int_jz[sfzmesh];
	float rhos_phi_sf_int_jx[sfnxstep];
	float phi_sf_int_t[sfntmesh];
	float sfx[sfnxstep];
	float dsfx = (D/2.0)/sfnxstep;
	for (jx=0; jx<sfnxstep; jx++) {
		sfx[jx] = dsfx/2.0 + dsfx*jx;
	}
	//
	float xt,yt;
	float drad = M_PI/(sfntmesh-1); // radian
	float spr2 = M_PI*(dsfx/2.0)*(dsfx/2.0) / (2.0*M_PI*sfx[0]);
	//
	for (ix=0; ix<nxstep; ix++) {
		for (iz=0; iz<=(nzstep-2)/2; iz++){
			//
			rhos_phi_sf_int_ixiz[ix*nzstep+iz] = 0.0;
			for (jz=0; jz<sfzmesh; jz++) {
				rajz_under = (-float(jz)*dsfz-z[iz]);
				rajz_top   = ((H+float(jz)*dsfz)-z[iz]);
				for (jx=0; jx<sfnxstep; jx++) {
					for (t=0; t<sfntmesh; t++) {
						xt = sfx[jx]*std::cos(drad*float(t));
						yt = sfx[jx]*std::sin(drad*float(t));
						// under
						ra_under = (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz_under*rajz_under; // x, y, z
						ra_under = std::sqrt(ra_under);
						// top
						ra_top = (xt-x[ix])*(xt-x[ix]) + yt*yt + rajz_top*rajz_top; // x, y, z
						ra_top = std::sqrt(ra_top);
						//
						phi_sf_int_t[t]  = phi_att_sf(ra_under) + phi_att_sf(ra_top);
					}
					rhos_phi_sf_int_jx[jx] = 2.0*sfx[jx]*rho_ssq(float(jz)*dsfz)*delta*integral_simpson(phi_sf_int_t, sfntmesh-1, drad);
					//rhos_phi_sf_int_jx[jx] = 2.0*sfx[jx]*rho_ssq(float(jz)*dsfz)*delta*integral_trapezoidal(phi_sf_int_t, sfntmesh-1, drad);
				}
				//rhos_phi_sf_int_jz[jz] = integral_simpson(rhos_phi_sf_int_jx, nxstep-1, dsfx) + rhos_phi_sf_int_jx[0]/(2.0*M_PI*sfx[0])*M_PI*(dsfx/2.0)*(dsfx/2.0);
				rhos_phi_sf_int_jz[jz] = integral_simpson(rhos_phi_sf_int_jx, sfnxstep-1, dsfx) + rhos_phi_sf_int_jx[0]*spr2;
				//rhos_phi_sf_int_jz[jz] = integral_trapezoidal(rhos_phi_sf_int_jx, nxstep-1, dx) + rhos_phi_sf_int_jx[0]*spr2;
				rhos_phi_sf_int_ixiz[ix*nzstep+iz] = rhos_phi_sf_int_ixiz[ix*nzstep+iz] + rhos_phi_sf_int_jz[jz];
			}
			//
			//rhos_phi_sf_int_ixiz[ix*nzstep+iz] = integral_simpson(rhos_phi_sf_int_jz, sfzmesh-1, dsfz);
			//rhos_phi_sf_int_ixiz[ix*nzstep+iz] = integral_trapezoidal(rhos_phi_sf_int_jz, sfzmesh-1, dsfz);
			rhos_phi_sf_int_ixiz[ix*nzstep+((nzstep-1)-iz)] = rhos_phi_sf_int_ixiz[ix*nzstep+iz];
		}
	}
	//
	//for (ix=0; ix<nxstep; ix++) {
	//	for (iz=0; iz<nzstep; iz++) {
	//		std::cout << "x[ix]=" << x[ix] << ", z[iz]=" << z[iz] << ", rho_phi_sf="<< rhos_phi_sf_int_ixiz[ix*nzstep+iz] << std::endl;
	//	}
	//}
	return 0;
}

// xi include kb1*T*(std::log(rho_b)) type.
// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
float xi(float *rho, float *x, float *z, int ix0, int iz0, float rho_b, float *rho_s_jxjz, float *rho_s0_jxjz, float *rho_s1_jxjz, float *rho_s2_jxjz, float *phi_att_ff_int_ixizjxjz, float *rho_dfex_int_ixiz, float *rho_phi_ff_int_ixiz, float *rhos_phi_sf_int_ixiz){
	int jz;  // wall area
	int jx;  // radius on x-y plane
	int t;  // theta
	float ra;  // distance
	float rajz; // z axis
	//
	float rho_dfex_int_t[ntmesh];
	float rho_dfex_int_jx[nxstep];
	float rho_dfex_int_jz[nzstep];
	float rho_phi_ff_int_jx[nxstep];
	float rho_phi_ff_int_jz[nzstep];
	//
	float xt,yt;
	float drad = M_PI/(ntmesh-1); // radian
	float spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	// j, z
	for (jz=0; jz<nzstep; jz++) {
		rajz = (z[jz]-z[iz0]);
		// k, x
		for (jx=0; jx<nxstep; jx++) {
			for (t=0; t<ntmesh; t++) {
				xt = x[jx]*std::cos(drad*float(t));
				yt = x[jx]*std::sin(drad*float(t));
				ra = (xt-x[ix0])*(xt-x[ix0]) + yt*yt + rajz*rajz; // x, y, z
				ra = std::sqrt(ra);
				rho_dfex_int_t[t] = drhos_per_drho_j(ra, rho_s_jxjz[jx*nzstep+jz], rho_s1_jxjz[jx*nzstep+jz], rho_s2_jxjz[jx*nzstep+jz]);
			}
			rho_dfex_int_jx[jx] = 2.0*x[jx]*rho[jx*nzstep+jz]*dfex_per_drhos(rho_s_jxjz[jx*nzstep+jz])*integral_simpson(rho_dfex_int_t, ntmesh-1, drad);
			rho_phi_ff_int_jx[jx]  = rho[jx*nzstep+jz]*phi_att_ff_int_ixizjxjz[ix0*nzstep*nxstep*nzstep+iz0*nxstep*nzstep+jx*nzstep+jz];
		}
		//rho_dfex_int_jz[jz] = integral_simpson(rho_dfex_int_jx, nxstep-1, dx) + rho_dfex_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
		//rho_phi_ff_int_jz[jz]  = integral_simpson(rho_phi_ff_int_jx, nxstep-1, dx) + rho_phi_ff_int_jx[0]/(2.0*M_PI*x[0])*M_PI*(dx/2.0)*(dx/2.0);
		rho_dfex_int_jz[jz] = integral_simpson(rho_dfex_int_jx, nxstep-1, dx) + rho_dfex_int_jx[0]*spr2;
		rho_phi_ff_int_jz[jz]  = integral_simpson(rho_phi_ff_int_jx, nxstep-1, dx) + rho_phi_ff_int_jx[0]*spr2;
	}
	//integral_simpson(float *f, int n, float dx)
	rho_dfex_int_ixiz[ix0*nzstep+iz0] = integral_simpson(rho_dfex_int_jz, nzstep-1, dz);
	rho_phi_ff_int_ixiz[ix0*nzstep+iz0]  = integral_simpson(rho_phi_ff_int_jz, nzstep-1, dz);
	//
	float xi_out;
	xi_out = ( - rho_b*alpha - rho_dfex_int_ixiz[ix0*nzstep+iz0] - f_ex(rho_s_jxjz[ix0*nzstep+iz0]) ) + ( mu_ex(rho_b) - rho_phi_ff_int_ixiz[ix0*nzstep+iz0] ) + ( kb1*T*std::log(rho_b) -rhos_phi_sf_int_ixiz[ix0*nzstep+iz0]);
	// debug
	//std::cout << "xi_out, -rho_b*alpha, -rho_dfex_int_ixiz[ix0*nzstep+iz0], -f_ex(rho_s_jxjz[ix0*nzstep+iz0]), mu_ex(rho_b),  -rho_phi_ff_int_ixiz[ix0*nzstep+iz0], kb1*T*std::log(rho_b)" << std::endl;
	//std::cout << xi_out << ", " << -rho_b*alpha << ", " << -rho_dfex_int_ixiz[ix0*nzstep+iz0] << ", " <<  -f_ex(rho_s_jxjz[ix0*nzstep+iz0]) << ", " << mu_ex(rho_b) << ", " <<  -rho_phi_ff_int_ixiz[ix0*nzstep+iz0] << ", " << kb1*T*std::log(rho_b) << ", " << -rhos_phi_sf_int_ixiz[ix0*nzstep+iz0] << std::endl;
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
	int iter_max_drhob0 = 250000;
	int iter_max_dmue = 1500;
	float drhob0 = 0.0001;
	float dmue = 0.01;
	float threshold_diff = 0.3;
	float threshold_find = 0.3;
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
float omega(float *rho, float *x, float *z, float *rho_dfex_int_ixiz, float *rho_phi_ff_int_ixiz){
	int ix;
	int iz;
	int omega_nstep = (nzstep-2)/2;
	// z
	float rho_int_iz[omega_nstep+1];
	float rho_x_rho_dfex_int_iz[omega_nstep+1];
	float rho_x_rho_phi_ff_int_iz[omega_nstep+1];
	// x
	float rho_int_ix[nxstep];
	float rho_x_rho_dfex_int_ix[nxstep];
	float rho_x_rho_phi_ff_int_ix[nxstep];
	//
	float omega1, omega2, omega3;
	float omega_out;
	//
	float tpi = 2.0*M_PI;
	float spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
	//
	for (iz=0; iz<=omega_nstep; iz++){
		for (ix=0; ix<nxstep; ix++){
			rho_int_ix[ix] = tpi * x[ix] * rho[ix*nzstep+iz];
			rho_x_rho_dfex_int_ix[ix] = tpi * x[ix] * rho[ix*nzstep+iz] * rho_dfex_int_ixiz[ix*nzstep+iz];
			rho_x_rho_phi_ff_int_ix[ix]  = tpi * x[ix] * rho[ix*nzstep+iz] * rho_phi_ff_int_ixiz[ix*nzstep+iz];
		}
		rho_int_iz[iz] = integral_simpson(rho_int_ix, nxstep, dx) + rho_int_ix[0]*spr2;
		rho_x_rho_dfex_int_iz[iz] = integral_simpson(rho_x_rho_dfex_int_ix, nxstep, dx) + rho_x_rho_dfex_int_ix[0]*spr2;
		rho_x_rho_phi_ff_int_iz[iz] = integral_simpson(rho_x_rho_phi_ff_int_ix, nxstep, dx) + rho_x_rho_phi_ff_int_ix[0]*spr2;
	}
	omega1 = -(kb1*T) * integral_simpson(rho_int_iz, omega_nstep, dz);
	omega2 = -integral_simpson(rho_x_rho_dfex_int_iz, omega_nstep, dz);
	omega3 = -0.5 * integral_simpson(rho_x_rho_phi_ff_int_iz, omega_nstep, dz);
	omega_out = (omega1 + omega2 + omega3) * 2.0 / epsilon_ff;
	return omega_out;
}

int main(){
	//
	float rho_b;
	// from parameters.txt
	read_parameters();
	//
	float z[nzstep];
	float x[nxstep]; // x = y
	//float rho[nxstep][nzstep]; // [(nzstep+1)*nxstep], a[x][z]= a[x*n+z] for a[][n]
	float *rho     = (float *)malloc(sizeof(float)*((nxstep+1)*nzstep));
	float *rho_new = (float *)malloc(sizeof(float)*((nxstep+1)*nzstep));
	//
	int ix;
	for (ix=0; ix<nxstep; ix++){
		x[ix] = dx/2.0 + dx*float(ix); // dx = (D/2.0)/float(nxstep);
	}
	int iz;
	for (iz=0; iz<nzstep; iz++){
		z[iz] = sigma_ss/2.0 + dz*float(iz); // dz = (H-sigma_ss)/float(nzstep+1);
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
	
	std::cout << "--------------------------------------------------" << std::endl;
	//float rho[nxstep][nzstep]; // [(nxstep+1)*nzstep)], a[x][z]= a[x*n+z] for a[][n]
	float *rho_s_ixiz  = (float *)malloc(sizeof(float)*((nxstep+1)*nzstep)); // rho_s_jxjz in xi()
	float *rho_s0_ixiz = (float *)malloc(sizeof(float)*((nxstep+1)*nzstep));
	float *rho_s1_ixiz = (float *)malloc(sizeof(float)*((nxstep+1)*nzstep)); // rho_s1_jxjz in xi()
	float *rho_s2_ixiz = (float *)malloc(sizeof(float)*((nxstep+1)*nzstep)); // rho_s2_jxjz in xi()
	//
	//float phi_att_ff_int_ij[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//Memo: a[i][j][k]= a[i*n*o+j*n+k] for a[][n][o], a[i][j][k][l]= a[i*n*o*p+j*o*p+k*p+l] for a[][n][o][p]
	float *phi_att_ff_int_ixizjxjz = (float *)malloc(sizeof(float)*(nxstep*nzstep*nxstep*nzstep+nzstep*nxstep*nzstep+nxstep*nzstep+nzstep));
	phi_att_ff_int(x, z, phi_att_ff_int_ixizjxjz); // calculate integral phi_att_ff at r[i]
	std::cout << "phi_att_ff_int calculation was finished" << std::endl;
	//
	//float rho[nxstep][nzstep]; // [(nzstep+1)*nxstep], a[x][z]= a[x*n+z] for a[][n]
	float *rhos_phi_sf_int_ixiz  = (float *)malloc(sizeof(float)*(nxstep*nzstep+nzstep));
	phi_att_sf_int(x, z, rhos_phi_sf_int_ixiz); // calculate integral phi_att_sf at r[i] -> rhos * phi_att_sf
	std::cout << "phi_att_sf_int calculation was finished" << std::endl;
	//
	//std::cout << rho_b0 << std::endl;
	// initialization
	float pre_rho;
	for (ix=0; ix<nxstep; ix++){
		for (iz=0; iz<nzstep; iz++){
			pre_rho = rho_b0*-rhos_phi_sf_int_ixiz[ix*nzstep+iz]/2000.0;
			if ( pre_rho <= 0.0 ) {
				rho[iz*nxstep+ix] = 0.0;
			} else {
				rho[iz*nxstep+ix] = pre_rho;
			}
			rho_new[ix*nzstep+iz] = 0.0;
		}
	}
	//
	//float rho[nxstep][nzstep]; // [(nzstep+1)*nxstep], a[x][z]= a[x*n+z] for a[][n]
	float *rho_dfex_int_ixiz  = (float *)malloc(sizeof(float)*((nxstep+1)*nzstep));
	float *rho_phi_ff_int_ixiz = (float *)malloc(sizeof(float)*((nxstep+1)*nzstep));
	//
	// rho_si_int_t_iixizjxjz[i*nxstep*nzstep*nxstep*nzstep+ix*nzstep*nxstep*nzstep+iz*nxstep*nzstep+jx*nzstep+jz]
	float *rho_si_int_t_iixizjxjz = (float *)malloc(sizeof(float)*(2*nxstep*nzstep*nxstep*nzstep+nxstep*nzstep*nxstep*nzstep+nzstep*nxstep*nzstep+nxstep*nzstep+nzstep));
	rho_si_int_t(rho_si_int_t_iixizjxjz, x, z);
	std::cout << "rho_si_int_t calculation was finished" << std::endl;
	//
	float diff_old2 = 1.0;
	float diff_old1 = 1.0;
	float diff;
	float diff0;
	float mixing;
	float threshold = 0.5/100*nzstep*nxstep;
	float trho = 0.0;
	//
	float rho_int_ix[nxstep];
	float rho_int_iz[((nzstep-2)/2)];
	//
	float v_gamma;
	float press_b, press_b0, pp0;
	float press_b_Pa, press_b0_Pa;
	float v_mmol_per_cm3;
	float v_cm3STP_per_cm3;
	float grand_potential;
	//
	int j,k;
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
	ofsppov_vs << "# P/P0, P[Pa], V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	std::cout << "P/P0, P[Pa], V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega/epsilon_ff[1/nm2]" << std::endl;
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
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dz*nzstep/2. 
			rho_s(rho, x, z, rho_s_ixiz, rho_s0_ixiz, rho_s1_ixiz, rho_s2_ixiz, rho_si_int_t_iixizjxjz);
			for (ix=0; ix<nxstep; ix++){
				for (iz=0; iz<=(nzstep-2)/2; iz++){
					//rho_new[ix*nzstep+iz] = std::exp(xi(rho, x, z, ix, iz, rho_b, rho_s_ixiz, rho_s0_ixiz, rho_s1_ixiz, rho_s2_ixiz, phi_att_ff_int_ixizjxjz, rho_dfex_int_ixiz, rho_phi_ff_int_ixiz, rhos_phi_sf_int_ixiz)/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
					xio = xi(rho, x, z, ix, iz, rho_b, rho_s_ixiz, rho_s0_ixiz, rho_s1_ixiz, rho_s2_ixiz, phi_att_ff_int_ixizjxjz, rho_dfex_int_ixiz, rho_phi_ff_int_ixiz, rhos_phi_sf_int_ixiz)/(kb1*T);
					//std::cout << "ix=" << ix << ", iz=" << iz << ", " << rho_new[ix*nzstep+iz] << ", " << -rhos_phi_sf_int_ixiz[ix*nzstep+iz]/(kb1*T) << std::endl;
					//
					if (-14 < xio && xio < 12){
						rho_new[ix*nzstep+iz] = std::exp(xio); // xi include kb1*T*(std::log(rho_b)) type.
					} else if (xio < -14){
						rho_new[ix*nzstep+iz] = 1e-6;
					} else {
						// overflow about std::exp(730)
				    	// to avoid overflow
						rho_new[ix*nzstep+iz] = rho[ix*nzstep+iz] / 10.0;
					}
				}
			}
			diff_old2 = diff_old1;
			diff_old1 = diff;
			diff = 0.0;
			trho = 0.0;
			for (ix=0; ix<nxstep; ix++){
				for (iz=0; iz<=(nzstep-2)/2; iz++){
					trho = trho + rho[ix*nzstep+iz];
					diff0 = std::abs(rho_new[ix*nzstep+iz] - rho[ix*nzstep+iz]);
					diff = diff + diff0;
					mixing = wmixing + wmixing/(0.5+diff0);
					rho[ix*nzstep+iz] = mixing*rho_new[ix*nzstep+iz] + (1.0-mixing)*rho[ix*nzstep+iz];
					rho[ix*nzstep+((nzstep-1)-iz)] = rho[ix*nzstep+iz]; // The rest is filled with mirror symmetry. 
				}
			}
			//if (diff/nstep < 0.005 && diff_old/nstep < 0.005 && j >= 20) {
			//float threshold = 0.5/100*nstep;
			//if (diff < threshold && diff_old1 < threshold && diff_old2 < threshold) {
			if (diff < threshold && diff_old1 < threshold) {
				break;
			}
			//
			//std::cout << "j=" << j << ", ix=" << int(nxstep/2) << ", rho=" << rho[int(nxstep/2)*nzstep+int(nzstep/2)] << ", mixing=" << mixing << ", diff=" << diff << ", trho=" << trho << std::endl;
			//for (ix=0; ix<nxstep; ix++){
			//	std::cout << "j=" << j << ", ix=" << ix << ", rho=" << rho[ix*nzstep+int(nzstep/2)] << ", mixing=" << mixing << ", tdiff=" << (diff/(nxstep*nzstep/2)*100.0) << std::endl;
			//}
		}
		//
		float spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
		for (iz=0; iz<nzstep; iz++){
			for (ix=0; ix<nxstep; ix++){
				rho_int_ix[ix] = 2.0*M_PI * x[ix] * rho[ix*nzstep+iz];
			}
			rho_int_iz[iz] = integral_simpson(rho_int_ix, nxstep, dx) + rho_int_ix[0]*spr2;
		}
		v_gamma = integral_simpson(rho_int_iz, nzstep-1, dz);
		v_gamma = v_gamma/((H-sigma_ss)*(M_PI*(D/2.0)*(D/2.0))) - rho_b;
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
		// 1 [K/nm3] = 1.38064878e-23/(10^-9)^3 [Nm/m3] = 13806.4878 [Pa]
		press_b_Pa = press_b*13806.4878; // [Pa]
		//press_b0_Pa = press_b0*13806.4878; // [Pa]
		//
		pp0 = press_b/press_b0;
		grand_potential = omega(rho, x, z, rho_dfex_int_ixiz, rho_phi_ff_int_ixiz);
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_vs << pp0 << ", " << press_b_Pa << ", " << v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", " << press_b_Pa << ", " << v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
	}
	// reverse
	std::ofstream ofsppov_ls("./PP0_vs_Vgamma_data_ls.txt");
	ofsppov_ls << "# w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_ls << "# P/P0, P[Pa], V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	//std::cout << "w = (H-sigma_ss) = pore width = " << w_pw << " [nm]" << std::endl;
	//std::cout << "P/P0, V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/g], Omega/epsilon_ff[1/nm2]" << std::endl;
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
			// Since it is mirror-symmetric with respect to the z-axis, this routine calculates up to z/2 = dz*nzstep/2. 
			rho_s(rho, x, z, rho_s_ixiz, rho_s0_ixiz, rho_s1_ixiz, rho_s2_ixiz, rho_si_int_t_iixizjxjz);
			for (ix=0; ix<nxstep; ix++){
				for (iz=0; iz<=(nzstep-2)/2; iz++){
					//rho_new[ix*nzstep+iz] = std::exp(xi(rho, x, z, ix, iz, rho_b, rho_s_ixiz, rho_s0_ixiz, rho_s1_ixiz, rho_s2_ixiz, phi_att_ff_int_ixizjxjz, rho_dfex_int_ixiz, rho_phi_ff_int_ixiz, rhos_phi_sf_int_ixiz)/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
					xio = xi(rho, x, z, ix, iz, rho_b, rho_s_ixiz, rho_s0_ixiz, rho_s1_ixiz, rho_s2_ixiz, phi_att_ff_int_ixizjxjz, rho_dfex_int_ixiz, rho_phi_ff_int_ixiz, rhos_phi_sf_int_ixiz)/(kb1*T);
					//std::cout << "ix=" << ix << ", iz=" << iz << ", " << rho_new[ix*nzstep+iz] << ", " << -rhos_phi_sf_int_ixiz[ix*nzstep+iz]/(kb1*T) << std::endl;
					//
					if (-14 < xio && xio < 12){
						rho_new[ix*nzstep+iz] = std::exp(xio); // xi include kb1*T*(std::log(rho_b)) type.
					} else if (xio < -14){
						rho_new[ix*nzstep+iz] = 1e-6;
					} else {
						// overflow about std::exp(730)
				    	// to avoid overflow
						rho_new[ix*nzstep+iz] = rho[ix*nzstep+iz] / 10.0;
					}
				}
			}
			diff_old2 = diff_old1;
			diff_old1 = diff;
			diff = 0.0;
			trho = 0.0;
			for (ix=0; ix<nxstep; ix++){
				for (iz=0; iz<=(nzstep-2)/2; iz++){
					trho = trho + rho[ix*nzstep+iz];
					diff0 = std::abs(rho_new[ix*nzstep+iz] - rho[ix*nzstep+iz]);
					diff = diff + diff0;
					mixing = wmixing + wmixing/(0.5+diff0);
					rho[ix*nzstep+iz] = mixing*rho_new[ix*nzstep+iz] + (1.0-mixing)*rho[ix*nzstep+iz];
					rho[ix*nzstep+((nzstep-1)-iz)] = rho[ix*nzstep+iz]; // The rest is filled with mirror symmetry. 
				}
			}
			//if (diff/nstep < 0.005 && diff_old/nstep < 0.005 && j >= 20) {
			//float threshold = 0.5/100*nstep;
			//if (diff < threshold && diff_old1 < threshold && diff_old2 < threshold) {
			if (diff < threshold && diff_old1 < threshold) {
				break;
			}
		}
		//
		float spr2 = M_PI*(dx/2.0)*(dx/2.0) / (2.0*M_PI*x[0]);
		for (iz=0; iz<nzstep; iz++){
			for (ix=0; ix<nxstep; ix++){
				rho_int_ix[ix] = 2.0*M_PI * x[ix] * rho[ix*nzstep+iz];
			}
			rho_int_iz[iz] = integral_simpson(rho_int_ix, nxstep, dx) + rho_int_ix[0]*spr2;
		}
		v_gamma = integral_simpson(rho_int_iz, nzstep-1, dz);
		v_gamma = v_gamma/((H-sigma_ss)*(M_PI*(D/2.0)*(D/2.0))) - rho_b;
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
		// 1 [K/nm3] = 1.38064878e-23/(10^-9)^3 [Nm/m3] = 13806.4878 [Pa]
		press_b_Pa = press_b*13806.4878; // [Pa]
		//press_b0_Pa = press_b0*13806.4878; // [Pa]
		//
		pp0 = press_b/press_b0;
		grand_potential = omega(rho, x, z, rho_dfex_int_ixiz, rho_phi_ff_int_ixiz);
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_ls << pp0 << ", " << press_b_Pa << ", " << v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", " << press_b_Pa << ", " << v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
	}
	free(rho_s_ixiz);
	free(rho_s0_ixiz);
	free(rho_s1_ixiz);
	free(rho_s2_ixiz);
	free(phi_att_ff_int_ixizjxjz);
	free(rhos_phi_sf_int_ixiz);
	free(rho_dfex_int_ixiz);
	free(rho_phi_ff_int_ixiz);
	free(rho_si_int_t_iixizjxjz);
	return 0;
}
