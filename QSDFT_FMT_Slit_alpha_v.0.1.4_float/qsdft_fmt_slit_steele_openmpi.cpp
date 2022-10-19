#include <fstream>   // for file in and out
#include <iostream>  // for cout
#include <cmath>     // for log, exp
#include <sstream>   // for read parameters
#include <omp.h>     // OpenMP (c++ nldft.cpp -fopenmp) (set OMP_NUM_THREADS=4)
#include <mpi.h>     // OpenMPI (mpic++ nldft.cpp) (set OMP_NUM_THREADS=1)

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
//float H = 1.00; //distace of slit [nm]
//float sigma_ss = 0.34; // [nm]
//int nstep = 100;
//float w_pw = (H-sigma_ss); // pore width [nm]
//float dr = w_pw/float(nstep);
float H;
float sigma_ss;
//#define nstep=1001;
//constexpr int nstep = 1001;
int nstep;
float w_pw;
float dr;
// ---------- ----------- ------------ ------------
int nsteps; // solid
// ---------- ----------- ------------ ------------
// assume rho is same value in x-y plane.
// cylinder and normalization, because of cut off (rc).
//int nrmesh = 20; //rho_si and xi function
int nrmesh;
//int ndmesh = d_hs*nrmesh/rc
int ndmesh;
//float drc = rc/float(nrmesh-1);
float drc;
//float dd = 2.0*d_hs/float(ndmesh-1);
float dd;
// ---------- ----------- ------------ ------------
// iteration of rho
//int cycle_max = 50;
int cycle_max;
//float wmixing = 0.005;
float wmixing;
// ---------- ----------- ------------ ------------
// fluid-fluid
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
float rc;  // for fluid
float rcsf; // for solid-fluid
// ---------- ----------- ------------ ------------
//float rm = std::pow(2.0,1.0/6.0)*sigma_ff; //minimum position of LJ
//float rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
float rm;   // fluid
float rmsf; // solid-fluid
// ---------- ----------- ------------ ------------
// fluid-solid
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
float kb1 = 1.0;  // use all, except thermal de Broglie wavelength calculation.
// ---------- ----------- ------------ ------------
// thermal de Broglie wavelength
//extern float lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm], Maxwell_construction()
float kb = 1.38e-23; //[J/K] (8.61733262e-5 [eV/K])
// ----------
//float m = 14.0067*2.0/(6.02214076e23)/1000; // N2 = 4.65173e-26 [kg]
//float m = 4.65173e-26; //[kg] (N2) (e.g., Ar = 6.63e-26 [kg])
float m;  // for fluid
float ms; // for solid
//extern float T = 77.347; //[K]
// ----------
float T;
// ----------
float h = 6.63e-34; //[Js] (4.135667696e-15 [eVs])
// ----------
float lam;  // for fluid
float lams; // for solid
// Ref: https://www1.doshisha.ac.jp/~bukka/lecture/statistic/pdftext/std-07.pdf
// ---------- ----------- ------------ ------------
// alpha = integal phi_att_ff * -1.0
//extern float alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
//	( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
float alpha;  //fluid
// ---------- ----------- ------------ ------------
// alpha = integal phi_att_ss * -1.0
float alphas; //solid
// ---------- ----------- ------------ ------------
// rho_b0 is related with P0
float rho_b0;
// ---------- ----------- ------------ ------------
float p0;
// ---------- ----------- ------------ ------------
// The edge position ze of the solid wall
float h0;
float ze;
// ---------- ----------- ------------ ------------
// Ris [nm] is the hard-sphere radius of solid (for QSDFT)
float Ris;
// ---------- ----------- ------------ ------------
// graphite wall (Steele)
float epsilon_sfs;
float sigma_sfs;
float deltas;
// ---------- ----------- ------------ ------------
float epsilon_s = 52.84; //[K] (solid)
float sigma_s = 0.343; //[nm] (solid)
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

//Barker-Henderson (BH) perturbation theory
//float d_bh_calc(float epsilon, float sigma){
//	//float epsilon = 94.45;
//	//float sigma = 0.3575;
//	//Lstoskie et al.,
//	float d_bh_out;
//	float Ts = kb1*T/epsilon;
//	d_bh_out = (1.0+0.2977*Ts)/(1.0+0.331637*Ts+0.00104771*Ts*Ts)*sigma;
//	std::cout << "--------------------------------------------------" << std::endl;
//	std::cout << "d = d_hs = " << d_bh_out << " [nm] at " << T << " [K] from Barker-Henderson (BH) perturbation theory" << std::endl;
//	return d_bh_out;
//}

float rho_ssq(float z){
	float rho_ssq_out;
	//float h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	//float rho_ss = 114.0; // [molecules/nm3] (the density of bulk carbon)
	//float delta = 0.13; // [nm] (the roughness parameter) (the half-width of the density ramp)
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
float calc_ze(int ze_nstep){
	int i;
	float ze_out;
	//float h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	//float rho_ss = 114.0; // [molecules/nm3] (the density of bulk carbon)
	//float delta = 0.13; // [nm] (the roughness parameter) (the half-width of the density ramp)
	float rho_ssi[ze_nstep];
	float dss = (2.0*delta)/(ze_nstep-1);
	for (i=0; i<ze_nstep; i++){
		rho_ssi[i] = rho_ssq(h0+float(i)*dss);
	}
	//integral_trapezoidal(float *f, int n, float dx)
	//ze_out = integral_trapezoidal(rho_ssi, ze_nstep-1, dss)/rho_ss + h0;
	//integral_simpson(float *f, int n, float dx)
	ze_out = integral_simpson(rho_ssi, ze_nstep-1, dss)/rho_ss + h0;
	return ze_out;
}

void read_parameters(void){
	std::ifstream ifs("parameters.txt");
	std::string str;
	float num[25];
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
	if ( nstep == 0 ) {
		nstep = int((H-sigma_ss)/0.005 + 0.5) + 20;
		if ( nstep%2 == 1 ){
			nstep = nstep + 1;
		}
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "autoset nstep = " << nstep << std::endl;
	}
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
	if ( rc == 0.0 ) { 
		rc = 5.0*sigma_ff;
		std::cout << "autoset (cut off) rc = " << rc << " [nm]" << std::endl;
	}
	// move below(sigma_sf)
	// ---------- ----------- ------------ ------------
	nrmesh = int(num[9]);
	if ( nrmesh == 0 ) {
		nrmesh = int(rc/0.02 + 0.5);
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
	m = num[14]; //[g/mol], fluid, H2=2.01568, Ar=39.948, N2=28.0134, CO2=44.01, O2=31.998
	m = m/(6.02214076e23)/1000; //[kg]
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
	p0 = num[17]; // [Pa]
	// ---------- ----------- ------------ ------------
	h0 = num[18]; // [nm]
	// ---------- ----------- ------------ ------------
	Ris = num[19]; // [nm], the hard-sphere radius of solid
	std::cout << "The hard-sphere radius of solid =" << Ris << " [nm]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	// ---------- ----------- ------------ ------------
	ms = num[20]; // [kg], solid
	ms = ms/(6.02214076e23)/1000; //C=12.0107
	// ---------- ----------- ------------ ------------
	epsilon_sfs = num[21]; //[K]  graphite wall (Steele)
	sigma_sfs = num[22];   //[nm] graphite wall (Steele)
	deltas = num[23];      //[nm] graphite wall (Steele)
	// ---------- ----------- ------------ ------------
	
	ze = calc_ze(2000);
	std::cout << "The edge position of the solid wall (one side), ze = " << ze << " [nm]" << std::endl;
	float ze_ssf;
	ze_ssf = (ze+sigma_sf);
	std::cout << "ze+sigma_sf = " << ze_ssf << " [nm]" << std::endl;
	std::cout << "(2.0*ze) = " << (2.0*ze) << " [nm]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	
	if ( nstep == 0 ) {
		//nstep = int((H-2.0*h0)/0.02 + 0.5);
		nstep = int((H-(2.0*ze))/0.005 + 0.5);
		if ( nstep%2 == 1 ){
			nstep = nstep + 1;
		}
		std::cout << "autoset nstep = " << nstep << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
	}
	
	// ---------- ----------- ------------ ------------
	
	w_pw = (H-(2.0*ze)); // pore width [nm]
	dr = (H-(2.0*ze))/float(nstep-1);
	rm = 1.12246205*sigma_ff; // 2^(1/6)=1.12246205
	rmsf = 1.12246205*sigma_sf; // 2^(1/6)=1.12246205
	
	// ---------- ----------- ------------ ------------
	
	ndmesh = nrmesh;
	drc = rc/float(nrmesh-1); // xi()
	dd = drc;
	
	// ---------- ----------- ------------ ------------
	
	// thermal de Broglie wavelength of fluid
	//lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm], Maxwell_construction()
	lam = h/std::sqrt(2.0*M_PI*m*kb*T)*1e9; //[nm], Maxwell_construction()
	
	// thermal de Broglie wavelength of solid
	//float ms = 12.0107/(6.02214076e23)/1000; // 1.99442366e-26 [kg]
	lams = h/std::pow((2.0*M_PI*ms*kb*T),0.5)*1e9;
	
	// fluid-fluid
	// alpha = integal phi_att_ff * -1.0
	alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
		( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );
	// rm = rm when the potential is split according to the WCA schem and rm = simga_ff when the LJ potential is split according to the BH decomposition.
	
	// solid-solid
	// alphas = integal phi_att_ss * -1.0
	alphas = (32.0/9.0)*M_PI*epsilon_s*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_s*std::pow(sigma_s,3.0)*
		( 3.0*std::pow((sigma_s/rc),3.0) - std::pow((sigma_s/rc),9.0) );
	// rm = rm when the potential is split according to the WCA schem and rm = simga_ss when the LJ potential is split according to the BH decomposition.
	
	std::cout << "thermal de Broglie wavelength of fluid = " << lam << " [nm]" << std::endl;
	std::cout << "thermal de Broglie wavelength of solid = " << lams << " [nm]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "integal phi_att_ff * -1.0 = alpha = " << alpha << std::endl;
	std::cout << "integal phi_att_ss * -1.0 = alphas = " << alphas << std::endl;
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
	//}else {
	//	e = 0.0;
	}
	//std::cout << e << std::endl;
	return e;
}

// The attractive potentials of solid-fluid interactions.
float phi_att_sf(float r){
	float e;
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

// The attractive potentials of solid-solid interactions.
float phi_att_ss(float r){
	float e;
	// WCA (Weeks-Chandler-Anderson) type
	if (r < rmsf){
		e = - epsilon_s;
	}else if (rmsf <= r && r <= rcsf){
		// Lennard-Jones（LJ) potential
		//e = 4.0*epsilon_s*( std::pow((sigma_s/r),12.0) - std::pow((sigma_s/r),6.0) );
		e = std::pow((sigma_s/r),6.0);
		e = 4.0*epsilon_s*( e*e - e );
	//}else {
	//	e = 0.0;
	}
	//std::cout << e << std::endl;
	return e;
}

// Steele 10-4-3 potential
float phi_sf(float z){
	float phi_sf_out;
	//float epsilon_sfs = 53.72; // [K]
	//float sigma_sfs = 0.3508;  // [nm]
	float sigma_sfs2 = sigma_sfs*sigma_sfs;
	float sfspz = (sigma_sfs/z);
	float sfspz2 = sfspz*sfspz;
	//float deltas = 0.335; // For NLDFT, graphite: 0.6708/2.0=0.3354 [nm]
	float dez = (0.61*deltas+z);
	//phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sfs*std::pow(sigma_sfs,2.0)*deltas*
	//			( (2.0/5.0)*std::pow((sigma_sfs/z),10.0)-std::pow((sigma_sfs/z),4.0)-std::pow(sigma_sfs,4.0)/
	//			(3.0*deltas*std::pow((0.61*deltas+z),3.0)) );
	phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sfs*(sigma_sfs2)*deltas*
				( (2.0/5.0)*std::pow(sfspz2,5.0)-(sfspz2*sfspz2)-(sigma_sfs2*sigma_sfs2)/
				(3.0*deltas*(dez*dez*dez)) );
	return phi_sf_out;
}

// e.g., wall potential (Carbon slit)
float phi_ext(float z){
	float phi_ext_out;
	phi_ext_out = phi_sf(z) + phi_sf(H-z);
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

// The excess hard sphere chemical potential (mus_ex) in the bulk solid.
// mus_ex is calculated by the PY equation.
float mus_ex(float rho_sb){
	float y, mus_ex_out;
	float d_hss = sigma_ss;
	//eta = M_PI*rho_sb*std::pow(d_hss,3.0)/6.0;
	//eta = M_PI*rho_sb*(d_hss*d_hss*d_hss)/6.0;
	y = M_PI*rho_sb*(d_hss*d_hss*d_hss)/6.0;
	//float den1e = (1.0-eta);
	float den1y = (1.0-y);
	//mus_ex_out = kb1*T*(-std::log(1-eta) + eta*(14.0 - 13.0*eta + 5.0*eta*eta)/(2.0*std::pow((1.0-eta),3.0)));
	//mus_ex_out = kb1*T*(-std::log(den1e) + eta*(14.0 - 13.0*eta + 5.0*eta*eta)/(2.0*(den1e*den1e*den1e)));
	mus_ex_out = kb1*T*(-std::log(den1y) + y*(14.0 - 13.0*y + 5.0*y*y)/(2.0*den1y*den1y*den1y));
	return mus_ex_out;
}

float mu_b(float rho_b){
	float mu_id, mu_hs, mu_b_out;
	//mu_id = kb1*T*std::log(std::pow(lam,3.0)*rho_b);
	mu_id = kb1*T*std::log((lam*lam*lam)*rho_b);
	mu_hs = mu_id + mu_ex(rho_b);
	mu_b_out = mu_hs - rho_b*alpha;
	return mu_b_out;
}

float ni_wall(float *r, float *n0_wall_i, float *n1_wall_i, float *n2_wall_i, float *n3_wall_i, float *nv1_wall_i, float *nv2_wall_i) {
	//
	int i;
	int w;
	float rai;
	float xs, xs2;
	float n0, n1, n2, n3, nv1, nv2;
	//
	int nwstep = 500;
	float dw = (ze)/nwstep;
	//
	float n0_wall_w[nwstep];
	float n1_wall_w[nwstep];
	float n2_wall_w[nwstep];
	float n3_wall_w[nwstep];
	float nv1_wall_w[nwstep];
	float nv2_wall_w[nwstep];
	//
	// left
	for (i=0; i<nstep; i++) {
		for (w=0; w<nwstep; w++) {
			rai = (dw*float(w)-r[i]);
			//
			xs2 = (Ris*Ris-rai*rai);
			if ( xs2 >= 0.0 ){
				xs = std::sqrt(xs2);
			} else{
				xs = 0.0;
			}
			//
			n0_wall_w[w] = rho_ssq(dw*float(w))/(2.0*Ris*Ris)*xs;
			n1_wall_w[w] = rho_ssq(dw*float(w))/(2.0*Ris)*xs;
			n2_wall_w[w] = rho_ssq(dw*float(w))*(2.0*M_PI*xs);
			n3_wall_w[w] = rho_ssq(dw*float(w))*(M_PI*xs*xs);
			nv1_wall_w[w] = rho_ssq(dw*float(w))/(2.0*Ris)*(rai/Ris)*xs;
			nv2_wall_w[w] = rho_ssq(dw*float(w))*(rai/Ris)*(2.0*M_PI*xs);
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
			rai = ((H-dw*float(w))-r[i]);
			//
			xs2 = (Ris*Ris-rai*rai);
			if ( xs2 >= 0.0 ){
				xs = std::sqrt(xs2);
			} else{
				xs = 0.0;
			}
			//
			n0_wall_w[w] = rho_ssq(H-dw*float(w))/(2.0*Ris*Ris)*xs;
			n1_wall_w[w] = rho_ssq(H-dw*float(w))/(2.0*Ris)*xs;
			n2_wall_w[w] = rho_ssq(H-dw*float(w))*(2.0*M_PI*xs);
			n3_wall_w[w] = rho_ssq(H-dw*float(w))*(M_PI*xs*xs);
			nv1_wall_w[w] = rho_ssq(H-dw*float(w))/(2.0*Ris)*(rai/Ris)*xs;
			nv2_wall_w[w] = rho_ssq(H-dw*float(w))*(rai/Ris)*(2.0*M_PI*xs);
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

float ni(float *rho, float *r, int i, float *n0_j, float *n1_j, float *n2_j, float *n3_j, float *nv1_j, float *nv2_j,
		  float *n0, float *n1, float *n2, float *n3, float *nv1, float *nv2,
		  float *n0_wall_i, float *n1_wall_i, float *n2_wall_i, float *n3_wall_i, float *nv1_wall_i, float *nv2_wall_i){
	int j;
	float raj;
	float xf, xs, xf2, xs2;
	float Rif;
	Rif = d_hs/2.0; // [nm], Rif is the hard-sphere radius of fluid
	//float Ris;
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
    //integral_trapezoidal(float *f, int n, float dx)
	//n0[i] = integral_trapezoidal(n0_j, nstep-1, dr) + n0_wall_i[i];
	//n1[i] = integral_trapezoidal(n1_j, nstep-1, dr) + n1_wall_i[i];
	//n2[i] = integral_trapezoidal(n2_j, nstep-1, dr) + n2_wall_i[i];
	//n3[i] = integral_trapezoidal(n3_j, nstep-1, dr) + n3_wall_i[i];
	//nv1[i] = integral_trapezoidal(nv1_j, nstep-1, dr) + nv1_wall_i[i];
	//nv2[i] = integral_trapezoidal(nv2_j, nstep-1, dr) + nv2_wall_i[i];
	//
	//integral_simpson(float *f, int n, float dx)
	n0[i] = integral_simpson(n0_j, nstep-1, dr) + n0_wall_i[i];
	n1[i] = integral_simpson(n1_j, nstep-1, dr) + n1_wall_i[i];
	n2[i] = integral_simpson(n2_j, nstep-1, dr) + n2_wall_i[i];
	n3[i] = integral_simpson(n3_j, nstep-1, dr) + n3_wall_i[i];
	nv1[i] = integral_simpson(nv1_j, nstep-1, dr) + nv1_wall_i[i];
	nv2[i] = integral_simpson(nv2_j, nstep-1, dr) + nv2_wall_i[i];
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
float dfex(float *r, int i, float *n0, float *n1, float *n2, float *n3, float *nv1, float *nv2){
	int j;
	float raj;
	float x, x2;
	float dfex_out;
	float dphi_per_n0, dphi_per_n0_j[nstep];
	float dphi_per_n1, dphi_per_n1_j[nstep];
	float dphi_per_n2, dphi_per_n2_j[nstep];
	float dphi_per_n3, dphi_per_n3_j[nstep];
	float dphi_per_nv1, dphi_per_nv1_j[nstep];
	float dphi_per_nv2, dphi_per_nv2_j[nstep];
	float sxi;
	float sign;
	float Rif;
	Rif = d_hs/2.0; // [nm] Rif is the hard-sphere radius of fluid
	//float Ris;
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
    //integral_trapezoidal(float *f, int n, float dx)
	//dphi_per_n0  = integral_trapezoidal(dphi_per_n0_j, nstep-1, dr);
	//dphi_per_n1  = integral_trapezoidal(dphi_per_n1_j, nstep-1, dr);
	//dphi_per_n2  = integral_trapezoidal(dphi_per_n2_j, nstep-1, dr);
	//dphi_per_n3  = integral_trapezoidal(dphi_per_n3_j, nstep-1, dr);
	//dphi_per_nv1 = integral_trapezoidal(dphi_per_nv1_j, nstep-1, dr);
	//dphi_per_nv2 = integral_trapezoidal(dphi_per_nv2_j, nstep-1, dr);
	//
	//integral_simpson(float *f, int n, float dx)
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

float calc_alpha(float *r){
	int i,j,k;
	float ra;
	float raj;
	float rak;
	//float drc = rc/float(nrmesh-1);
	float tpidrc = 2.0*M_PI*drc;
	float alpha_other_method;
	float alpha_int_j[nstep];
	float alpha_int_k[nrmesh];
	alpha_int_k[0] = 0.0;
	//
	for (i=0; i<=(nstep-2)/2; i++){
		for (j=0; j<nstep; j++) {
			raj = (r[i]-r[j]);
			for (k=1; k<nrmesh; k++) {
				rak = drc*float(k);
				//ra = std::pow((r[i]-r[j]),2.0) + std::pow((float(k)*drc),2.0);
				//ra = (r[i]-r[j])*(r[i]-r[j]) + (float(k)*drc)*(float(k)*drc);
				ra = raj*raj + rak*rak;
				//ra = std::pow(ra,0.5);
				ra = std::sqrt(ra);
				//std::cout << ra << std::endl;
				//alpha_int_k[k]  = -phi_att_ff(ra)*(2.0*M_PI*(float(k)*drc));
				alpha_int_k[k]  = -phi_att_ff(ra)*(tpidrc*float(k));
			}
			//integral_simpson(float *f, int n, float dx)
			alpha_int_j[j]  = integral_simpson(alpha_int_k, nrmesh-1, drc);
		}
		//integral_simpson(float *f, int n, float dx)
		//alpha_other_method  = alpha_other_method + integral_simpson(alpha_int_j, nstep-1, dr)*2.0*dr;
		alpha_other_method  = alpha_other_method + integral_simpson(alpha_int_j, nstep-1, dr);
	}
	alpha_other_method  = alpha_other_method * 2.0 * dr / (H-(2.0*ze));
	//std::cout << "--------------------------------------------------" << std::endl;
	//std::cout << "average alpha of other method = " << alpha_other_method << " in (carbon) slit" << std::endl;
	return alpha_other_method;
}

// fluid-fluid
float phi_att_ff_int(float *r, float *phi_att_ff_int_ij){
	int i,j,k;
	float ra;
	float raj;
	float rak;
	//float drc = rc/float(nrmesh-1);
	//dd = drc;
	float phi_ff_int_k[nrmesh];
	float tpidrc = 2.0*M_PI*drc;
	phi_ff_int_k[0] = 0.0;
	//
	for (i=0; i<nstep; i++) {
		for (j=0; j<nstep; j++) {
			raj = (r[i]-r[j]);
			for (k=1; k<nrmesh; k++) {
				rak = drc*float(k);
				//ra = std::pow((r[i]-r[j]),2.0) + std::pow((float(k)*drc),2.0);
				//ra = (r[i]-r[j])*(r[i]-r[j]) + (float(k)*drc)*(float(k)*drc);
				ra = raj*raj + rak*rak;
				//ra = std::pow(ra,0.5);
				ra = std::sqrt(ra);
				//std::cout << ra << std::endl;
				//rho_phi_ff_int_k[k]  = rho[j]*phi_att_ff(ra)*(2.0*M_PI*(float(k)*drc)); // old ver.1.1.0
				//phi_ff_int_k[k]  = phi_att_ff(ra)*(2.0*M_PI*(float(k)*drc));
				phi_ff_int_k[k]  = phi_att_ff(ra)*(tpidrc*float(k));
			}
		phi_att_ff_int_ij[i*nstep+j] = integral_simpson(phi_ff_int_k, nrmesh-1, drc);
		}
	}
	return 0;
}

// solid-fluid
float phi_att_sf_int(float *r, float *rhos_phi_sf_int_i){
	int i,j,k;
	float ra_left;
	float ra_right;
	float raj_left;
	float raj_right;
	float rak;
	//dd = drc = rc/float(nrmesh-1);
	//
	int sfmesh = 500;
	float dsf = (h0+2.0*delta)/(sfmesh-1);
	float rhos_phi_sf_int_j[sfmesh];
	//
	int sfnrmesh = 500;
	float drcsf = rcsf/(sfnrmesh-1);
	float phi_sf_int_k[sfnrmesh];
	//
	float tpi = 2.0*M_PI;
	phi_sf_int_k[0] = 0.0;
	//
	//for (i=0; i<nstep; i++) {
	for (i=0; i<=(nstep-2)/2; i++){
		//
		for (j=0; j<sfmesh; j++) {
			raj_left  = (float(j)*dsf - r[i]);
			raj_right = ((H-float(j)*dsf) - r[i]);
			for (k=1; k<sfnrmesh; k++) {
				rak = drcsf*float(k);
				//
				ra_left  = rak*rak + raj_left*raj_left;
				ra_left  = std::sqrt(ra_left);
				//
				ra_right = rak*rak + raj_right*raj_right;
				ra_right = std::sqrt(ra_right);
				//
				phi_sf_int_k[k]  = ( phi_att_sf(ra_left) + phi_att_sf(ra_right) ) * (tpi*rak);
			}
			rhos_phi_sf_int_j[j] = rho_ssq(float(j)*dsf)*integral_simpson(phi_sf_int_k, sfnrmesh-1, drcsf);
			//rhos_phi_sf_int_j[j] = rho_ssq(float(j)*dsf)*integral_trapezoidal(phi_sf_int_k, sfnrmesh-1, drcsf);
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

// solid-solid
float phi_att_ss_int(float *r, float *rhos_phi_ss_int_i){
	int i,j,k;
	float ra_left;
	float ra_right;
	float raj_left;
	float raj_right;
	float rak;
	//dd = drc = rc/float(nrmesh-1);
	//
	int sfmesh = 500;
	float dsf = (h0+2.0*delta)/(sfmesh-1);
	float rhos_phi_ss_int_j[sfmesh];
	//
	int sfnrmesh = 500;
	float drcsf = rcsf/(sfnrmesh-1);
	float phi_ss_int_k[sfnrmesh];
	//
	float tpi = 2.0*M_PI;
	phi_ss_int_k[0] = 0.0;
	//
	//for (i=0; i<nstep; i++) {
	for (i=0; i<=(nstep-2)/2; i++){
		//
		for (j=0; j<sfmesh; j++) {
			raj_left  = (float(j)*dsf - r[i]);
			raj_right = ((H-float(j)*dsf) - r[i]);
			for (k=1; k<sfnrmesh; k++) {
				rak = drcsf*float(k);
				//
				ra_left  = rak*rak + raj_left*raj_left;
				ra_left  = std::sqrt(ra_left);
				//
				ra_right = rak*rak + raj_right*raj_right;
				ra_right = std::sqrt(ra_right);
				//
				phi_ss_int_k[k]  = ( phi_att_ss(ra_left) + phi_att_ss(ra_right) ) * (tpi*rak);
			}
			rhos_phi_ss_int_j[j] = rho_ssq(float(j)*dsf)*integral_simpson(phi_ss_int_k, sfnrmesh-1, drcsf);
			//rhos_phi_ss_int_j[j] = rho_ssq(float(j)*dsf)*integral_trapezoidal(phi_ss_int_k, sfnrmesh-1, drcsf);
		}
		//
		rhos_phi_ss_int_i[i] = integral_simpson(rhos_phi_ss_int_j, sfmesh-1, dsf);
		//rhos_phi_ss_int_i[i] = integral_trapezoidal(rhos_phi_ss_int_j, sfmesh-1, dsf);
		rhos_phi_ss_int_i[(nstep-1)-i] = rhos_phi_ss_int_i[i];
	}
	//for (i=0; i<nstep; i++) {
	//	std::cout << i << ", " << r[i] << ", " << rhos_phi_ss_int_i[i] << std::endl;
	//}
	return 0;
}

// xi include kb1*T*(std::log(rho_b)) type.
// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
//float xi(float *rho, float *r, int i, float rho_b, float *rho_sj, float *rho_s0j, float *rho_s1j, float *rho_s2j, float *phi_att_ff_int_ij){
float xi(float *rho, float *r, int i, float rho_b, float *phi_att_ff_int_ij, float *rho_phi_ff_int_i, float *phi_ext_i){
	int j;
	float rho_phi_ff_int_j[nstep];
	for (j=0; j<nstep; j++) {
		rho_phi_ff_int_j[j]  = rho[j]*phi_att_ff_int_ij[i*nstep+j];
	}
	//integral_simpson(float *f, int n, float dx)
	rho_phi_ff_int_i[i]  = integral_simpson(rho_phi_ff_int_j, nstep-1, dr);
	//
	float xi_out;
	xi_out = ( - rho_b*alpha ) + ( mu_ex(rho_b) - rho_phi_ff_int_i[i] ) + ( kb1*T*std::log(rho_b) - phi_ext_i[i] );
	//
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

// For FMT, from Percus Yevick (PY) equation
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
	float dmue = 0.005;
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
	std::cout << "Bulk pressure, P0 = " << press_b0*kb*1e27 << " [Pa] = " << press_b0*kb*1e27/101325.0 << " [atm], (rho_b0 = " << rho_b0_out << ")" <<std::endl;
	std::cout << std::endl;
	std::cout << "gas phase   : rho_b0_gas        = " << rho_b0_gas        << ", rho_b0_gas*d_hs^3        = " << rho_b0_gas*std::pow(d_hs,3.0)        << std::endl;
	std::cout << "metastable  : rho_b0_metastable = " << rho_b0_metastable << ", rho_b0_metastable*d_hs^3 = " << rho_b0_metastable*std::pow(d_hs,3.0) << std::endl;
	std::cout << "liquid phase: rho_b0_liquid     = " << rho_b0_liquid     << ", rho_b0_liquid*d_hs^3     = " << rho_b0_liquid*std::pow(d_hs,3.0)     << std::endl;
	return rho_b0_out;
}

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
float omega(float *rho, float *r, float *fex_i, float *rho_phi_ff_int_i, float *rhos_phi_sf_int_i, float rho_b, float *rhos, float *rhos_phi_ss_int_i){
	float omega_out;
	omega_out = 1.0;
	float omega1_ff, omega1_ss, omega2, omega3_ff, omega3_sf, omega3_ss, omega4_ff, omega4_ss;
	int i;
	float fidf[nstep];
	float rho_x_rho_phi_ff_int[nstep];
	float rho_x_rhos_phi_sf_int[nstep];
	float rho_x_muf[nstep];
	float muf = (kb1*T)*std::log(rho_b*lam*lam*lam) + mu_ex(rho_b) - rho_b*alpha;
	//
	float fids[nsteps];
	float rho_x_rhos_phi_ss_int[nsteps];
	float rho_x_mus[nsteps];
	float rho_sb = rho_b; // dummy
	float mus = (kb1*T)*std::log(rho_sb*lams*lams*lams) + mus_ex(rho_sb) - rho_sb*alphas;
	//float h0 = 2.0*0.34; // [nm] (the thickness of the solid wall)
	//float rho_ss = 114.0; // [molecules/nm3] (the density of bulk carbon)
	//float delta = 0.13; // [nm] (the roughness parameter) (the half-width of the density ramp)
	//
	// fluid
#pragma omp parallel for
	for (i=0; i<nstep; i++){
		fidf[i] = rho[i]*(std::log(rho[i]*lam*lam*lam)-1.0); //fluid
		rho_x_rho_phi_ff_int[i] = rho[i] * rho_phi_ff_int_i[i];
		rho_x_rhos_phi_sf_int[i] = rho[i] * rhos_phi_sf_int_i[i];
		rho_x_muf[i] = rho[i] * - muf;
		//
	}
	// solid
#pragma omp parallel for
	for (i=0; i<nsteps; i++){
		fids[i] = 2.0*rhos[i]*(std::log(rhos[i]*lams*lams*lams)-1.0);
		rho_x_rhos_phi_ss_int[i] = 2.0*(rhos[i] * rhos_phi_ss_int_i[i]);
		rho_x_mus[i] = 2.0*(rhos[i] * - mus);
	}
	//
	omega1_ff = (kb1*T) * integral_simpson(fidf, nstep-1, dr); //fluid
	omega1_ss = (kb1*T) * integral_simpson(fids, nsteps-1, dr); //solid
	//
	omega2 = (kb1*T) * integral_simpson(fex_i, nstep-1, dr); // Fex (fluid + solid)
	//
	omega3_ff = 0.5 * integral_simpson(rho_x_rho_phi_ff_int, nstep-1, dr); //fluid-fluid
	omega3_sf = 0.5 * integral_simpson(rho_x_rhos_phi_sf_int, nstep-1, dr); //solid-fluid
	omega3_ss = 0.5 * integral_simpson(rho_x_rhos_phi_ss_int, nsteps-1, dr); //solid-solid
	//
	omega4_ff = integral_simpson(rho_x_muf, nstep-1, dr); //fluid
	omega4_ss = integral_simpson(rho_x_mus, nsteps-1, dr); //solid
	//
	//Omega[rho_f(r);rho_s(r)] = Fint[rho_f(r);rho_s(r)] + omega4_ff + omega4_ss
	//Fint[rho_f(r);rho_s(r)] = Fid[rho_f(r);rho_s(r)] + Fex[rho_f(r);rho_s(r)]
	//Fid[rho_f(r);rho_s(r)] = Fid[rho_f(r)] + Fid[rho_s(r)]
	//Fid[rho_f(r)] = omega1_ff, Fid[rho_s(r)] = omega1_ss
	//Fex[rho_f(r);rho_s(r)] = omega2 + omega3_ff + omega3_sf + omega3_ss
	//Fex_hs[rho_f(r);rho_s(r)] = omega2
	omega_out = (omega1_ff + omega1_ss + omega2 + omega3_ff + omega3_sf + omega3_ss + omega4_ff + omega4_ss) / epsilon_ff;
	//std::cout << "omega1_ff=" << omega1_ff << ", omega1_ss=" << omega1_ss << std::endl;
	//std::cout << "   omega2=" << omega2 << std::endl;
	//std::cout << "omega3_ff=" << omega3_ff << ", omega3_sf=" << omega3_sf << ", omega3_ss=" << omega3_ss << std::endl;
	//std::cout << "omega4_ff=" << omega4_ff << ", omega4_ss=" << omega4_ss << std::endl;
	return omega_out;
}

//int main(int argc, char **argv){
//MPI::Init(argc,argv);
int main(){
MPI::Init();
	int i,j,k;
	float v_gamma;
	float press_b, press_b0, pp0;
	float press_b_Pa, press_b0_Pa;
	float rho_b;
	float v_mmol_per_cm3;
	float v_cm3STP_per_cm3;
	float grand_potential;
	//
	read_parameters();
	float r[nstep];
	float rho[nstep], rho_new[nstep];
	float rhos[nstep];
	//
#pragma omp parallel for
	for (i=0; i<nstep; i++){
		r[i] = (2.0*ze)/2.0 + dr*float(i);
		//std::cout << i << ", " << r[i] << std::endl;
	}
	
	// show alpha
	//calc_alpha(r);
	// alpha = calc_alpha(r);
	
	// set rho_b0
	float y, a, b, c;
	float flag_P; flag_P = 0.0;
	float rho_b1; rho_b1 = 0.0;
	if ( rho_b0 != 0.0 ){
		std::cout << "rho_b0 = " << rho_b0 << std::endl;
	} else if ( rho_b0 == 0.0 ) {
		rho_b0 = Maxwell_construction();
	} else {
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
		rho[i] = rho_b0/(nstep*dr)*3.91276e-08;
		rho_new[i] = 0.0;
		//
		// solid (one side)
		if ( 0.0 <= r[i] && r[i] < h0 ){
			rhos[i] = rho_ss;
		} else if ( h0 <= r[i] && r[i] < h0+2.0*delta ){
			rhos[i] = 0.75*rho_ss * (1.0 - (r[i] - h0)/(2.0*delta));
			nsteps = i;
		} else {
			rhos[i] = 0.0;
		}
	}
	if (nsteps%2 == 0) {
		nsteps += 1;
		rhos[i] = 0.0;
	}
	
	std::cout << "--------------------------------------------------" << std::endl;
	float *phi_att_ff_int_ij = (float *)malloc(sizeof(float)*((nstep+1)*nstep));
	if (phi_att_ff_int_ij == NULL) {
		printf("Memory cannot be allocated.");
		std::exit(1);
	} else {
		printf("Memory has been allocated. The address is %p\n", phi_att_ff_int_ij);
	}
	phi_att_ff_int(r, phi_att_ff_int_ij); // calculate integral phi_att_ff at r[i]
	//float phi_att_sf_int_i[nstep];
	float rho_phi_ff_int_i[nstep];
	float rhos_phi_sf_int_i[nstep];
	phi_att_sf_int(r, rhos_phi_sf_int_i); // calculate integral phi_att_sf at r[i] -> rhos * phi_att_sf
	//
	float rhos_phi_ss_int_i[nstep];
	phi_att_ss_int(r, rhos_phi_ss_int_i); // calculate integral phi_att_ss at r[i] -> rhos * phi_att_ss
	//
	float n0_wall_i[nstep];
	float n1_wall_i[nstep];
	float n2_wall_i[nstep];
	float n3_wall_i[nstep];
	float nv1_wall_i[nstep];
	float nv2_wall_i[nstep];
	ni_wall(r, n0_wall_i, n1_wall_i, n2_wall_i, n3_wall_i, nv1_wall_i, nv2_wall_i);
	//
	float phi_ext_i[nstep];
	for (i=0; i<nstep; i++){
		phi_ext_i[i] = phi_ext(r[i]);
	}
	std::cout << "phi_ext_i calculation was finished" << std::endl;
	//
	float n0_j[nstep], n0[nstep];
	float n1_j[nstep], n1[nstep];
	float n2_j[nstep], n2[nstep];
	float n3_j[nstep], n3[nstep];
	float nv1_j[nstep], nv1[nstep];
	float nv2_j[nstep], nv2[nstep];
	float c1;
	float fex_i[nstep]; // For grand potential, Omega
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
	ofsppov_vs << "# w = (H-(2.0*ze)) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_vs << "# P[" << Punit << "], V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Relative_Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "w = (H-(2.0*ze)) = pore width = " << w_pw << " [nm]" << std::endl;
	std::cout << "P[" << Punit << "], V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Relative_Omega/epsilon_ff[1/nm2]" << std::endl;
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
			for (i=0; i<nstep; i++){
				ni(rho, r, i, n0_j, n1_j, n2_j, n3_j, nv1_j, nv2_j, n0, n1, n2, n3, nv1, nv2, n0_wall_i, n1_wall_i, n2_wall_i, n3_wall_i, nv1_wall_i, nv2_wall_i);
			}
			for (i=0; i<nstep; i++){
				c1 = dfex(r, i, n0, n1, n2, n3, nv1, nv2);
				//rho_new[i] = std::exp(c1+xi(rho,r,i,rho_b, phi_att_ff_int_ij,rho_phi_ff_int_i)/(kb1*T)-rhos_phi_sf_int_i[i]/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				xio = c1+xi(rho,r,i,rho_b, phi_att_ff_int_ij,rho_phi_ff_int_i, phi_ext_i)/(kb1*T)-rhos_phi_sf_int_i[i]/(kb1*T);
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
			diff_old1 = diff;
			diff = 0.0;
#pragma omp parallel for
			for (i=0; i<=(nstep-2)/2; i++){
				diff0 = std::abs((rho_new[i]-rho[i])/rho[i]);
				diff = diff + 2.0*diff0;
				rho[i] = wmixing*rho_new[i] + (1.0-wmixing)*rho[i];
				rho[(nstep-1)-i] = rho[i]; // The rest is filled with mirror symmetry. 
			}
			//std::cout << "diff=" << diff << std::endl;
			if (diff < threshold && diff_old1 < threshold) {
				break;
			}
		}
		//for (i=0; i<nstep; i++){
		//	std::cout << "i=" << i << ", r=" << r[i] << ", rho_new=" << rho_new[i] << ", rho=" << rho[i] << ", (diff/nstep*100.0)=" << (diff/nstep*100.0) << std::endl;
		//}
		//
		v_gamma = integral_simpson(rho, nstep-1, dr);
		//v_gamma = v_gamma/(H-sigma_ss) - rho_b; // for NLDFT
		v_gamma = v_gamma/(H-(2.0*ze)) - rho_b;
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076 ) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10; // [mmol/cm3]
		//v_cm3STP_per_cm3 = v_mmol_per_cm3 / 22.414 / (rho_ss*12.0107*(1e7*1e7*1e7)/(6.02214076*1e23)); // [cm3(STP)/g], 2.226 [g/cm3]
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
		if(flag_P==0.0){
			pp0 = press_b/press_b0;
		} else if (flag_P<=-10.0){
			pp0 = press_b*kb*1e27/p0;
		} else {
			// kb1=1, kb = 1.38e-23 [J/K], T [K], rho_b [N/nm^3], 1 [atm] = 101325 [Pa]
			pp0 = press_b*kb*1e27;
		}
		// 1 [K/nm3] = 1.38064878e-23/(10^-9)^3 [Nm/m3] = 13806.4878 [Pa]
		// grand ppotential, Omega(ff+sf part)
		for (i=0; i<nstep; i++){
			fex_i[i] = fex(i, n0, n1, n2, n3, nv1, nv2);
		}
		grand_potential = omega(rho, r, fex_i, rho_phi_ff_int_i, rhos_phi_sf_int_i, rho_b, rhos, rhos_phi_ss_int_i);
		//grand_potential = 1.0;		//grand_potential = 1.0;
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_vs << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
	}
	// reverse
	// P/P0, V[molecules/nm^3], Omega/epsilon_ff[nm^-2]
	std::ofstream ofsppov_ls("./PP0_vs_Vgamma_data_ls.txt");
	ofsppov_ls << "# w = (H-(2.0*ze)) = pore width = " << w_pw << " [nm]" << std::endl;
	ofsppov_ls << "# P/P0, P[Pa], V[molecules/nm3], V[mmol/cm3], V[cm3(STP)/cm3], Relative_Omega/epsilon_ff[1/nm2]" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	//
	for (k=181; k>=0; k--){
		//rho_b = rho_b0 * rho_b_k[k];
		if(flag_P<=-100.0){
			rho_b = (rho_b0 - rho_b1) * (rho_b_k[k] - 3.91276e-08) + rho_b1;
		} else {
			rho_b = rho_b0 * rho_b_k[k];
		}
		//
		for (j=0; j<cycle_max; j++){
			for (i=0; i<nstep; i++){
				ni(rho, r, i, n0_j, n1_j, n2_j, n3_j, nv1_j, nv2_j, n0, n1, n2, n3, nv1, nv2, n0_wall_i, n1_wall_i, n2_wall_i, n3_wall_i, nv1_wall_i, nv2_wall_i);
			}
			for (i=0; i<nstep; i++){
				c1 = dfex(r, i, n0, n1, n2, n3, nv1, nv2);
				//rho_new[i] = std::exp(c1+xi(rho,r,i,rho_b, phi_att_ff_int_ij,rho_phi_ff_int_i)/(kb1*T)-rhos_phi_sf_int_i[i]/(kb1*T)); // xi include kb1*T*(std::log(rho_b)) type.
				xio = c1+xi(rho,r,i,rho_b, phi_att_ff_int_ij,rho_phi_ff_int_i, phi_ext_i)/(kb1*T)-rhos_phi_sf_int_i[i]/(kb1*T);
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
			diff_old1 = diff;
			diff = 0.0;
#pragma omp parallel for
			for (i=0; i<=(nstep-2)/2; i++){
				diff0 = std::abs((rho_new[i]-rho[i])/rho[i]);
				diff = diff + 2.0*diff0;
				rho[i] = wmixing*rho_new[i] + (1.0-wmixing)*rho[i];
				rho[(nstep-1)-i] = rho[i]; // The rest is filled with mirror symmetry. 
			}
			if (diff < threshold && diff_old1 < threshold) {
				break;
			}
		}
		//
		v_gamma = integral_simpson(rho, nstep-1, dr);
		//v_gamma = v_gamma/(H-sigma_ss) - rho_b; // for NLDFT
		v_gamma = v_gamma/(H-(2.0*ze)) - rho_b;
		//v_mmol_per_cm3 = v_gamma * (1e7 * 1e7 * 1e7) / (6.02214076 * 1e23) * 1e3; // [mmol/cm3]
		//v_mmol_per_cm3 = (v_gamma / 6.02214076 ) * (1e24 / 1e23); // [mmol/cm3]
		v_mmol_per_cm3 = (v_gamma / 6.02214076) * 10; // [mmol/cm3]
		//v_cm3STP_per_cm3 = v_mmol_per_cm3 / 22.414 / (rho_ss*12.0107*(1e7*1e7*1e7)/(6.02214076*1e23)); // [cm3(STP)/g], 2.226 [g/cm3]
		v_cm3STP_per_cm3 = v_mmol_per_cm3 * 22.414;
		if (v_gamma < 0) { v_gamma = 0.0; }
		//v_gamma = v_gamma * (0.8064/28.0134/1e21*6.02214e23)/rho_b;
		// N2(77K): 0.8064 g/mL, 0.8064/28.0134 mol/mL, 0.8064/28.0134/1e21 mol/nm3, 0.8064/28.0134/1e21*6.02214e23 molecules/nm3
		//std::cout << "V= " << v_gamma << std::endl;
		//
		// press_hs(rho_b) from Carnahan-Starling (CS) equation of state
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		if(flag_P==0.0){
			pp0 = press_b/press_b0;
		} else if (flag_P<=-10.0){
			pp0 = press_b*kb*1e27/p0;
		} else {
			// kb1=1, kb = 1.38e-23 [J/K], T [K], rho_b [N/nm^3], 1 [atm] = 101325 [Pa]
			pp0 = press_b*kb*1e27;
		}
		// 1 [K/nm3] = 1.38064878e-23/(10^-9)^3 [Nm/m3] = 13806.4878 [Pa]
		// grand ppotential, Omega(ff+sf part)
		for (i=0; i<nstep; i++){
			fex_i[i] = fex(i, n0, n1, n2, n3, nv1, nv2);
		}
		grand_potential = omega(rho, r, fex_i, rho_phi_ff_int_i, rhos_phi_sf_int_i, rho_b, rhos, rhos_phi_ss_int_i);
		//grand_potential = 1.0;
		//std::cout << "P/P0= " << pp0 << std::endl;
		ofsppov_ls << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
		std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
	}
	free(phi_att_ff_int_ij);
MPI::Finalize();
	return 0;
}
