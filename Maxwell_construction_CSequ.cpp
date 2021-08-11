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

// compiling: c++ nldft.cpp -O2
// usage: ./a.out

// debag mode
// compiling: c++ nldft.cpp -g -Wall -O0
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
	rho_b0 = num[16];
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

// The excess hard sphere chemical potential (mu_ex) in the bulk fulid.
// mu_ex is calculated by the PY equation.
//double mu_ex(double rho_b){
//	double y, mu_ex_out;
//	//eta = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
//	//eta = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
//	y = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
//	//double den1e = (1.0-eta);
//	double den1y = (1.0-y);
//	//mu_ex_out = kb1*T*(-std::log(1-eta) + eta*(14.0 - 13.0*eta + 5.0*eta*eta)/(2.0*std::pow((1.0-eta),3.0)));
//	//mu_ex_out = kb1*T*(-std::log(den1e) + eta*(14.0 - 13.0*eta + 5.0*eta*eta)/(2.0*(den1e*den1e*den1e)));
//	mu_ex_out = kb1*T*(-std::log(den1y) + y*(14.0 - 13.0*y + 5.0*y*y)/(2.0*den1y*den1y*den1y));
//	return mu_ex_out;
//}

double mu_b(double rho_b){
	double mu_id, mu_hs, mu_b_out;
	//mu_id = kb1*T*std::log(std::pow(lam,3.0)*rho_b);
	mu_id = kb1*T*std::log((lam*lam*lam)*rho_b);
	mu_hs = mu_id + mu_ex(rho_b);
	mu_b_out = mu_hs - rho_b*alpha;
	return mu_b_out;
}

// For SDA, from Carnahan-Starling (CS) equation
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

// For FMT, from Percus-Yevick (PY) equation
// http://www.sklogwiki.org/SklogWiki/index.php/Exact_solution_of_the_Percus_Yevick_integral_equation_for_hard_spheres
//double press_hs(double rho_b){
//	double eta, press_hs_out;
//	//eta = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
//	eta = M_PI*rho_b*(d_hs*d_hs*d_hs)/6.0;
//	double den1e = (1.0-eta);
//	//press_hs_out = rho_b*kb1*T* (1.0 + eta + eta*eta)/std::pow((1.0-eta),3.0);
//	//press_hs_out = rho_b*kb1*T* (1.0 + eta + eta*eta)/((1.0-eta)*(1.0-eta)*(1.0-eta));
//	press_hs_out = rho_b*kb1*T* (1.0 + eta + eta*eta)/(den1e*den1e*den1e);
//	return press_hs_out;
//}

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
	}
	std::cout << "--------------------------------------------------" << std::endl;
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
	std::cout << "--------------------------------------------------" << std::endl;
	return rho_b0_out;
}

int main(){
	
	read_parameters();
	
	rho_b0 = Maxwell_construction();
	
	return 0;
}
