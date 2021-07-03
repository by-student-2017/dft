#include <iostream>
#include <cmath>

//using namespace std;

//non-local smoothed density approximation：SDA
//non-local density functional theory（NLDFT)
//Reference: https://www.j-ad.org/adsorption_news/30_1.pdf

// compiling: c++ nldft.cpp
// usage: ./a.out

// Adsorbent 
double H = 1.00; //distace of slit [nm]
double dH = 0.01; //dH = dr = dz in this case, because of one dimension calculation.
unsigned int nstep = H/dH;
double sigma_ss = 0.34; // [nm]

// iteration of rho
unsigned int cycle_max = 20;

//Carbon dioxide 253.9  [K](epsilon), 0.3454 [nm](sigma), 0.3495 [nm](d_hs)
//Argon          118.05 [K](epsilon), 0.3305 [nm](sigma), 0.3390 [nm](d_hs)
//Nitrogen        94.45 [K](epsilon), 0.3575 [nm](sigma), 0.3575 [nm](d_hs)
double epsilon_ff = 94.45;
double sigma_ff = 0.3575;
double d_hs = 0.3575;
double rc = 12.8; // cut off
double rm = std::pow(2.0,1.0/6.0)*sigma_ff; //minimum position of LJ

// Carbon dioxide/Carbon slit 81.5  [K](epsilon), 0.3430 [nm](sigma)
// Nitrogen/Carbon slit       53.72 [K](epsilon), 0.3508 [nm](sigma)
double epsilon_sf = 53.72; // [K] 
double sigma_sf = 0.3508; // [nm]

// slit pore (graphite)
double delta = 0.335; // [nm]
double rho_ss = 11.4; // [nm^-3]

//double m = 14.0067*2.0/(6.02214076e23)/1000; // N2 = 4.65173e-26 [kg]
double m = 4.65173e-26; //[kg] (N2) (e.g., Ar = 6.63e-26 [kg])
double k = 1.0;
double kb = 1.38e-23; //[J/K] (8.61733262e-5 [eV/K])
double T = 77.347; //[K]
double h = 6.63e-34; //[Js] (4.135667696e-15 [eVs])
// thermal de Broglie wavelength
double lam = h/std::pow((2.0*M_PI*m*kb*T),0.5)*1e9; //[nm]
// Ref: https://www1.doshisha.ac.jp/~bukka/lecture/statistic/pdftext/std-07.pdf

// alpha = integal phi_att * -1.0
double alpha = (32.0/9.0)*M_PI*epsilon_ff*std::pow(rm,3.0) - (16.0/9.0)*M_PI*epsilon_ff*std::pow(sigma_ff,3.0)*
	( 3.0*std::pow((sigma_ff/rc),3.0) - std::pow((sigma_ff/rc),9.0) );

//Barker-Henderson model
double d_bh(void){
	double epsilon = 94.45; //dummy
	double sigma = 0.3575; //dummy
	//Lstoskie et al.,
	double xi1 = 0.3837;
	double xi2 = 1.035;
	double xi3 = 0.4249;
	double xi4 = 1.0;
	double d_bh_out;
	d_bh_out = (xi1*k*T/epsilon+xi2)/(xi3*k*T/epsilon+xi4)*sigma;
	return d_bh_out;
}

double phi_att(double r){
        double e;
	// WCA (Weeks-Chandler-Anderson) type
	if (r < rm){
		e = -1.0*epsilon_ff;
	}else if (rm <= r && r < rc){
		// Lennard-Jones（LJ) potential
		e = 4.0*epsilon_ff*( std::pow((sigma_ff/r),12.0) - std::pow((sigma_ff/r),6.0) );
	}else if (r < rc){
		e = 0.0;
	}
	return e;
}

// Tarazona theory
double wi(double r, int i){
	double wi_out;
	switch(i){
		case 0:
			if (r <= d_hs){
				wi_out = 3.0/(4.0*M_PI*std::pow(d_hs,3.0));
			} else {
				wi_out = 0.0;
			}
			break;
		case 1:
			if (r <= d_hs) {
				wi_out = 0.475-0.648*(r/d_hs)+0.113*std::pow((r/d_hs),2.0);
			} else if (d_hs < r && r <= 2.0*d_hs) {
				wi_out = 0.288*(d_hs/r)-0.924+0.764*(r/d_hs)-0.187*std::pow((r/d_hs),2.0);
			} else {
				wi_out = 0.0;
			}
			break;
		case 2:
			if (r <= d_hs) {
				wi_out = 5.0*M_PI*std::pow(d_hs,3.0)/144.0 * (6.0-12.0*(r/d_hs)+5.0*std::pow((r/d_hs),2.0));
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

double rho_si(double *rho, double r1, double *r, int i){
	unsigned int j;
	double rho_s_out;
	rho_s_out = 0.0;
	for (j=0; j<nstep; j++) {
		rho_s_out = rho_s_out + rho[j]*wi(std::abs(r1-r[j]),i);
	}
	return rho_s_out;
}

double rho_s(double *rho, double r1, double *r){
	double rho_den1, rho_den2, rho_s_out;
	rho_den1 = std::pow((1.0 - rho_si(rho,r1,r,1)),2.0);
	rho_den2 = std::pow((rho_den1 - 4.0*rho_si(rho,r1,r,0)*rho_si(rho,r1,r,2)),0.5);
	rho_s_out = 2.0*rho_si(rho,r1,r,0)/(1.0 - rho_si(rho,r1,r,1)+rho_den2);
	return rho_s_out;
}

double f_ex(double rho_s){
        double eta, f_ex_out;
	eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
	f_ex_out = k*T*eta*(4.0-3.0*eta)/std::pow((1.0-eta),2.0);
	return f_ex_out;
}

//Steele 10-4-3 potential
double phi_sf(double z){
	double phi_sf_out;
	phi_sf_out = 2.0*M_PI*rho_ss*epsilon_sf*std::pow(sigma_sf,2.0)*
		delta*( (2.0/5.0)*std::pow((sigma_sf/z),10.0)-std::pow((sigma_sf/z),4.0)-std::pow(sigma_sf,4.0)/
			(3.0*delta*std::pow((0.61*delta+z),3.0)) );
	return phi_sf_out;
}

//e.g., Carbon slit
double phi_ext(double z){
	double phi_ext_out;
	phi_ext_out = phi_sf(z) + phi_sf(H-z);
	return phi_ext_out;
}

double mu_ex(double rho_b){
	double y, mu_ex_out;
	y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	mu_ex_out = k*T*(8.0*y-9.0*y*y+3.0*y*y*y)/std::pow((1.0-y),3.0);
	return mu_ex_out;
}

double mu_b(double rho_b){
	double mu_id, mu_hs, mu_b_out;
	mu_id = k*T*std::log(std::pow(lam,3.0)*rho_b);
	mu_hs = mu_id + mu_ex(rho_b);
	mu_b_out = mu_hs + rho_b*-1.0*alpha; 
	return mu_b_out;
}

// d(f_ex)/d(rho_s)
double dfex_per_drhos(double rho_s){
        double dfex_per_drhos_out;
	double eta;
	eta = M_PI*rho_s*std::pow(d_hs,3.0)/6.0;
	dfex_per_drhos_out = k*T*(4.0-18.0*eta+12.0*eta*eta)*M_PI*std::pow(d_hs,3.0)/6.0;
	return dfex_per_drhos_out;
}

// d(rho_s)/d(rho)
double drhos_per_drho(double *rho, double r1, double r2, double *r){
	double w, drhos_per_drho_out;
	w = wi(std::abs(r1-r2),0) + wi(std::abs(r1-r2),1)*rho_s(rho,r1,r) + wi(std::abs(r1-r2),2)*std::pow(rho_s(rho,r1,r),2.0);
	drhos_per_drho_out = w/(1.0-rho_si(rho,r2,r,1)-2.0*rho_si(rho,r2,r,2)*rho_s(rho,r2,r));
	return drhos_per_drho_out;
}

// Grand potential Omega
// Euler-Lagrange equation d(Omega)/d(rho) = 0 at mu = mu_b
double xi(double *rho, int i, double rho_b, double *r){
	unsigned int j;
	double rho_fex_int, rho_phi_int;
	double xi_out;
	rho_fex_int = 0.0;
	rho_phi_int = 0.0;
	for (j=0; j<nstep; j++) {
		// d(f_ex)/d(rho) = d(f_ex)/d(rho_s) * d(rho_s)/d(rho)
        	rho_fex_int = rho_fex_int + rho[j]*dfex_per_drhos(rho_s(rho,r[j],r))*drhos_per_drho(rho,r[i],r[j],r)*dH;
		rho_phi_int = rho_phi_int + rho[j]*phi_att(std::abs(r[i]-r[j]))*dH;
		//std::cout << dfex_per_drhos(rho_s(rho,r[j],r)) << ", " << drhos_per_drho(rho,r[i],r[j],r) << std::endl;
	}
	xi_out = mu_ex(rho_b) - rho_b*alpha - phi_ext(r[i]) - f_ex(rho_s(rho,r[i],r)) - rho_fex_int - rho_phi_int;
	//std::cout << f_ex(rho_s(rho,r[i],r)) << ", " << rho_fex_int << ", " << rho_phi_int << std::endl;
	return xi_out;
}

double press_hs(double rho_b){
	double y, press_hs_out;
	y = M_PI*rho_b*std::pow(d_hs,3.0)/6.0;
	press_hs_out = rho_b*k*T * (1.0 + y + y*y - y*y*y)/std::pow((1.0-y),3.0);
	return press_hs_out;
}

double Maxwell_equal_area_rule(void){
	unsigned int i,j;
	unsigned int iter_max_drhob0 = 1000000;
	unsigned int iter_max_dmue = 5000;
	double drhob0 = 0.000001;
	double dmue = 0.01;
	double threshold_diff = 0.1;
	double threshold_find = 0.1;
	//
	double mu_b_per_epsilon_ff[iter_max_drhob0];
	double mu_e_per_epsilon_ff,diff,diffp;
	int flag;
	double press_b0;
	double rho_b0;
	// rho_b vs. mu_b/epsilon_ff
	for (i=0; i<iter_max_drhob0; i++){
		rho_b0 = drhob0*double(i+1.0);
		mu_b_per_epsilon_ff[i] = mu_b(rho_b0)/epsilon_ff;
		//std::cout << "rho_b0 = "<< rho_b0 << ", mu_b/epsilon_ff = " << mu_b_per_epsilon_ff[i] << std::endl;
	}
	// Maxwell equal area rule
	for (j=0; j<iter_max_dmue; j++){
		mu_e_per_epsilon_ff = dmue*double(j+1.0) - 25.0;
		diff = 0.0;
		flag = 0;
		for (i=0; i<iter_max_drhob0; i++){
			diffp = mu_b_per_epsilon_ff[i] - mu_e_per_epsilon_ff;
			if (diffp > 0.0 && flag != 1) {
				diff = diff + diffp*drhob0;
			} else if (diffp <= 0.0) {
				diff = diff + diffp*drhob0;
				flag = 1;
			}
			//std::cout << diffp << std::endl;
		}
		//std::cout << "mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << ", diff = " << diff << std::endl;
		if (std::abs(diff) <= threshold_diff) {
			//std::cout << "mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << ", diff = " << diff << std::endl;
			//break;
		}
	}
	// find rho_b0
	for (i=0; i<iter_max_drhob0; i++){
		rho_b0 = drhob0*double(i+1.0);
		if (std::abs(mu_b(rho_b0)/epsilon_ff - mu_e_per_epsilon_ff) <= threshold_find) {
			//std::cout << rho_b0 << std::endl;
			//break;
		}
	}
	//std::cout << "mu_e/epsilon_ff = " << mu_e_per_epsilon_ff << std::endl;
	//std::cout << "rho_b0 = " << rho_b0 << std::endl;
	return rho_b0;
}

int main(){
	unsigned int i,j,k;
	double w = 0.3;
	double r[nstep];
	double mu_e,diff,diffp;
	int flag;
	double rho[nstep], rho_old[nstep];
	double v_gamma;
	double press_b, press_b0, pp0;
	double rho_b, rho_b0;
	// set dr
	for (i=0; i<nstep; i++){
		r[i] = sigma_ss/2.0 + (H-sigma_ss/2.0)/nstep*i;
		//std::cout << i << ", " << r[i] << std::endl;
	}
	// set rho_b0
	//rho_b0 = Maxwell_equal_area_rule();
	rho_b0 = 0.50;
	// initialization
	for (i=0; i<nstep; i++){
		rho[i] = 1e-6/nstep;
		rho_old[i] = 1e-6/nstep;
	}
	// volume and pressure
	for (k=1; k<=100; k++){
		rho_b = rho_b0 * std::exp(-double(20-2*k/10));
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "rho_b = " << rho_b << std::endl;
		for (j=0; j<=cycle_max; j++){
			for (i=0; i<nstep; i++){
				rho[i] = rho_b*std::exp(xi(rho,i,rho_b,r)/(k*T));
				rho[i] = w*rho[i] + (1.0-w)*rho_old[i];
				rho_old[i] = rho[i];
			}
			//std::cout << rho_b << ", " << j << ", " << rho[nstep/2] << std::endl;
		}
		//
		v_gamma = 0.0;
		for (i=0; i<nstep; i++){
			//std::cout << r[i] << ", " << rho[i] << std::endl;
			v_gamma = v_gamma + rho[i];
		}
		v_gamma = v_gamma/(H-sigma_ss) - rho_b;
		std::cout << "V= " << v_gamma << std::endl;
		//
		press_b = press_hs(rho_b) - 0.5*std::pow(rho_b,2.0)*alpha;
		press_b0 = press_hs(rho_b0) - 0.5*std::pow(rho_b0,2.0)*alpha;
		std::cout << "P= " << press_b << std::endl;
		pp0 = press_b/press_b0;
		//std::cout << "P/P0= " << pp0 << std::endl;
	}
        return 0;
}
