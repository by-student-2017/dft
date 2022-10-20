#include <stdio.h>  
#include <fstream>   // for file in and out (std::ofstream)
#include <string.h>  // strtok
#include <cmath>     // for log, exp
#include <iostream>  // for cout

// c++ make_eq.cpp -o make_eq
// ./make_eq

// global parametrs
//float kb = 1.38e-23; //[J/K] (8.61733262e-5 [eV/K])
float kb1 = 1.0;
float T  = 77.347; //[K]
int nsize = 2000; // array size

int main(void) {
	int ndata; // number of data
	//
	//ads
    FILE *ads_fin = fopen("./data_vs.txt", "rt");
    if (!ads_fin) {
        perror("fopen");
        return 1;
    }

	char str[256 + 1];
	float data1, data2, data3, data4, data5;
    float ads_pp0[nsize];
	float ads_v_gamma[nsize];
	float ads_v_mmol_per_cm3[nsize];
	float ads_v_cm3STP_per_cm3[nsize];
	float ads_grand_potential[nsize];

	fgets(str, 256, ads_fin); // 1st line
	fgets(str, 256, ads_fin); // 2nd line
	int i=0;
	while (fgets(str, 256, ads_fin) != NULL) {
		//printf("%s",str);
		if (sscanf(str, "%e,%f,%f,%f,%f", &data1, &data2, &data3, &data4, &data5) != 5) {
			continue;
		}
		//printf("%e,%f,%f,%f,%f \n", data1, data2, data3, data4, data5);
		//
		ads_pp0[i] = data1;
		ads_v_gamma[i] = data2;
		ads_v_mmol_per_cm3[i] = data3;
		ads_v_cm3STP_per_cm3[i] = data4;
		ads_grand_potential[i] = data5;
		//printf("%e,%f,%f,%f,%f \n", ads_pp0[i], ads_v_gamma[i], ads_v_mmol_per_cm3[i], ads_v_cm3STP_per_cm3[i], ads_grand_potential[i]);
		//
		i++;
	}
	ndata = i;
	//printf("Number of data = %d \n", ndata);
    fclose(ads_fin);
	//-----------------------------------------------------------------------------------
	//
	//des
    FILE *des_fin = fopen("./data_ls.txt", "rt");
    if (!des_fin) {
    	
    	
        perror("fopen");
        return 1;
    }

	char str1[256 + 1];
	char str2[256 + 1];
	//float data1, data2, data3, data4, data5;
    float des_pp0[nsize];
	float des_v_gamma[nsize];
	float des_v_mmol_per_cm3[nsize];
	float des_v_cm3STP_per_cm3[nsize];
	float des_grand_potential[nsize];

	fgets(str1, 256, des_fin); // 1st line
	fgets(str2, 256, des_fin); // 2nd line
	strtok(str2, "\n\0");
	//int i=0;
	i=ndata-1;
	while (fgets(str, 256, des_fin) != NULL) {
		//printf("%s",str);
		if (sscanf(str, "%e,%f,%f,%f,%f", &data1, &data2, &data3, &data4, &data5) != 5) {
			continue;
		}
		//printf("%e,%f,%f,%f,%f \n", data1, data2, data3, data4, data5);
		//
		des_pp0[i] = data1;
		des_v_gamma[i] = data2;
		des_v_mmol_per_cm3[i] = data3;
		des_v_cm3STP_per_cm3[i] = data4;
		des_grand_potential[i] = data5;
		//printf("%e,%f,%f,%f,%f \n", des_pp0[i], des_v_gamma[i], des_v_mmol_per_cm3[i], des_v_cm3STP_per_cm3[i], des_grand_potential[i]);
		//
		i--;
	}
    fclose(des_fin);
	//-----------------------------------------------------------------------------------
	//
	//equilibrium
	float Ap;
	std::ofstream ofsppov_eq("./data_eq.txt");
	ofsppov_eq << str1 << str2 << std::endl;
	float pp0, v_gamma, v_mmol_per_cm3, v_cm3STP_per_cm3, grand_potential;
	for (i=0; i<ndata; i++){
		if (ads_grand_potential[i] <= des_grand_potential[i]) {
			//pp0 = ads_pp0[i];
			//v_gamma = ads_v_gamma[i];
			//v_mmol_per_cm3 = ads_v_mmol_per_cm3[i];
			//v_cm3STP_per_cm3 = ads_v_cm3STP_per_cm3[i];
			//grand_potential = ads_grand_potential[i];
			Ap = 1.0/(1.0+std::exp(-(des_grand_potential[i]-ads_grand_potential[i])/(kb1*T)));
			pp0 = ads_pp0[i];
			v_gamma = ads_v_gamma[i]*Ap + des_v_gamma[i]*(1.0-Ap);
			v_mmol_per_cm3 = ads_v_mmol_per_cm3[i]*Ap + des_v_mmol_per_cm3[i]*(1.0-Ap);
			v_cm3STP_per_cm3 = ads_v_cm3STP_per_cm3[i]*Ap + des_v_cm3STP_per_cm3[i]*(1.0-Ap);
			grand_potential = ads_grand_potential[i]*Ap + des_grand_potential[i]*(1.0-Ap);
			Ap = 1.0;
		} else {
			//pp0 = des_pp0[i];
			//v_gamma = des_v_gamma[i];
			//v_mmol_per_cm3 = des_v_mmol_per_cm3[i];
			//v_cm3STP_per_cm3 = des_v_cm3STP_per_cm3[i];
			//grand_potential = des_grand_potential[i];
			Ap = 1.0/(1.0+std::exp(-(ads_grand_potential[i]-des_grand_potential[i])/(kb1*T)));
			pp0 = des_pp0[i];
			v_gamma = des_v_gamma[i]*Ap + ads_v_gamma[i]*(1.0-Ap);
			v_mmol_per_cm3 = des_v_mmol_per_cm3[i]*Ap + ads_v_mmol_per_cm3[i]*(1.0-Ap);
			v_cm3STP_per_cm3 = des_v_cm3STP_per_cm3[i]*Ap + ads_v_cm3STP_per_cm3[i]*(1.0-Ap);
			grand_potential = des_grand_potential[i]*Ap + ads_grand_potential[i]*(1.0-Ap);
			Ap = 0.0;
		}
		ofsppov_eq << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << std::endl;
		//std::cout << pp0 << ", "<< v_gamma << ", " << v_mmol_per_cm3 << ", " <<  v_cm3STP_per_cm3 << ", " << grand_potential << ", " << Ap << std::endl;
	}

    return 0;
}