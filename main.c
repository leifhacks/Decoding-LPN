#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MaxDepth 5
int ErrorFlag = 0, HDFlag = 0, NaiveFlag = 0, TypeFlag = 2, WorstCaseFlag = 1, NNFlag = 1, LPNFlag = 0, PreciseFlag = 1, MMT_Flag = 0, SingleFlag = 0;

#include "tools.h"
#include "decode.h"
#include <time.h>

int main(int argc, char *argv[]) {

	// initialization
	initH();
    clock_t tStart = clock(), tRound;
	time_t tNow;
	srand(time(NULL));
	FILE *file;
	int contFlag = 0;
	char *filename = "output.txt";
	double pmin=0, pmax=0, psteps=0.1;
    double kmin=0.0, kmax=0.5, ksteps=0.1;
	double wmin=0, wmax=0, wsteps=0.1;
    double emin[MaxDepth], emax[MaxDepth], esteps[MaxDepth];
    for (int i=0;i<MaxDepth;i++) {emin[i] = 0.0, emax[i] = 0.0, esteps[i] = 1.0;}
    double wwmin[MaxDepth][MaxDepth], wwmax[MaxDepth][MaxDepth], wwsteps[MaxDepth][MaxDepth];
    for (int i=0;i<MaxDepth;i++) {
		for (int j=0;j<MaxDepth;j++) {
			wwmin[i][j] = 0.0, wwmax[i][j] = 0.0, wwsteps[i][j] = 1.0;
		}
	}
    double lmin[MaxDepth], lmax[MaxDepth], lsteps[MaxDepth];
    for (int i=0;i<MaxDepth;i++) {lmin[i] = 0.0, lmax[i] = 0.0, lsteps[i] = 1.0;}
    double k=0, w=0;
    double spaceBound=1.0, Tmax=0.0, Tmin=1.0, MAXvalues[MaxDepth*8-1], MINvalues[MaxDepth*8-1], MINParams[3*MaxDepth-3], MAXParams[3*MaxDepth-3];
	for (int i=0;i<MaxDepth*8-1;i++) MAXvalues[i] = 0;
	for (int i=0;i<MaxDepth*8-1;i++) MINvalues[i] = 0;
	for (int i=0;i<3*MaxDepth-3;i++) MINParams[i] = 0;
	for (int i=0;i<3*MaxDepth-3;i++) MAXParams[i] = 0;
	
	// standard settings
	int depth = 2, LPN_k = 512;
	double LPN_t = 0.25;	
	
	// read input
	int wset=0;
	char *item[2];
    for (int i=1; i<argc; i++) {
		item[0] = strtok(argv[i], "=");
		item[1] = strtok(NULL, "=");
		
		if (strcmp(item[0], "m")==0) {
            depth = atoi(item[1]);
        } else if (strcmp(item[0],"kmin")==0) {
            kmin = atof(item[1]);
        } else if (strcmp(item[0],"kmax")==0) {
            kmax = atof(item[1]);
        } else if (strcmp(item[0],"ksteps")==0) {
            ksteps = atof(item[1]);
        } else if (strcmp(item[0], "Precise")==0) {
            PreciseFlag = atoi(item[1]);
        } else if (strcmp(item[0], "WC")==0) {
            WorstCaseFlag = atoi(item[1]);
        } else if (strcmp(item[0], "HD")==0) {
            HDFlag = atoi(item[1]);
        } else if (strcmp(item[0], "Naive")==0) {
            NaiveFlag = atoi(item[1]);
        } else if (strcmp(item[0], "NN")==0) {
            NNFlag = atoi(item[1]);
        } else if (strcmp(item[0], "Type")==0) {
            TypeFlag = atoi(item[1]);
        } else if (strcmp(item[0], "LPN")==0) {
            LPNFlag = atoi(item[1]);
        } else if (strcmp(item[0], "LPNk")==0) {
            LPN_k = atoi(item[1]);
        } else if (strcmp(item[0], "LPNtau")==0) {
            LPN_t = atof(item[1]);
        } else if (strcmp(item[0],"wmin")==0) {
            wmin = atof(item[1]);
			wset=1;
        } else if (strcmp(item[0],"wmax")==0) {
            wmax = atof(item[1]);
			wset=1;
        } else if (strcmp(item[0], "wsteps")==0) {
            wsteps = atof(item[1]);
			wset=1;
        } else if (strcmp(item[0], "Output")==0) {
			filename = item[1];
        } else if (strcmp(item[0], "Space")==0) {
			spaceBound = atof(item[1]);
        } else if (strcmp(item[0], "MMT")==0) {
			MMT_Flag = atoi(item[1]);
		} else if (strcmp(item[0], "Single")==0) {	
			SingleFlag = atoi(item[1]);
        }
    }
	
	// create file
	file = fopen(filename, "wb");
	fclose(file);
	
	// initialization 2
	if (WorstCaseFlag == 0) Tmax=1.0;
	if (LPNFlag == 1) {
		Tmax *= LPN_k / kmin;
		if (wset == 0) wmin = wmax = LPN_t;
		if (spaceBound <= 1.0) spaceBound *= LPN_k;
		spaceBound *= k / LPN_k;
	}
	
	// Optimization for every k, w
	for (w=wmin; w<=wmax + wsteps*0.1; w += wsteps) {
	for (k=kmin; k<=kmax + ksteps*0.1; k += ksteps) {
		
		// Our algorithm
		if (TypeFlag==2) {
			if (depth==2) {
				pmin=0.00, pmax=0.06, psteps=0.01;
				wwmin[1][0]=0.01, wwmax[1][0]=0.07, wwsteps[1][0]=0.01;
				lmin[1]=0.13, lmax[1]=0.19, lsteps[1]=0.01;
				emin[1]=0.004, emax[1]=0.010, esteps[1]=0.001;
			// p[0]=0.0340 w[1]=0.0071 w[2]=0.0098 l[1]=0.0818 l[2]=0.0517 e[1]=0.00393 e[2]=0.00030	
			} else if (depth==3) {
				pmin=0.00, pmax=0.06, psteps=0.01;
				wwmin[1][0]=0.00, wwmax[1][0]=0.06, wwsteps[1][0]=0.01;
				wwmin[2][0]=0.00, wwmax[2][0]=0.06, wwsteps[2][0]=0.01;
				lmin[1]=0.00, lmax[1]=0.13, lsteps[1]=0.01;
				lmin[2]=0.00, lmax[2]=0.08, lsteps[2]=0.01;
				emin[1]=0.000, emax[1]=0.007, esteps[1]=0.001;
				emin[2]=0.000, emax[2]=0.006, esteps[2]=0.001;
			// p[0]=0.0347 w[1]=0.0119 w[2]=0.0099 w[3]=0.0066 l[1]=0.0932 l[2]=0.0549 l[3]=0.0371 e[1]=0.00296 e[2]=0.00060 e[3]=0.00024 
			} else if (depth==4) {
				pmin=0.00, pmax=0.06, psteps=0.01;
				wwmin[1][0]=0.00, wwmax[1][0]=0.06, wwsteps[1][0]=0.01;
				wwmin[2][0]=0.00, wwmax[2][0]=0.06, wwsteps[2][0]=0.01;
				wwmin[3][0]=0.00, wwmax[3][0]=0.06, wwsteps[3][0]=0.01;
				lmin[1]=0.06, lmax[1]=0.12, lsteps[1]=0.01;
				lmin[2]=0.02, lmax[2]=0.08, lsteps[2]=0.01;
				lmin[3]=0.00, lmax[3]=0.06, lsteps[3]=0.01;
				emin[1]=0.000, emax[1]=0.006, esteps[1]=0.001;
				emin[2]=0.000, emax[2]=0.006, esteps[2]=0.001;
				emin[3]=0.000, emax[3]=0.006, esteps[3]=0.001;
			} else if (depth==5) {
				pmin=0.00, pmax=0.06, psteps=0.01;
				wwmin[1][0]=0.00, wwmax[1][0]=0.06, wwsteps[1][0]=0.01;
				wwmin[2][0]=0.00, wwmax[2][0]=0.06, wwsteps[2][0]=0.01;
				wwmin[3][0]=0.00, wwmax[3][0]=0.06, wwsteps[3][0]=0.01;
				wwmin[4][0]=0.00, wwmax[4][0]=0.06, wwsteps[4][0]=0.01;
				lmin[1]=0.06, lmax[1]=0.12, lsteps[1]=0.01;
				lmin[2]=0.02, lmax[2]=0.08, lsteps[2]=0.01;
				lmin[3]=0.00, lmax[3]=0.06, lsteps[3]=0.01;
				lmin[4]=0.00, lmax[4]=0.06, lsteps[4]=0.01;
				emin[1]=0.000, emax[1]=0.006, esteps[1]=0.001;
				emin[2]=0.000, emax[2]=0.006, esteps[2]=0.001;
				emin[3]=0.000, emax[3]=0.006, esteps[3]=0.001;
				emin[4]=0.000, emax[4]=0.006, esteps[4]=0.001;
			}
			
		// BJMM(d)
		} else if (TypeFlag==1) {
			if (depth==2) {
				pmin=0.00, pmax=0.07, psteps=0.01;
				lmin[1]=0.00, lmax[1]=0.12, lsteps[1]=0.01;
				emin[1]=0.000, emax[1]=0.010, esteps[1]=0.001;	   
			} else if (depth==3) {
				pmin=0.00, pmax=0.09, psteps=0.01;
				lmin[1]=0.00, lmax[1]=0.30, lsteps[1]=0.01;
				emin[1]=0.000, emax[1]=0.025, esteps[1]=0.001;
				emin[2]=0.000, emax[2]=0.007, esteps[2]=0.001;
			} else if (depth>=4) {
				pmin=0.00, pmax=0.10, psteps=0.01;
				lmin[1]=0.00, lmax[1]=0.28, lsteps[1]=0.01;
				emin[1]=0.000, emax[1]=0.035, esteps[1]=0.001;
				emin[2]=0.000, emax[2]=0.018, esteps[2]=0.001;
				emin[3]=0.000, emax[3]=0.006, esteps[3]=0.001;
			}
		}
		
		if (LPNFlag == 1) {
			double xx=0.1;
			pmin*=xx, pmax*=xx, psteps*=xx;
			wwmin[1][0]*=xx, wwmax[1][0]*=xx, wwsteps[1][0]*=xx;
			wwmin[2][0]*=xx, wwmax[2][0]*=xx, wwsteps[2][0]*=xx;
			wwmin[3][0]*=xx, wwmax[3][0]*=xx, wwsteps[3][0]*=xx;
			wwmin[4][0]*=xx, wwmax[4][0]*=xx, wwsteps[4][0]*=xx;
			lmin[1]*=xx, lmax[1]*=xx, lsteps[1]*=xx;
			lmin[2]*=xx, lmax[2]*=xx, lsteps[2]*=xx;
			lmin[3]*=xx, lmax[3]*=xx, lsteps[3]*=xx;
			lmin[4]*=xx, lmax[4]*=xx, lsteps[4]*=xx;
			emin[1]*=xx, emax[1]*=xx, esteps[1]*=xx;
			emin[2]*=xx, emax[2]*=xx, esteps[2]*=xx;
			emin[3]*=xx, emax[3]*=xx, esteps[3]*=xx;
			emin[4]*=xx, emax[4]*=xx, esteps[4]*=xx;
		}
		
		if (MMT_Flag == 1) {
			for (int i=0;i<6;i++) emin[i]=emax[i]=0.0;
		}
		
/* 		// 0.46: T:0.087284 p[0]=0.0360 w[1]=0.0180 w[2]=0.0089 w[3]=0.0100 w[4]=0.0110 l[1]=0.1260 l[2]=0.0499 l[3]=0.0590 l[4]=0.0600 
		// e[1]=0.00400 e[2]=0.00100 e[3]=0.00100 e[4]=0.00300 llww
		if (k == 0.46) {
			pmin=0.036, pmax=0.036, psteps=0.0001;
			wwmin[1][0]=0.0180, wwmax[1][0]=0.0180, wwsteps[1][0]=0.0001;
			wwmin[2][0]=0.0089, wwmax[2][0]=0.0089, wwsteps[2][0]=0.0001;
			wwmin[3][0]=0.0100, wwmax[3][0]=0.0100, wwsteps[3][0]=0.0001;
			wwmin[4][0]=0.0110, wwmax[4][0]=0.0110, wwsteps[4][0]=0.0001;
			lmin[1]=0.1260, lmax[1]=0.1260, lsteps[1]=0.0001;
			lmin[2]=0.0499, lmax[2]=0.0499, lsteps[2]=0.0001;
			lmin[3]=0.0590, lmax[3]=0.0590, lsteps[3]=0.0001;
			lmin[4]=0.0600, lmax[4]=0.0600, lsteps[4]=0.0001;
			emin[1]=0.0040, emax[1]=0.0040, esteps[1]=0.0001;
			emin[2]=0.0010, emax[2]=0.0010, esteps[2]=0.0001;
			emin[3]=0.0010, emax[3]=0.0010, esteps[3]=0.0001;
			emin[4]=0.0030, emax[4]=0.0030, esteps[4]=0.0001;
		
		// 0.45: T:0.087155 p[0]=0.0360 w[1]=0.0180 w[2]=0.0095 w[3]=0.0110 w[4]=0.0110 l[1]=0.1240 l[2]=0.0502 l[3]=0.0590 l[4]=0.0600 
		// e[1]=0.00400 e[2]=0.00100 e[3]=0.00100 e[4]=0.00300 llw
		} else if (k == 0.45) {
			pmin=0.036, pmax=0.036, psteps=0.0001;
			wwmin[1][0]=0.0180, wwmax[1][0]=0.0180, wwsteps[1][0]=0.0001;
			wwmin[2][0]=0.0095, wwmax[2][0]=0.0095, wwsteps[2][0]=0.0001;
			wwmin[3][0]=0.0110, wwmax[3][0]=0.0110, wwsteps[3][0]=0.0001;
			wwmin[4][0]=0.0110, wwmax[4][0]=0.0110, wwsteps[4][0]=0.0001;
			lmin[1]=0.1240, lmax[1]=0.1240, lsteps[1]=0.0001;
			lmin[2]=0.0502, lmax[2]=0.0502, lsteps[2]=0.0001;
			lmin[3]=0.0590, lmax[3]=0.0590, lsteps[3]=0.0001;
			lmin[4]=0.0600, lmax[4]=0.0600, lsteps[4]=0.0001;
			emin[1]=0.0040, emax[1]=0.0040, esteps[1]=0.0001;
			emin[2]=0.0010, emax[2]=0.0010, esteps[2]=0.0001;
			emin[3]=0.0010, emax[3]=0.0010, esteps[3]=0.0001;
			emin[4]=0.0030, emax[4]=0.0030, esteps[4]=0.0001;
		// 0.47: T:0.086675 p[0]=0.0350 w[1]=0.0130 w[2]=0.0089 w[3]=0.0100 w[4]=0.0110 l[1]=0.0970 l[2]=0.0526 l[3]=0.0570 l[4]=0.0590 
		// e[1]=0.00300 e[2]=0.00100 e[3]=0.00100 e[4]=0.00300 e+
		} else if (k == 0.47) {
			pmin=0.035, pmax=0.035, psteps=0.0001;
			wwmin[1][0]=0.0130, wwmax[1][0]=0.0130, wwsteps[1][0]=0.0001;
			wwmin[2][0]=0.0089, wwmax[2][0]=0.0089, wwsteps[2][0]=0.0001;
			wwmin[3][0]=0.0100, wwmax[3][0]=0.0100, wwsteps[3][0]=0.0001;
			wwmin[4][0]=0.0110, wwmax[4][0]=0.0110, wwsteps[4][0]=0.0001;
			lmin[1]=0.0970, lmax[1]=0.0970, lsteps[1]=0.0001;
			lmin[2]=0.0526, lmax[2]=0.0526, lsteps[2]=0.0001;
			lmin[3]=0.0570, lmax[3]=0.0570, lsteps[3]=0.0001;
			lmin[4]=0.0590, lmax[4]=0.0590, lsteps[4]=0.0001;
			emin[1]=0.0030, emax[1]=0.0030, esteps[1]=0.0001;
			emin[2]=0.0010, emax[2]=0.0010, esteps[2]=0.0001;
			emin[3]=0.0010, emax[3]=0.0010, esteps[3]=0.0001;
			emin[4]=0.0030, emax[4]=0.0030, esteps[4]=0.0001;
		}
		
		// T:0.086659 p[0]=0.0360 w[1]=0.0168 w[2]=0.0087 w[3]=0.0108 w[4]=0.0109 l[1]=0.1242 l[2]=0.0485 l[3]=0.0608 l[4]=0.0596 
		// e[1]=0.00400 e[2]=0.00080 e[3]=0.00100 e[4]=0.00260
		pmin=0.036, pmax=0.036, psteps=0.0001;
		wwmin[1][0]=0.0168, wwmax[1][0]=0.0168, wwsteps[1][0]=0.0001;
		wwmin[2][0]=0.0087, wwmax[2][0]=0.0087, wwsteps[2][0]=0.0001;
		wwmin[3][0]=0.0108, wwmax[3][0]=0.0108, wwsteps[3][0]=0.0001;
		wwmin[4][0]=0.0109, wwmax[4][0]=0.0109, wwsteps[4][0]=0.0001;
		lmin[1]=0.1242, lmax[1]=0.1242, lsteps[1]=0.0001;
		lmin[2]=0.0485, lmax[2]=0.0485, lsteps[2]=0.0001;
		lmin[3]=0.0608, lmax[3]=0.0608, lsteps[3]=0.0001;
		lmin[4]=0.0596, lmax[4]=0.0596, lsteps[4]=0.0001;
		emin[1]=0.0040, emax[1]=0.0040, esteps[1]=0.0001;
		emin[2]=0.0008, emax[2]=0.0008, esteps[2]=0.0001;
		emin[3]=0.0010, emax[3]=0.0010, esteps[3]=0.0001;
		emin[4]=0.0026, emax[4]=0.0026, esteps[4]=0.0001; */
												 
		file = fopen(filename, "a+");
		fprintf(file, "\n========\n w = %.4f k = %.4f \n========", w, k);
		fclose(file);
		printf("\n========\n w = %.4f k = %.4f \n========", w, k);
		while (1) {
			tRound = clock();
			// optimize parameters
			if (TypeFlag == 0) {
				Prange(k,w,MINvalues,&Tmin);
			} else if (TypeFlag == 1) {
				BJMMPlus(k,w,depth,spaceBound,pmin,pmax,psteps,lmin[1],lmax[1],lsteps[1],emin,emax,esteps,MINvalues,MINParams,&Tmin);
			} else if (TypeFlag == 2) {
				NewV3(k,w,depth,spaceBound,pmin,pmax,psteps,lmin,lmax,lsteps,emin,emax,esteps,wwmin,wwmax,wwsteps,MINvalues,MINParams,&Tmin);
			}
			file = fopen(filename, "a+");
			time(&tNow);
			fprintf(file, "\n\n Time taken: %.2fs / %s", (double)(clock() - tRound)/CLOCKS_PER_SEC, asctime(localtime(&tNow)));
			printf("\n\n Time taken: %.2fs / %s", (double)(clock() - tRound)/CLOCKS_PER_SEC, asctime(localtime(&tNow)));
			// print parameters
			if (LPNFlag == 1) Tmin *= LPN_k / k;
			fprintf(file, "T:%f", Tmin);
			printf("T:%f", Tmin);
			fclose(file);
			if (TypeFlag == 0) { // break if Prange
				break;
			}
			file = fopen(filename, "a+");
			fprintf(file, " p[0]=%.4f ", MINParams[0]);
			printf(" p[0]=%.4f ", MINParams[0]);
			if (TypeFlag == 1) {
				fprintf(file, "l=%.4f ", MINParams[depth]);
				printf("l=%.4f ", MINParams[depth]);
				for(int i=1; i<depth; i++) fprintf(file, "e[%d]=%.5f ", i, MINParams[2*depth-2+i]);
				for(int i=1; i<depth; i++) printf("e[%d]=%.5f ", i, MINParams[2*depth-2+i]);
		
			} else {
				for(int i=1; i<depth; i++) fprintf(file, "w[%d]=%.4f ", i, MINParams[i]);
				for(int i=1; i<depth; i++) fprintf(file, "l[%d]=%.4f ", i, MINParams[depth-1+i]);
				for(int i=1; i<depth; i++) fprintf(file, "e[%d]=%.5f ", i, MINParams[2*depth-2+i]);			
				for(int i=1; i<depth; i++) printf("w[%d]=%.4f ", i, MINParams[i]);
				for(int i=1; i<depth; i++) printf("l[%d]=%.4f ", i, MINParams[depth-1+i]);
				for(int i=1; i<depth; i++) printf("e[%d]=%.5f ", i, MINParams[2*depth-2+i]);			
			}
			// enable to try just one parameter set
			if (SingleFlag == 1) break;
			// adjust parameter ranges
			contFlag = 0;
			if (((int) round (pmin/psteps) != fmax(0,(int) round (MINParams[0]/psteps)-3)) || ((int) round (pmax/psteps) != (int) round (MINParams[0]/psteps)+3)) {
				contFlag = 1;
				fprintf(file, "p");
				printf("p");
			}
			for(int i=1; i<depth; i++) {
				if (MMT_Flag == 0 && (((int) round (emin[i]/esteps[i]) != fmax(0,(int) round (MINParams[2*depth-2+i]/esteps[i])-3)) || ((int) round (emax[i]/esteps[i]) != (int) round (emin[i]/esteps[i])+6))) {
					contFlag = 1;
					fprintf(file, "e");
					printf("e");
				}
				if (TypeFlag == 1 && i > 1) continue;
				if (((int) round (lmin[i]/lsteps[i]) != fmax(0,(int) round (MINParams[depth-1+i]/lsteps[i])-3)) || ((int) round (lmax[i]/lsteps[i]) != (int) round (lmin[i]/lsteps[i])+6)) {
					contFlag = 1;
					fprintf(file, "l");
					printf("l");
				}
				if (TypeFlag == 1) continue;
				if (((int) round (wwmin[i][0]/wwsteps[i][0]) != fmax(0,(int) round (MINParams[i]/wwsteps[i][0])-3)) || ((int) round (wwmax[i][0]/wwsteps[i][0]) != (int) round (wwmin[i][0]/wwsteps[i][0])+6)) {
					contFlag = 1;
					fprintf(file, "w");
					printf("w");
				}
			}
			
			if (!contFlag) {
				if (psteps >= 0.01) {
					psteps *= 0.1;
					contFlag = 1;
					fprintf(file, "p+");
					printf("p+");
				} else if (TypeFlag == 2 && wwsteps[1][0] >= 0.01) {
					for(int i=1; i<depth; i++) {
						wwsteps[i][0] *= 0.1;
					}
					contFlag = 1;
					fprintf(file, "w+");
					printf("w+");
				} else if (lsteps[1] >= 0.01) {
					for(int i=1; i<depth; i++) {
						lsteps[i] *= 0.1;
						if (TypeFlag == 1) break;
					}
					contFlag = 1;
					fprintf(file, "l+");
					printf("l+");
				} else if (MMT_Flag == 0 && esteps[1] >= 0.001) {
					for(int i=1; i<depth; i++) {
						esteps[i] *= 0.1;
					}
					contFlag = 1;
					fprintf(file, "e+");
					printf("e+");
				} else if (PreciseFlag == 1) {
					if (psteps >= 0.001) {
						psteps *= 0.1;
						contFlag = 1;
						fprintf(file, "p+");
						printf("p+");
					} else if (TypeFlag == 2 && wwsteps[1][0] >= 0.001) {
						for(int i=1; i<depth; i++) {
							wwsteps[i][0] *= 0.1;
						}
						contFlag = 1;
						fprintf(file, "w+");
						printf("w+");
					} else if (lsteps[1] >= 0.001) {
						for(int i=1; i<depth; i++) {
							lsteps[i] *= 0.1;
							if (TypeFlag == 1) break;
						}
						contFlag = 1;
						fprintf(file, "l+");
						printf("l+");
					} else if (MMT_Flag == 0 && esteps[1] >= 0.0001) {
						for(int i=1; i<depth; i++) {
							esteps[i] *= 0.1;
						}
						contFlag = 1;
						fprintf(file, "e+");
						printf("e+");
					}
				}
			}
			fclose(file);
			if (contFlag) {
				pmin = fmax(0,MINParams[0]-psteps*3);
				pmax = pmin+psteps*6;
				for(int i=1; i<depth; i++) {
					if (MMT_Flag == 0) {
						emin[i] = fmax(0,MINParams[2*depth-2+i]-esteps[i]*3);
						emax[i] = emin[i]+esteps[i]*6;
					}
					if (TypeFlag == 1 && i > 1) continue;
					lmin[i] = fmax(0,MINParams[depth+i-1]-lsteps[i]*3);
					lmax[i] = lmin[i]+lsteps[i]*6;
					if (TypeFlag == 1) continue;
					wwmin[i][0] = fmax(0,MINParams[i]-wwsteps[i][0]*3);
					wwmax[i][0] = wwmin[i][0]+wwsteps[i][0]*6;
					
				}
				continue;
			}
			break;
		}
		
		// save T
		if (((WorstCaseFlag == 1) && (Tmax < Tmin)) || ((WorstCaseFlag == 0) && (Tmax > Tmin))) {
			Tmax=Tmin;
			for (int i=0; i<MaxDepth*8-1; i++) MAXvalues[i]=MINvalues[i];
			for (int i=0; i<MaxDepth*3-3; i++) MAXParams[i]=MINParams[i];
		} else;// break;
	}
	}
	
    // Output
	file = fopen(filename, "a+");
	if (TypeFlag == 0) {
		fprintf(file, "\n\nk=%f w=%f \n", MAXvalues[0], MAXvalues[1]);
		fprintf(file, "T=%.5f \n", Tmax);
		printf("\n\nk=%f w=%f \n", MAXvalues[0], MAXvalues[1]);
		printf("T=%.5f \n", Tmax);
		
	} else if (TypeFlag == 1) {
		fprintf(file, "\n\n=======\nDepth %d: \n=======\n", depth);
		fprintf(file, " k=%f w=%f l=%.4f ", MAXvalues[0], MAXvalues[1], MAXvalues[2]);
		for(int i=1; i<depth; i++) fprintf(file, "R[%d]=%.4f ", i, MAXvalues[i+2]);
		for(int i=0; i<=depth; i++) fprintf(file, "p[%d]=%.6f ", i, MAXvalues[depth+2+i]);
		for(int i=1; i<=depth; i++) fprintf(file, "S[%d]=%.5f ", i, MAXvalues[2*depth+2+i]);
		for(int i=1; i<=depth; i++) fprintf(file, "C[%d]=%.5f ", i, MAXvalues[3*depth+2+i]);
		if (LPNFlag==1) {
			fprintf(file, "\nT=%.5f/%.5f (P=%.5f) \n", Tmax, MAXvalues[4*depth+3], MAXvalues[4*depth+4]);
		} else {
			fprintf(file, "\nT=%.5f (P=%.5f) \n", MAXvalues[4*depth+3], MAXvalues[4*depth+4]);
		}
		printf("\n\n=======\nDepth %d: \n=======\n", depth);
		printf(" k=%f w=%f l=%.4f ", MAXvalues[0], MAXvalues[1], MAXvalues[2]);
		for(int i=1; i<depth; i++) printf("R[%d]=%.4f ", i, MAXvalues[i+2]);
		for(int i=0; i<=depth; i++) printf("p[%d]=%.6f ", i, MAXvalues[depth+2+i]);
		for(int i=1; i<=depth; i++) printf("S[%d]=%.5f ", i, MAXvalues[2*depth+2+i]);
		for(int i=1; i<=depth; i++) printf("C[%d]=%.5f ", i, MAXvalues[3*depth+2+i]);
		if (LPNFlag==1) {
			printf("\nT=%.5f/%.5f (P=%.5f) \n", Tmax, MAXvalues[4*depth+3], MAXvalues[4*depth+4]);
		} else {
			printf("\nT=%.5f (P=%.5f) \n", MAXvalues[4*depth+3], MAXvalues[4*depth+4]);
		}
		
		fprintf(file, "\n p[0]=%.4f ", MAXParams[0]);
		fprintf(file, "l=%.4f ", MAXParams[depth]);
		for(int i=1; i<depth; i++) fprintf(file, "e[%d]=%.5f ", i, MAXParams[2*depth-2+i]);
		printf("\n p[0]=%.4f ", MAXParams[0]);
		printf("l=%.4f ", MAXParams[depth]);
		for(int i=1; i<depth; i++) printf("e[%d]=%.5f ", i, MAXParams[2*depth-2+i]);
		
	} else if (TypeFlag == 2) {
		fprintf(file, "\n\n=======\nDepth %d: \n=======\n", depth);
		fprintf(file, "k=%f w=%f \n", MAXvalues[0], MAXvalues[1]);
		for(int i=1; i<depth; i++) fprintf(file, "R[%d]=%.4f ", i, MAXvalues[2*depth+i-1]);
		fprintf(file, "\n");
		for(int i=1; i<depth; i++) fprintf(file, "w[%d][%d]=%.6f ", i, i, MAXvalues[3*depth+i-2]);
		for(int i=2; i<depth; i++) fprintf(file, "w[%d][%d]=%.6f ", i, i-1, MAXvalues[4*depth+i-4]);
		fprintf(file, "\n");
		for(int i=1; i<=depth; i++) fprintf(file, "p[%d]=%.6f ", i, MAXvalues[5*depth+i-4]);
		fprintf(file, "\n");
		for(int i=1; i<=depth; i++) fprintf(file, "S[%d]=%.5f ", i, MAXvalues[6*depth+i-4]);
		for(int i=1; i<=depth; i++) fprintf(file, "C[%d]=%.5f ", i, MAXvalues[7*depth+i-4]);
		fprintf(file, "\n");
		if (LPNFlag==1) {
			fprintf(file, "T=%.5f/%.5f (P=%.5f) \n", Tmax, MAXvalues[8*depth-3], MAXvalues[8*depth-2]);
		} else {
			fprintf(file, "T=%.5f (P=%.5f) \n", MAXvalues[8*depth-3], MAXvalues[8*depth-2]);
		}
		printf("\n\n=======\nDepth %d: \n=======\n", depth);
		printf("k=%f w=%f \n", MAXvalues[0], MAXvalues[1]);
		for(int i=1; i<depth; i++) printf("R[%d]=%.4f ", i, MAXvalues[2*depth+i-1]);
		printf("\n");
		for(int i=1; i<depth; i++) printf("w[%d][%d]=%.6f ", i, i, MAXvalues[3*depth+i-2]);
		for(int i=2; i<depth; i++) printf("w[%d][%d]=%.6f ", i, i-1, MAXvalues[4*depth+i-4]);
		printf("\n");
		for(int i=1; i<=depth; i++) printf("p[%d]=%.6f ", i, MAXvalues[5*depth+i-4]);
		printf("\n");
		for(int i=1; i<=depth; i++) printf("S[%d]=%.5f ", i, MAXvalues[6*depth+i-4]);
		for(int i=1; i<=depth; i++) printf("C[%d]=%.5f ", i, MAXvalues[7*depth+i-4]);
		printf("\n");
		if (LPNFlag==1) {
			printf("T=%.5f/%.5f (P=%.5f) \n", Tmax, MAXvalues[8*depth-3], MAXvalues[8*depth-2]);
		} else {
			printf("T=%.5f (P=%.5f) \n", MAXvalues[8*depth-3], MAXvalues[8*depth-2]);
		}

		fprintf(file, "\n p[0]=%.4f ", MAXParams[0]);
		for(int i=1; i<depth; i++) fprintf(file, "w[%d]=%.4f ", i, MAXParams[i]);
		for(int i=1; i<depth; i++) fprintf(file, "l[%d]=%.4f ", i, MAXParams[depth+i-1]);
		for(int i=1; i<depth; i++) fprintf(file, "e[%d]=%.5f ", i, MAXParams[2*depth-2+i]);
		printf("\n p[0]=%.4f ", MAXParams[0]);
		for(int i=1; i<depth; i++) printf("w[%d]=%.4f ", i, MAXParams[i]);
		for(int i=1; i<depth; i++) printf("l[%d]=%.4f ", i, MAXParams[depth+i-1]);
		for(int i=1; i<depth; i++) printf("e[%d]=%.5f ", i, MAXParams[2*depth-2+i]);
    }
	
    // Output time
	fprintf(file, "\n Time taken: %.2fs\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    printf("\n Time taken: %.2fs\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	fclose(file);
    return 0;

}