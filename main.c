#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int GLOBAL_STOP=0;
#include "tools.h"
#include "decode.h"
#include <time.h>

int main(int argc, char *argv[]) {
    // running time counter
    clock_t tStart = clock(), tRound;
	srand(time(NULL));
	FILE *file;
	file = fopen("output.txt", "wb");
    
    // Precompute H, H^-1
    static double H1[scale+1];
    for (int i=0; i<=scale; i++) {
        H1[i]=H_(i/(2*(double)scale));
        //printf("%d %f %.20f \n",i, i/(2*(double)scale), H1[i]);
    }
    
    // initialization
    int depth = 0, cont_flag = 0, HD_flag = 0, Naive_flag = 0, Type_flag = 0, WorstCase_flag = 1, NN_flag = 1, LPN_flag = 0, Precise_flag = 1;
    double pmin, pmax, psteps;
    double kmin, kmax, ksteps;
	double wmin=0, wmax=0, wsteps=0.1;
    double emin[6], emax[6], esteps[6];
    for (int i=0;i<6;i++) {emin[i] = 0.0, emax[i] = 0.0, esteps[i] = 1.0;}
    double dwmin[6][6], dwmax[6][6], dwsteps[6][6];
    for (int i=0;i<6;i++) {
		for (int j=0;j<6;j++) {
			dwmin[i][j] = 0.0, dwmax[i][j] = 0.0, dwsteps[i][j] = 1.0;
		}
	}
    double wwmin[6], wwmax[6], wwsteps[6];
    for (int i=0;i<6;i++) {wwmin[i] = 0.0, wwmax[i] = 0.0, wwsteps[i] = 1.0;}
    double lmin[6], lmax[6], lsteps[6];
    for (int i=0;i<6;i++) {lmin[i] = 0.0, lmax[i] = 0.0, lsteps[i] = 1.0;}
    double k=0, w=0, LPN_k, LPN_t;
    double Tmax=0.0, Tmin=1.0, MAXvalues[32], MINvalues[32], MINParams[20], MAXParams[20];
	
	// standard settings
	Precise_flag = 1;
	WorstCase_flag = 1;
	HD_flag = 0;
	Naive_flag = 0;
	NN_flag = 1; // only BJMM+
	Type_flag = 2; // 0: Prange, 1: BJMM+, 2: BJMM v3
	depth=2;
	LPN_flag = 0;
	LPN_k = 512;
	LPN_t = 0.25;
	ksteps=0.05, kmin=0.40, kmax=0.50;
	wsteps=0.01	, wmin=0.00, wmax=0.00;
	
	// read input
	int kset=0, wset=0;
	char *item[2];
    for (int i=1; i<argc; i++) {
		item[0] = strtok(argv[i], "=");
		item[1] = strtok(NULL, "=");
		if (strcmp(item[0], "m")==0) {
            depth = atoi(item[1]);
        }
		if (strcmp(item[0],"kmin")==0) {
            kmin = atof(item[1]);
			kset=1;
        }
		if (strcmp(item[0],"kmax")==0) {
            kmax = atof(item[1]);
			kset=1;
        }
		if (strcmp(item[0],"ksteps")==0) {
            ksteps = atof(item[1]);
			kset=1;
        }
		
		if (strcmp(item[0], "Precise")==0) {
            Precise_flag = atoi(item[1]);
        }
        if (strcmp(item[0], "WC")==0) {
            WorstCase_flag = atoi(item[1]);
        }
		if (strcmp(item[0], "HD")==0) {
            HD_flag = atoi(item[1]);
        }
		if (strcmp(item[0], "Naive")==0) {
            Naive_flag = atoi(item[1]);
        }
		if (strcmp(item[0], "NN")==0) {
            NN_flag = atoi(item[1]);
        }
		if (strcmp(item[0], "Type")==0) {
            Type_flag = atoi(item[1]);
        }
		if (strcmp(item[0], "LPN")==0) {
            LPN_flag = atoi(item[1]);
        }
		if (strcmp(item[0], "LPNk")==0) {
            LPN_k = atoi(item[1]);
        }
        if (strcmp(item[0], "LPNtau")==0) {
            LPN_t = atof(item[1]);
        }
		if (strcmp(item[0],"wmin")==0) {
            wmin = atof(item[1]);
			wset=1;
        }
		if (strcmp(item[0],"wmax")==0) {
            wmax = atof(item[1]);
			wset=1;
        }
		if (strcmp(item[0], "wsteps")==0) {
            wsteps = atof(item[1]);
			wset=1;
        }
    }
	
	//standard settings part 2
	if (WorstCase_flag == 0) Tmax=1.0;
	if (LPN_flag == 1) {
		Tmax *= LPN_k / kmin;
		if (wset==0) wmin=wmax=LPN_t;
	}
	
	// Optimization for every k, w
	for (w=wmin; w<=wmax+wsteps*0.1; w+=wsteps) {
	for (k=kmin; k<=kmax+ksteps*0.1; k+=ksteps) {
		
		// Our algorithm
		if (Type_flag==2) {
			if (depth==2) {
				pmin=0.00, pmax=0.06, psteps=0.01;
				wwmin[1]=0.01, wwmax[1]=0.07, wwsteps[1]=0.01;
				lmin[1]=0.13, lmax[1]=0.19, lsteps[1]=0.01;
				emin[1]=0.004, emax[1]=0.010, esteps[1]=0.001;
			} else if (depth==3) {
				pmin=0.02, pmax=0.06, psteps=0.01;
				wwmin[1]=0.01, wwmax[1]=0.05, wwsteps[1]=0.01;
				wwmin[2]=0.01, wwmax[2]=0.05, wwsteps[2]=0.01;
				lmin[1]=0.11, lmax[1]=0.15, lsteps[1]=0.01;
				lmin[2]=0.03, lmax[2]=0.07, lsteps[2]=0.01;
				emin[1]=0.001, emax[1]=0.005, esteps[1]=0.001;
				emin[2]=0.001, emax[2]=0.005, esteps[2]=0.001;
			} else if (depth==4) {
				pmin=0.00, pmax=0.06, psteps=0.01;
				wwmin[1]=0.00, wwmax[1]=0.06, wwsteps[1]=0.01;
				wwmin[2]=0.00, wwmax[2]=0.06, wwsteps[2]=0.01;
				wwmin[3]=0.00, wwmax[3]=0.06, wwsteps[3]=0.01;
				dwmin[3][2]=0.00, dwmax[3][2]=0.06, dwsteps[3][2]=0.01;
				lmin[1]=0.15, lmax[1]=0.21, lsteps[1]=0.01;
				lmin[2]=0.05, lmax[2]=0.11, lsteps[2]=0.01;
				lmin[3]=0.00, lmax[3]=0.06, lsteps[3]=0.01;
				emin[1]=0.000, emax[1]=0.006, esteps[1]=0.001;
				emin[2]=0.000, emax[2]=0.006, esteps[2]=0.001;
				emin[3]=0.000, emax[3]=0.006, esteps[3]=0.001;
			}
			
		// BJMM(d)
		} else if (Type_flag==1) {
			if (depth==2) {
				pmin=0.00, pmax=0.07, psteps=0.01;
				lmin[1]=0.00, lmax[1]=0.12, lsteps[1]=0.01;
				emin[1]=0.000, emax[1]=0.010, esteps[1]=0.001;	   
			} else if (depth==3) {
				pmin=0.00, pmax=0.09, psteps=0.01;
				lmin[1]=0.00, lmax[1]=0.30, lsteps[1]=0.01;
				emin[1]=0.000, emax[1]=0.025, esteps[1]=0.001;
				emin[2]=0.000, emax[2]=0.007, esteps[2]=0.001;
			} else if (depth==4) {
				pmin=0.00, pmax=0.10, psteps=0.01;
				lmin[1]=0.00, lmax[1]=0.28, lsteps[1]=0.01;
				emin[1]=0.000, emax[1]=0.035, esteps[1]=0.001;
				emin[2]=0.000, emax[2]=0.018, esteps[2]=0.001;
				emin[3]=0.000, emax[3]=0.006, esteps[3]=0.001;
			}
		}		
												 
		printf("\n========\n w = %.4f k = %.4f \n========\n", w, k);
		while (1) {
			tRound = clock();
			
			// optimize parameters
			if (Type_flag == 0) {
				Prange(k,w,HD_flag,H1,MINvalues,&Tmin);
			} else if (Type_flag == 1) {
				BJMMPlus(k,w,depth,HD_flag,NN_flag,Naive_flag,pmin,pmax,psteps,lmin[1],lmax[1],lsteps[1],emin,emax,esteps,H1,MINvalues,MINParams,&Tmin);
			} else if (Type_flag == 2) {
				NewV3(k,w,depth,HD_flag,Naive_flag,pmin,pmax,psteps,lmin,lmax,lsteps,emin,emax,esteps,dwmin,dwmax,dwsteps,wwmin,wwmax,wwsteps,H1,MINvalues,MINParams,&Tmin);
			}
			printf("\nTime taken: %.2fs", (double)(clock() - tRound)/CLOCKS_PER_SEC);
			
			// print parameters
			if (LPN_flag == 1) Tmin *= LPN_k / k;
			printf("\nT:%f", Tmin);
			if (Type_flag == 0) { // break if Prange
				break;
			}
			printf(" p[0]=%.4f ", MINParams[0]);
			if (Type_flag == 1) {
				printf("l=%.4f ", MINParams[4]);
				for(int i=1; i<depth; i++) printf("e[%d]=%.5f ", i, MINParams[6+i]);
		
			} else {
				for(int i=1; i<depth; i++) printf("w[%d]=%.4f ", i, MINParams[i]);
				for(int i=1; i<depth; i++) printf("l[%d]=%.4f ", i, MINParams[3+i]);
				for(int i=1; i<depth; i++) printf("e[%d]=%.5f ", i, MINParams[6+i]);			
				if (depth==4) printf("dw[3][2]=%.4f ", MINParams[10]);
			}

			// adjust parameter ranges
			cont_flag = 0;
			if (((int) round (pmin/psteps) != fmax(0,(int) round (MINParams[0]/psteps)-3)) || ((int) round (pmax/psteps) != (int) round (pmin/psteps)+6)) {
				cont_flag = 1;
				printf("p");
			}
			for(int i=1; i<depth; i++) {
				if (((int) round (emin[i]/esteps[i]) != fmax(0,(int) round (MINParams[6+i]/esteps[i])-3)) || ((int) round (emax[i]/esteps[i]) != (int) round (emin[i]/esteps[i])+6)) {
					cont_flag = 1;
					printf("e");
				}
				if (Type_flag == 1 && i > 1) continue;
				if (((int) round (lmin[i]/lsteps[i]) != fmax(0,(int) round (MINParams[3+i]/lsteps[i])-3)) || ((int) round (lmax[i]/lsteps[i]) != (int) round (lmin[i]/lsteps[i])+6)) {
					cont_flag = 1;
					printf("l");
				}
				if (Type_flag == 1) continue;
				if (((int) round (wwmin[i]/wwsteps[i]) != fmax(0,(int) round (MINParams[i]/wwsteps[i])-3)) || ((int) round (wwmax[i]/wwsteps[i]) != (int) round (wwmin[i]/wwsteps[i])+6)) {
					cont_flag = 1;
					printf("w");
				}
			}
			if (Type_flag == 2 && depth==4) {
				if (((int) round (dwmin[3][2]/dwsteps[3][2]) != fmax(0,(int) round (MINParams[10]/dwsteps[3][2])-3)) || ((int) round (dwmax[3][2]/dwsteps[3][2]) != (int) round (dwmin[3][2]/dwsteps[3][2])+6)) {
					cont_flag = 1;
					printf("d");
				}
			}
			
			if (!cont_flag) {
				if (psteps >= 0.01) {
					psteps *= 0.1;
					cont_flag = 1;
					printf("p+\n");
				} else if (Type_flag == 2 && wwsteps[1] >= 0.01) {
					for(int i=1; i<depth; i++) {
						wwsteps[i] *= 0.1;
					}
					cont_flag = 1;
					printf("w+\n");
				} else if (lsteps[1] >= 0.01) {
					for(int i=1; i<depth; i++) {
						lsteps[i] *= 0.1;
						if (Type_flag == 1) break;
					}
					cont_flag = 1;
					printf("l+\n");
				} else if (esteps[1] >= 0.001) {
					for(int i=1; i<depth; i++) {
						esteps[i] *= 0.1;
					}
					cont_flag = 1;
					printf("e+\n");
				} else if (Type_flag == 2 && depth==4 && dwsteps[3][2] >= 0.01) {
					dwsteps[3][2] *= 0.1;
					cont_flag = 1;
					printf("d+\n");
				} else if (Precise_flag == 1) {
					if (psteps >= 0.001) {
						psteps *= 0.1;
						cont_flag = 1;
						printf("p+\n");
					} else if (Type_flag == 2 && wwsteps[1] >= 0.001) {
						for(int i=1; i<depth; i++) {
							wwsteps[i] *= 0.1;
						}
						cont_flag = 1;
						printf("w+\n");
					} else if (lsteps[1] >= 0.001) {
						for(int i=1; i<depth; i++) {
							lsteps[i] *= 0.1;
							if (Type_flag == 1) break;
						}
						cont_flag = 1;
						printf("l+\n");
					} else if (esteps[1] >= 0.0001) {
						for(int i=1; i<depth; i++) {
							esteps[i] *= 0.1;
						}
						cont_flag = 1;
						printf("e+\n");
					} else if (Type_flag == 2 && depth==4 && dwsteps[3][2] >= 0.001) {
						dwsteps[3][2] *= 0.1;
						cont_flag = 1;
						printf("d+\n");
					} 
				}
			}
			
			if (cont_flag) {
				pmin = fmax(0,MINParams[0]-psteps*3);
				pmax = pmin+psteps*6;
				for(int i=1; i<depth; i++) {
					emin[i] = fmax(0,MINParams[6+i]-esteps[i]*3);
					emax[i] = emin[i]+esteps[i]*6;
					if (Type_flag == 1 && i > 1) continue;
					lmin[i] = fmax(0,MINParams[3+i]-lsteps[i]*3);
					lmax[i] = lmin[i]+lsteps[i]*6;
					if (Type_flag == 1) continue;
					wwmin[i] = fmax(0,MINParams[i]-wwsteps[i]*3);
					wwmax[i] = wwmin[i]+wwsteps[i]*6;
					
				}
				if (Type_flag == 2 && depth==4) {
					dwmin[3][2] = fmax(0,MINParams[10]-dwsteps[3][2]*3);
					dwmax[3][2] = dwmin[3][2]+dwsteps[3][2]*6;
				}
				continue;
			}
			break;
		}
		
		// save T
		if (((WorstCase_flag == 1) && (Tmax < Tmin)) || ((WorstCase_flag == 0) && (Tmax > Tmin))) {
			Tmax=Tmin;
			for (int i=0; i<32; i++) MAXvalues[i]=MINvalues[i];
			for (int i=0; i<20; i++) MAXParams[i]=MINParams[i];
		} else;// break;
		if (file != NULL) {
			char x[23];
			sprintf(x, "%.2f;%.5f;%3.5f\n", w, k, fabs(Tmin));
			fwrite(x, sizeof(x), 1, file);
		}
	}
	}
	
    // Output
	if (Type_flag == 0) {
		printf("\n\nk=%f w=%f \n", MAXvalues[0], MAXvalues[1]);
		printf("T=%.5f \n", Tmax);
		
	} else if (Type_flag == 1) {
		printf("\n=======\nTiefe %d: \n=======\n", depth);
		printf(" k=%f w=%f l=%.4f ", MAXvalues[0], MAXvalues[1], MAXvalues[2]);
		for(int i=1; i<depth; i++) printf("R[%d]=%.4f ", i, MAXvalues[i+2]);
		for(int i=0; i<=depth; i++) printf("p[%d]=%.6f ", i, MAXvalues[depth+2+i]);
		for(int i=1; i<=depth; i++) printf("S[%d]=%.5f ", i, MAXvalues[2*depth+2+i]);
		for(int i=1; i<=depth; i++) printf("C[%d]=%.5f ", i, MAXvalues[3*depth+2+i]);
		if (LPN_flag==1) {
			printf("\nT=%.5f/%.5f (P=%.5f) \n", Tmax, MAXvalues[4*depth+3], MAXvalues[4*depth+4]);
		} else {
			printf("\nT=%.5f (P=%.5f) \n", MAXvalues[4*depth+3], MAXvalues[4*depth+4]);
		}
		
		printf("\n p[0]=%.4f ", MAXParams[0]);
		printf("l=%.4f ", MAXParams[4]);
		for(int i=1; i<depth; i++) printf("e[%d]=%.5f ", i, MAXParams[6+i]);
		
	} else if (Type_flag == 2) {
		printf("\n=======\nTiefe %d: \n=======\n", depth);
		printf("k=%f w=%f \n", MAXvalues[0], MAXvalues[1]);
		for(int i=1; i<depth; i++) printf("R[%d]=%.4f ", i, MAXvalues[2*depth+i-1]);
		printf("\n");
		for(int i=1; i<depth; i++) printf("dw[%d]=%.6f ", i, MAXvalues[3*depth+i-2]);
		for(int i=2; i<depth; i++) printf("dw[%d][1]=%.6f ", i, MAXvalues[4*depth+i-4]);
		printf("\n");
		for(int i=1; i<=depth; i++) printf("p[%d]=%.6f ", i, MAXvalues[5*depth+i-4]);
		printf("\n");
		for(int i=1; i<=depth; i++) printf("S[%d]=%.5f ", i, MAXvalues[6*depth+i-4]);
		for(int i=1; i<=depth; i++) printf("C[%d]=%.5f ", i, MAXvalues[7*depth+i-4]);
		printf("\n");
		if (LPN_flag==1) {
			printf("T=%.5f/%.5f (P=%.5f) \n", Tmax, MAXvalues[8*depth-3], MAXvalues[8*depth-2]);
		} else {
			printf("T=%.5f (P=%.5f) \n", MAXvalues[8*depth-3], MAXvalues[8*depth-2]);
		}

		printf("\n p[0]=%.4f ", MAXParams[0]);
		for(int i=1; i<depth; i++) printf("w[%d]=%.4f ", i, MAXParams[i]);
		for(int i=1; i<depth; i++) printf("l[%d]=%.4f ", i, MAXParams[3+i]);
		for(int i=1; i<depth; i++) printf("e[%d]=%.5f ", i, MAXParams[6+i]);
		if (depth==4) printf("dw[3][2]=%.4f\n", MAXParams[10]);
    }
	
    // Output time
    printf("\n Time taken: %.2fs\n\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    
	fclose(file);
    return 0;

}