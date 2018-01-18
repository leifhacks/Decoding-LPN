CONTENTS OF THIS FILE
---------------------
   
 * Introduction
 * Configuration
 * Authors
 
 
INTRODUCTION
------------

This program is part of our paper 
"Decoding Linear Codes with High Error Rate and its Impact for LPN Security" 
(https://eprint.iacr.org/2017/1139) . 
Its purpose is to give estimates of the security for instances of the
 - syndrome decoding problem
 - LPN problem

It calculates upper bounds for the running time exponents of 
 - Pranges' Information Set Decoding algorithm [Prange61]
 - the BJMM Algorithm with or without May-Ozerov Nearest Neighbor and depth=2,3,4 [BJMM12, MO15, BM17]
 - our new decoding algorithm ().
 
CONFIGURATION
-------------

The program can be run using the comment 

./main [Params]


WorstCase_flag = 1;
	ksteps=0.05, kmin=0.40, kmax=0.50;
	wsteps=0.01	, wmin=0.00, wmax=0.00;
	
where the following parameters are possible:

Type	// 1=Prange, 2=BJMM, 3=Our algorithm (Standard: 2)
NN		// 1=with Nearest Neighbor, 0=without.. (only BJMM, Standard= 1)
Naive	// 0=with May-Ozerov Nearest Neighbor 1=without..	(only BJMM + Our algorithm, Standard: 0)
m		// depth of the binary search tree (possible values: 2,3,4, Standard: 2)
WC		// 1=Worst Case over all omega,k 0=Best Case over all omega,k (Standard: 1)
Precise // 1=High Precision for paramerters, 2=Low Precision (only BJMM + Our algorithm, Standard: 1)
LPN		// 1=LPN instance, 0= No LPN (Standard: 0)
LPNk	// k for LPN (Standard: 512)
LPNtau	// tau for LPN (Standard: 0.25)
HD		// 1=Half Distance Decoding, 0=Full Distance Decoding (only if LPN=0, Standard: 0)
kmin	// minimum code rate k/n (Standard: 0.4)
kmax	// maximum code rate k/n (Standard: 0.5)
ksteps	// every iteration ksteps is added to the code rate k/n till kmax is reached (Standard: 0.05)
wmin	// minimum error rate omega/n (Standard: 0.0)
wmax	// maximum error rate omega/n (Standard: 0.0)
wsteps	// every iteration wsteps is added to the error rate omega/n till wmax is reached (Standard: 0.01)

Note: Do not set wmin, wmax if you want get running times for Full/Half Distance decoding as omega is set 
via the Gilbert Varshamov bound then.

EXAMPLES
-------------

Worst case running time of our algorithm with depth 3 for Half Distance Decoding over k/n=0.4 ,0.41,..,0.5:
 ./main Type=2 m=3 kmin=0.40 kmax=0.50 ksteps=0.01 HD=1
 
Worst case running time of BJMM without NN and depth 2 for Full Distance Decoding over k/n=0.4 ,0.41,..,0.5: 
 ./main Type=1 NN=0 m=2 kmin=0.40 kmax=0.50 ksteps=0.01
 
Running time of Prange for McEliece istance k/n=0.775, omega/n=0.02:
 ./main Type=0 kmin=0.775 kmax=0.775 wmin=0.02 wmax=0.02
 
Optimal running time for LPN(512,0.25) for our algorithm with depth 3 for k/n=0.0040 ,0.0041,..,0.0060:
 ./main Type=2 m=3 kmin=0.0040 kmax=0.0060 ksteps=0.0001 WC=0 LPN=1 LPNk=512 LPNtau=0.25
 
OUTPUT, Example
-------------
   
k=0.460000 w=0.061867									// k/n and omega/n
R[1]=0.0333												// # of representations on layer 1
w[1][1]=0.001407											// w_1^(1) (only our algorithm)
p[1]=0.008220 p[2]=0.004110								// p_1, p_2
S[1]=0.02865 S[2]=0.02974 C[1]=0.03184 C[2]=0.03184		// S[i]: Size of lists in layer i, C[i]: running 
														// time for step i
T=0.04871/0.04871 (P=0.01687)							// log(T) (only LPN), running time exponent c in 
														// 2^cn (Decoding) / 2^ck (LPN), # Iterations in 
														// outter loop (for good pi)

p[0]=0.0116 w[1]=0.0015 l[1]=0.0396 e[1]=0.00242		// optimized parameters: p_0, w_1, l_1, 
														// e_1 = p_1/2-p_0
Time taken: 0.68s										// time for calculation

Note that the notation for our algorithm in this program slightly differs from the one in our paper. The 
algorithm in our paper starts with lists L_1^{(0)}...L_{2^m}^{(0)} resulting in a final list L^{(m)}. This 
program enumerates the lists vice versa, i.e. starting with lists L_i^{(m)} and a final list L^{(0)}. 
Therefore all variables are enumerated in a different order.

   
AUTHORS
-----------

Leif Both http://cits.rub.de/personen/both.html
Alexander May http://cits.rub.de/personen/may.html