#define scale 1000000


// Binary Entropy Function
double H_(double x) {
    return (-x*log(x)-(1-x)*log(1-x))/log(2);
}

// Binary Entropy Function, Precomputed
double H(double H1[], double x) {
	if (x!=x) {
		GLOBAL_STOP=1;
		return 0.0;
	}
	if (x<0.000001 && x > -0.000001) x = 0.0;
	if (x>0.999999 && x < 1.000001) x = 1.0;
    if (x<0.0 || x>1.0) {
        //printf("-H Error %f-\n", x);
		GLOBAL_STOP=1;
        return 0.0;
    }
    if (x==0.0 || x==1.0) return 0;
	if (x>0.5) x=1.0-x;
    return H1[(int)round(x*2*scale)];//H_(x);
}

// Inverse Binary Entropy Function, Precomputed
double inverse(double H1[], double y) {
    int left=0, right=scale, mid;
    while(left < right-1) {
        mid= (left+right)/2;
        //printf("%d %d %d \n",left, right, mid);
        if(y < H1[mid]) right=mid-1;
        else left=mid;
    }
    
    //printf("%f \n",fabs(y-H(left/(2*(double)scale))));
    if(fabs(y-H(H1,left/(2*(double)scale))) < fabs(y-H(H1,right/(2*(double)scale))) ) return(left/(2*(double)scale));
    else return (right/(2*(double)scale));
}

// Nearest Neighbor
double NN(double gamma, double lambda, double H1[], int Naive_flag, double *space) {
	double res, AndreNN, AndreOPT;
	if (lambda<0.000001 && lambda > -0.000001) lambda = 0.0;
	if (gamma<0.000001 && gamma > -0.000001) gamma = 0.0;
	*space = fmax(*space, lambda);
	
    if (gamma<0.0 || gamma>=0.5 || lambda<0.0) {
        //printf("-NN Error %f %f-\n", gamma, lambda);
		GLOBAL_STOP=1;
        return 2*lambda;
    }
    if (gamma==0.0) return lambda;
	
	// May Ozerov
    if (!Naive_flag && lambda < 1-H(H1,gamma/2)) {
	   res = (1-gamma) * (1 - H(H1,(inverse(H1, (1-lambda) ) - gamma/2 ) / (1 - gamma) ) );
	   if (res!=res) {
			//printf("-NN Error 2 %f %f -", gamma, lambda);
			GLOBAL_STOP=1;
			return 2*lambda;
	   }
	   
	// Enumerate Pairs
    } else {
		res = 2*lambda;
	}
	//AndreOPT = fmin(lambda, 1-2*gamma);
	//AndreNN = fmax(lambda, 2*lambda-AndreOPT)+H(H1,AndreOPT)-(1-gamma)*H(H1,AndreOPT/(1-gamma));
	//return fmin(fmin(AndreNN,res), fmin(2*lambda, fmax(lambda+H(H1,gamma), 2*lambda-1+H(H1,gamma))));
	
	if (res <= fmax(lambda+H(H1,gamma/2), 2*lambda-1+H(H1,gamma))) {
		return res;
		
	// Meet-in-the-Middle
	} else {
		*space = fmax(*space, lambda+H(H1,gamma/2));
		return fmax(lambda+H(H1,gamma/2), 2*lambda-1+H(H1,gamma));
	}
}