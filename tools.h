#define Scale 1000000
const double ScaleDiff=0.000001;
double H1[Scale+1];

// Binary Entropy Function
double H_(double x) {
    return (-x*log(x)-(1-x)*log(1-x))/log(2);
}

// Precompute H, H^-1
void initH() {
    for (int i = 0; i <= Scale; i++) H1[i] = H_(i / (2 * (double) Scale));
}

// Binary Entropy Function, Precomputed
double H(double x) {
	if (x != x) {
		ErrorFlag = 1;
		return 0.0;
	}
	if (x < ScaleDiff && x > -ScaleDiff) x = 0.0;
	if (x > 1-ScaleDiff && x < 1+ScaleDiff) x = 1.0;
    if (x < 0.0 || x > 1.0) {
		ErrorFlag = 1;
        return 0.0;
    }
    if (x == 0.0 || x == 1.0) return 0;
	if (x > 0.5) x = 1.0 - x;
    return H1[(int) round(x * 2 * Scale)];
}

// Inverse Binary Entropy Function, Precomputed
double inverse(double y) {
    int left = 0, right = Scale, mid;
	if (y != y) {
		ErrorFlag = 1;
		return 0.0;
	}
    while (left < right-1) {
        mid = (left + right) / 2;
        if (y < H1[mid]) right = mid-1;
        else left = mid;
    }
    
    if (fabs(y - H(left / (2 * (double) Scale))) < fabs(y - H(right / (2 * (double) Scale))) ) return (left / (2 * (double)Scale));
    else return (right / (2 * (double) Scale));
}

// Nearest Neighbor
double NN(double gamma, double lambda, int Naive_flag, double *space) {
	double ep, mitm, mo;
	if (gamma != gamma || lambda != lambda) {
		ErrorFlag = 1;
		return 0.0;
	}
	if (lambda < ScaleDiff && lambda > -ScaleDiff) lambda = 0.0;
	if (gamma < ScaleDiff && gamma > -ScaleDiff) gamma = 0.0;
	*space = fmax(*space, lambda);
	
    if (gamma < 0.0 || gamma >= 0.5 || lambda < 0.0) {
		ErrorFlag = 1;
        return 2 * lambda;
    }
    if (gamma == 0.0) return lambda;
	
	// Enumerate Pairs
	ep = 2 * lambda;
	
	// Meet-in-the-Middle
	mitm = fmax(lambda + H( gamma/2), 2 * lambda-1 + H( gamma));
	
	// May Ozerov
    if (!Naive_flag && lambda < 1 - H( gamma / 2)) {
	   mo = (1 - gamma) * (1 - H((inverse( (1 - lambda)) - gamma / 2) / (1 - gamma)));
	   if (mo != mo) {
			ErrorFlag = 1;
			return 2 * lambda;
	   }
	   if (mo < fmin(ep, mitm)) return mo;
	}
	
	if (ep < mitm) return ep;
	
	*space = fmax(*space, lambda + H( gamma/2));
	return mitm;
}