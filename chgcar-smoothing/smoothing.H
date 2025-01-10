#ifndef __SMOOTHING_H
#define __SMOOTHING_H
#include <math.h>
#include <stdio.h>

const double RHO_MIN_SCALE = 1e-4;

double min_nn(double latt[9], int N, double** u);

inline double cube_root (double x) { return exp(log(x)/3.); }

void add_list (double* G_l, double* rho_l, int& Nl, double G, double rho);

void extract_filter_gaussian(double* rhoR, double* rhoI, 
                    int N[3], double latt[9], double& l0);

void apply_filter_gaussian(double* rhoR, double* rhoI, 
			   int N[3], double latt[9],
			   double l0);

void extract_filter(double* rhoR, double* rhoI, 
		    int N[3], double latt[9],
		    double& filter_a, double& filter_b);

void apply_filter(double* rhoR, double* rhoI, 
		  int N[3], double latt[9],
		  double filter_a, double filter_b);

// 3^3 * N*(N+1)/2 algorithm for finding shortest NN distance
double min_nn(double latt[9], int N, double** u) 
{
    double metric[9];
    square(latt, metric);
    double rmin2=metric[0];
    for (int m=0; m<N; ++m)
    for (int n=m; n<N; ++n) {
        double ud0[3] = {u[n][0]-u[m][0],u[n][1]-u[m][1],u[n][2]-u[m][2]}, ud[3];
        int shift[3];
        for (shift[0]=-1; shift[0]<=1; ++shift[0])
            for (shift[1]=-1; shift[1]<=1; ++shift[1])
                for (shift[2]=-1; shift[2]<=1; ++shift[2]) {
	                for (int d=0; d<3; ++d) ud[d] = ud0[d] + shift[d];
                    double rmagn2 = magnsq(metric, ud);
                    if (rmagn2 == 0) continue;
                    if (rmagn2 < rmin2) rmin2 = rmagn2;
                }
    }
    return sqrt(rmin2);
}

void add_list (double* G_l, double* rho_l, int& Nl, double G, double rho) {
    // I don't know why this was here... I've removed it
    //  if (zero(rho)) return;
    if (Nl==0) {
        Nl++; G_l[0] = G; rho_l[0] = rho; return;
    }
    
    int i;
    for (i=0; (i<Nl) && (G>G_l[i]); ++i) ;
    if (dcomp(G_l[i], G)) {
        if (rho_l[i] < rho) rho_l[i] = rho;
    }
    else {
        for (int j=Nl; j>i; --j) {
            G_l[j] = G_l[j-1]; rho_l[j] = rho_l[j-1];
        }
        G_l[i] = G; rho_l[i] = rho; ++Nl;
    }
}

// *************************** GAUSSIAN FILTERS ************************

void extract_filter_gaussian(double* rhoR, double* rhoI, int N[3], double latt[9],
                             double& l0) {
    double Gmetric[9];
    {
        double temp[9];
        careful_inverse(latt, temp);	// temp = [a]^-1
        self_transpose(temp); // Definition is changed here, to avoid a problem.
        mult(temp, 2*M_PI, temp);		// temp = 2pi[a]^-1
        square(temp, Gmetric);		// Gmetric = (2pi)^2 [a]^-T[a]^-1
    }
    
    int n[3], nc[3];
    int Nh[3] = {N[0]/2, N[1]/2, N[2]/2};
    
    int nind=0;
    int h[3], half[3];
    int MEMORY = (Nh[0]+1)*(Nh[1]+1)*(Nh[2]+1); // hardcoded, I know...
    double* G_l = new double[MEMORY];
    double* rho_l = new double[MEMORY];
    int Nl = 0;
    
    double rho_min = rhoR[0] * RHO_MIN_SCALE;

    for (n[2]=0; n[2]<N[2]; ++(n[2]))
        for (n[1]=0; n[1]<N[1]; ++(n[1]))
            for (n[0]=0; n[0]<N[0]; ++(n[0])) {
                for (int d=0; d<3; ++d) {
                    if (n[d]>Nh[d]) nc[d] = n[d]-N[d];
                    else            nc[d] = n[d];
                    if (n[d]==Nh[d]) half[d] = 1;
                    else             half[d] = 0;
                }
                double scaled_rho = hypot(rhoR[nind], rhoI[nind])
                /( (half[0]+1)*(half[1]+1)*(half[2]+1) );
                ++nind;
                if (scaled_rho < rho_min) continue;

                for (h[0]=0; h[0]<=half[0]; ++(h[0]))
                    for (h[1]=0; h[1]<=half[1]; ++(h[1]))
                        for (h[2]=0; h[2]<=half[2]; ++(h[2])) {
                            double nc_h[3];
                            for (int d=0; d<3; ++d) nc_h[d] = nc[d] - h[d]*N[d];
                            double Gabs2 = magnsq(Gmetric, nc_h);
                            add_list(G_l, rho_l, Nl, Gabs2, log(scaled_rho));
                            if (Nl >= MEMORY) {
                                fprintf(stderr, "** BAD MOJO: exceeded memory requirement in extract_filter_gaussian **\n");
                                exit(-1);
                            }
                        }
            }

    // now, fit to a Gaussian (which is a linear fit, in this case)
    double G_sum=0, G2_sum=0, rho_sum=0, rhoG_sum=0;
    for (int i=0; i<Nl; ++i) {
        G_sum    += G_l[i];
        G2_sum   += G_l[i] * G_l[i];
        rho_sum  += rho_l[i];
        rhoG_sum += rho_l[i] * G_l[i];
    }

    // garbage collection
    delete[] G_l; delete[] rho_l;
    
    l0 = (Nl*rhoG_sum - G_sum*rho_sum)/(Nl*G2_sum - G_sum*G_sum);
    if (l0>0) {
        fprintf(stderr, "Somehow did not get a _decreasing_ Gaussian in extract_filter_gaussian\n");
        l0 = 0;
    }
    l0 = sqrt(-l0)*(2.*M_PI);  // scale accordingly.
}

// assumes that l0 has been scaled by some factor
void apply_filter_gaussian(double* rhoR, double* rhoI, 
			   int N[3], double latt[9],
			   double l0) {

    double Gmetric[9];
    {
        double temp[9];
        careful_inverse(latt, temp);	// temp = [a]^-1
        self_transpose(temp); // Definition is changed here, to avoid a problem.
        mult(temp, 2*M_PI, temp);		// temp = 2pi[a]^-1
        square(temp, Gmetric);		// Gmetric = (2pi)^2 [a]^-T[a]^-1
    }
    
    int n[3], nc[3];
    int Nh[3] = {N[0]/2, N[1]/2, N[2]/2};

    int nind=0;
    double lsquare = (l0*l0)/(4.*M_PI*M_PI);

    for (n[2]=0; n[2]<N[2]; ++(n[2]))
        for (n[1]=0; n[1]<N[1]; ++(n[1]))
            for (n[0]=0; n[0]<N[0]; ++(n[0])) {
                for (int d=0; d<3; ++d) {
                    if (n[d]>Nh[d]) nc[d] = n[d]-N[d];
                    else            nc[d] = n[d];
                }
                double Gabs2 = magnsq(Gmetric, nc);
                double scale = exp(-Gabs2*lsquare);
                rhoR[nind] *= scale; rhoI[nind] *= scale; ++nind;
            }
}

//? *************************** POWER LAW TRUNCATION FILTERS ************************

void extract_filter(double* rhoR, double* rhoI, 
		    int N[3], double latt[9],
		    double& filter_a, double& filter_b) 
{
    double Gmetric[9];
    {
        double temp[9];
        careful_inverse(latt, temp);	// temp = [a]^-1
        self_transpose(temp); // Definition is changed here, to avoid a problem.
        mult(temp, 2*M_PI, temp);		// temp = 2pi[a]^-1
        square(temp, Gmetric);		// Gmetric = (2pi)^2 [a]^-T[a]^-1
    }
    
    int n[3], nc[3];
    int Nh[3] = {N[0]/2, N[1]/2, N[2]/2};
    
    int nind=0;
    int h[3], half[3];
    double Gabs;
    int MEMORY = (Nh[0]+1)*(Nh[1]+1)*(Nh[2]+1); // hardcoded, I know...
    double* G_l = new double[MEMORY];
    double* rho_l = new double[MEMORY];
    int Nl = 0;
    
    double rho_min = rhoR[0] * RHO_MIN_SCALE;
    
    for (n[2]=0; n[2]<N[2]; ++(n[2]))
        for (n[1]=0; n[1]<N[1]; ++(n[1]))
            for (n[0]=0; n[0]<N[0]; ++(n[0])) {
                for (int d=0; d<3; ++d) {
                    if (n[d]>Nh[d]) nc[d] = n[d]-N[d];
                    else            nc[d] = n[d];
                    if (n[d]==Nh[d]) half[d] = 1;
                    else             half[d] = 0;
                }
                double scaled_rho = hypot(rhoR[nind], rhoI[nind])
                /( (half[0]+1)*(half[1]+1)*(half[2]+1) );
                ++nind;
                
                if (zero_vect(nc)) continue;
                if (scaled_rho < rho_min) continue;
                
                for (h[0]=0; h[0]<=half[0]; ++(h[0]))
                    for (h[1]=0; h[1]<=half[1]; ++(h[1]))
                        for (h[2]=0; h[2]<=half[2]; ++(h[2])) {
                            double nc_h[3];
                            for (int d=0; d<3; ++d) nc_h[d] = nc[d] - h[d]*N[d];
                            Gabs = sqrt(magnsq(Gmetric, nc_h));
                            add_list(G_l, rho_l, Nl, log(Gabs), log(scaled_rho));
                            if (Nl >= MEMORY) {
                                fprintf(stderr, "** BAD MOJO: exceeded memory requirement in extract_filter **\n");
                                exit(-1);
                            }
                        }
            }
    
    // now, fit to a power law...
    double G_sum=0, G2_sum=0, rho_sum=0, rhoG_sum=0;
    for (int i=0; i<Nl; ++i) {
        G_sum    += G_l[i];
        G2_sum   += G_l[i] * G_l[i];
        rho_sum  += rho_l[i];
        rhoG_sum += rho_l[i] * G_l[i];
    }
    filter_b = (Nl*rhoG_sum - G_sum*rho_sum)/(Nl*G2_sum - G_sum*G_sum);
    
    filter_a = 0;
    for (int i=0; i<Nl; ++i) {
        double guess_a = exp(rho_l[i] - filter_b*G_l[i]);
        if (guess_a > filter_a) filter_a = guess_a;
    }
    
    // garbage collection
    delete[] G_l; delete[] rho_l;
}

void apply_filter(double* rhoR, double* rhoI, 
		  int N[3], double latt[9],
		  double filter_a, double filter_b) 
{
    double Gmetric[9];
    {
        double temp[9];
        careful_inverse(latt, temp);	// temp = [a]^-1
        self_transpose(temp); // Definition is changed here, to avoid a problem.
        mult(temp, 2*M_PI, temp);		// temp = 2pi[a]^-1
        square(temp, Gmetric);		// Gmetric = (2pi)^2 [a]^-T[a]^-1
    }
    
    int n[3], nc[3];
    int Nh[3] = {N[0]/2, N[1]/2, N[2]/2};
    
    int nind=0;
    double Gabs, rho_magn, rho_max;
    
    for (n[2]=0; n[2]<N[2]; ++(n[2]))
        for (n[1]=0; n[1]<N[1]; ++(n[1]))
            for (n[0]=0; n[0]<N[0]; ++(n[0])) {
                for (int d=0; d<3; ++d) {
                    if (n[d]>Nh[d]) nc[d] = n[d]-N[d];
                    else            nc[d] = n[d];
                }
                rho_magn = hypot(rhoR[nind], rhoI[nind]);
                ++nind;
                
                if (zero_vect(nc)) continue;
                Gabs = sqrt(magnsq(Gmetric, nc));
                rho_max = filter_a * exp(filter_b*log(Gabs));
                if (rho_magn > rho_max) {
                    rhoR[nind-1] *= rho_max/rho_magn;
                    rhoI[nind-1] *= rho_max/rho_magn;
                }
            }
}

#endif