//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include "io.H"
//#include <string.h>
//#include <math.h>
#include <unistd.h>
#include "chgcar.H"
#include "3dfft.H"
#include "matrix.H"
#include "smoothing.H"

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "[-rvt] CHGCAR.bulk CHGCAR.new [rcut]";

const char* ARGEXPL = 
"  CHGCAR.bulk: output from VASP; handles US, PAW, spin- and non-spin-polar\n\
  CHGCAR.new:  header portion of CHGCAR, including grid\n\
  rcut:        cutoff distance for being \"vacuum\"; default = WS radius\n\
\n\
  -k  keep charge density; no \"resetting\" to correct amount\n\
  -r  \"raw\" density -- no smoothing\n\
  -v  VERBOSE\n\
  -t  TESTING\n\
\n\
** WARNING: CHGCAR.bulk and CHGCAR.new should _not_ both be stdin";

const char* FILEEXPL =
"Both CHGCAR.bulk and CHGCAR.new are in VASP format; CHGCAR.new only requires\n\
the atom positions and the new grid dimensions, while CHGCAR.bulk must be a\n\
complete CHGCAR file, output from VASP to be meaningful.";

int VERBOSE = 0;	// The infamous verbose flag.
int TESTING = 0;	// Extreme verbosity (testing purposes)
int ERROR = 0;		// Analysis: Error flag (for analysis purposes)

int main(int argc, char **argv) {

    // ************************** INITIALIZATION ***********************
    char* progname = basename(argv[0]);
    int RESCALE = 1;	// rescale total charge
    int SMOOTH = 1;	// smooth total charge
    int GAUSSIAN = 1; // use Gaussian filter (1) or Power law truncation (0)

    char ch;
    while ((ch = getopt(argc, argv, "krvthp")) != -1) {
        switch (ch) {
            case 'k': RESCALE = 0; break;
            case 'r': SMOOTH = 0; break;
            case 'v': VERBOSE = 1; break;
            case 't': TESTING = 1; VERBOSE = 1; break;
            case 'h':
            case '?':
            case 'p': GAUSSIAN = 0; break;
            default: ERROR = 1;
        }
    }
    argc -= optind; argv += optind;
    if (argc < NUMARGS) ERROR = 2;
    // All hell broken loose yet?
    if (ERROR != 0) {
        fprintf(stderr, "%s %s\n", progname, ARGLIST);
        fprintf(stderr, "%s\n", ARGEXPL);
        fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
        exit(ERROR);
    }

    // ****************************** INPUT ****************************
    char dump[512];
    FILE* infile_bulk, *infile;
    char* chgcar_bulk_name = argv[0];
    char *chgcar_name = argv[1];
    double rcut = 0;
    if (argc > NUMARGS) rcut = strtod(argv[2], (char**)NULL);
    if (rcut < 0) {
        fprintf(stderr, "rcut=%g < 0?\n", rcut);
        exit(ERROR_BADFLAG);
    }

    infile_bulk = myopenr(chgcar_bulk_name);
    if (infile_bulk == NULL) {
        fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_bulk_name);
        exit(ERROR_NOFILE);
    }

    infile = myopenr(chgcar_name);
    if (infile == NULL) {
        fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_name);
        exit(ERROR_NOFILE);
    }

    if (infile_bulk == infile) {
        fprintf(stderr, "Cannot have both CHGCAR.bulk and CHGCAR.new be stdin.\n");
        exit(ERROR_NOFILE);
    }

    // ***************************** ANALYSIS **************************
    double a0;
    double latt_bulk[9], latt[9];
    int atomlist[100];
    int Natom_bulk=0, Natom=0;
    double **uatom_bulk = NULL, **uatom = NULL;
    double rsmall;

    // read header
    read_CHGCAR_header(infile_bulk, dump, a0, latt_bulk, atomlist, uatom_bulk);
    for (int i=0; atomlist[i]>=0; ++i) Natom_bulk += atomlist[i];
    for (int d=0; d<9; ++d) latt_bulk[d] *= a0;

    read_CHGCAR_header(infile, dump, a0, latt, atomlist, uatom);
    write_CHGCAR_header(stdout, dump, a0, latt, atomlist, uatom);
    for (int i=0; atomlist[i]>=0; ++i) Natom += atomlist[i];
    for (int d=0; d<9; ++d) latt[d] *= a0;

    // cutoff:
    if (rcut == 0) 
        // Wigner-Seitz radius:
        rcut = cube_root( 3./(4.*M_PI) *det(latt_bulk) );
    // minimum distance:
    rsmall = 0.5*min_nn(latt, Natom, uatom);

    // read grid values
    int N[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%d %d %d", N+0, N+1, N+2);
    int Ng=N[0]*N[1]*N[2];
    double* rho = new double[Ng]; // charge density
    read_grid(infile, rho, Ng);
    myclose(infile);

    // non-spin polarized case, hard-coded though...
    for (int ISPIN=0; ISPIN<1; ++ISPIN) {
        //++ grid values
        int Nb[3];
        fgets(dump, sizeof(dump), infile_bulk);
        sscanf(dump, "%d %d %d", Nb+0, Nb+1, Nb+2);

        int Nbg=Nb[0]*Nb[1]*Nb[2];
        double* rhob = new double[Nbg];

        // ++ read in our bulk grid
        read_grid(infile_bulk, rhob, Nbg);

        // Now, all of our analysis.
        // 1. determine our low-pass filter values (from bulk)
        if (TESTING) {
            double cell=0;
            for (int ng=0; ng<Nbg; ++ng) cell += rhob[ng];
            fprintf(stderr, "# sum rhob= %g x %d x %d\n", 
                cell/(Natom_bulk*Nbg), Natom_bulk, Nbg);
            int negcount=0;
            for (int ng=0; ng<Nbg; ++ng) if (rhob[ng]<0) ++negcount;
            fprintf(stderr, "# negcount= %d/%d = %g\n", 
                negcount, Nbg, (double)negcount/(double)Nbg);
        }
        double l0, filter_a, filter_b;
        double cell_charge;
        const double FILTER_MULT = 0.45;  // at some point, should make this a param...
        {
            double* rhobR = new double[Nbg];
            double* rhobI = new double[Nbg];
            forward_3dfft(rhob, Nb[0], Nb[1], Nb[2], rhobR, rhobI);
            cell_charge = rhobR[0]/Nbg; // need this for later
            if (GAUSSIAN)
                extract_filter_gaussian(rhobR, rhobI, Nb, latt_bulk, l0);
            else
                extract_filter(rhobR, rhobI, Nb, latt_bulk, filter_a, filter_b);
            if (VERBOSE) {
                fprintf(stderr, "# rsmall= %g\n", rsmall);
                fprintf(stderr, "# rcut=   %g\n", rcut);
                if (GAUSSIAN) {
                    fprintf(stderr, "# rho l0= %g\n", l0);
                    fprintf(stderr, "# mult_fact= %.8lf\n", FILTER_MULT);
                }
                else
                    fprintf(stderr, "# rho_filter(G)= %.8le |G|^ %.8lf\n", filter_a, filter_b);
                fprintf(stderr, "# scale_fact= %d\n", Natom);
            }
            // garbage collection
            delete[] rhob; delete[] rhobR; delete[] rhobI;
        }

        // 2. low-pass truncation filter
        double* rhoR = new double[Ng];
        double* rhoI = new double[Ng];
        double basescale = (double)Ng * (double)Natom; // needed due to size
        if (TESTING) {
            double cell=0;
            for (int ng=0; ng<Ng; ++ng) cell += rho[ng];
            fprintf(stderr, "# sum rho=  %g x %d x %d\n", 
                    cell/basescale, Natom, Ng);
            int negcount=0;
            for (int ng=0; ng<Ng; ++ng) if (rho[ng]<0) ++negcount;
            fprintf(stderr, "# negcount= %d/%d = %g\n", 
                    negcount, Ng, (double)negcount/(double)Ng);      
        }
        forward_3dfft(rho, N[0], N[1], N[2], rhoR, rhoI);
        
        if (VERBOSE) {
            fprintf(stderr, "# rhoR[G=0]= %g x %d x %d\n", 
                    rhoR[0]/basescale, Natom, Ng);
            fprintf(stderr, "# should be= %g\n", cell_charge);
            fprintf(stderr, "# RESCALE= %d\n", RESCALE);
            fprintf(stderr, "# SMOOTH= %d\n", SMOOTH);
        }
        // scale values to correct charge / spin:
        if (RESCALE) {
            double scale = cell_charge * basescale/ rhoR[0];
            for (int ng=0; ng<Ng; ++ng) {
                rhoR[ng] *= scale; rhoI[ng] *= scale;
            }
        }

        // smooth it out
        if (SMOOTH) {
            if (GAUSSIAN)
                apply_filter_gaussian(rhoR, rhoI, N, latt, l0/FILTER_MULT);
            else
                apply_filter(rhoR, rhoI, N, latt, filter_a, filter_b);
        }
        inverse_3dfft(rhoR, rhoI, N[0], N[1], N[2], rho);

        // garbage collection
        delete[] rhoR; delete[] rhoI;

        // quick sanity check against "negative" charge:
        if (ISPIN == 0) {
            for (int ng=0; ng<Ng; ++ng)
                if (rho[ng]<0) rho[ng]=0;
            if (RESCALE) {
                double cell=0;
                for (int ng=0; ng<Ng; ++ng) cell += rho[ng];
                if (TESTING) {
                    fprintf(stderr, "# sum rho=  %g x %d x %d\n", 
                    cell/basescale, Natom, Ng);
                }
                double scale = cell_charge * basescale/ cell;
                for (int ng=0; ng<Ng; ++ng) rho[ng] *= scale;
            }
        }

        // *** OUTPUT ***
        printf(" %4d %4d %4d\n", N[0], N[1], N[2]);
        write_grid(stdout, rho, Ng);
    
        if (feof(infile_bulk)) {
            myclose(infile_bulk);
            return ERROR;
        }

        //++ (augmentation charges?)
        fgets(dump, sizeof(dump), infile_bulk);
        if (dump[0] == 'a') {
            // augmentation charges...
            int aug_grid;
            sscanf(dump, "%*s %*s %*d %d", &aug_grid);
            double* aug = new double[aug_grid];
            read_grid(infile_bulk, aug, aug_grid);
            
            for (int n=0; n<Natom; ++n) {
                printf("augmentation occupancies %3d %3d\n", n+1, aug_grid);
                write_grid(stdout, aug, aug_grid);
            }
            fgets(dump, sizeof(dump), infile_bulk);
        }

        if (feof(infile_bulk)) {
            myclose(infile_bulk);
            return ERROR;
        }

    }
    myclose(infile_bulk);
    return ERROR;
}
