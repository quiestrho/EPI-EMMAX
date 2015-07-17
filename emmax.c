#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#include "kmacros.h"

#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
//#include <cblas.h>
//#include <clapack.h>
#include <zlib.h>
//#include "lapack_wrapper.h"

#include "mkl.h"

// Constants for I/O routines
#define DEFAULT_TFAM_NCOLS 6
#define DEFAULT_TPED_NCOLS 4
#define DEFAULT_NDIGITS 5
#define DEFAULT_DELIMS " \t\r\n"
#define SZ_LONG_BUF 1000000
#define DEFAULT_SIZE_MATRIX 1000000
#define DEFAULT_SIZE_HEADER 100000
#define MAX_NUM_MARKERS 20000000 // 20M max # snps
#define SZBUF 1024
#define DBL_MISSING -1e99 
#define DEFAULT_TPED_NUM_HEADER_COLS 4
#define DEFAULT_TFAM_NUM_HEADER_COLS 6
#define DEFAULT_TPED_SNPID_INDEX 1
#define DEFAULT_PHENO_NUM_HEADER_COLS 2
#define KINSHIP_IBS_MEAN 1
#define KINSHIP_IBS_RAND 2
#define KINSHIP_BALDING_NICHOLS 3

// Constants for numerical routines
#define FPMIN 1e-30
#define MAXIT 1000
#define EIGEN_EPS 1e-10
#define EIGEN_ADD 1.
#define TOL 1e-10
#define EPS 1e-10
#define DEFAULT_NGRIDS 100
#define DEFAULT_LLIM -10
#define DEFAULT_ULIM 10
#define XACCU 1e-4
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

struct HFILE {
  int gzflag;       // 1 if gz if used
  int wflag;        // r(0)/w(1) for plain, rb(0)/wb(1) for gz
  int nheadercols;  // # of header columns (0 if nrows=0)
  int nvaluecols;   // # of value cols (0 if nrows=0)
  int nrows;        // # of rows
  FILE* fp;         // plain file handle
  gzFile gzfp;      // gzip file handle
} g_logh;

static int g_verbose = 0;
static char ct = 'T';
static char cn = 'N';
static char cv = 'V';
static char cl = 'L';
static double onef = 1.0;
static double zerof = 0.0;
static double minusonef = -1.0;
static int onen = 1;
static int dim_int;
static int dim_f;

// Input routines
void close_file (struct HFILE* fhp);
struct HFILE open_file(char* filename, int gzflag, int wflag);
struct HFILE open_file_with_suffix(char* prefix, char* suffix, int gzflag, int wflag);
void read_matrix_with_col_headers( struct HFILE* fhp, int nheadercols, char* delims, int zero_miss_flag, int symmetric, int* p_nmiss, double** matrix, char*** headers);
double* tokenize_line_with_col_headers( struct HFILE* fhp, int nheadercols, char* delims, int zero_miss_flag, char* lbuf, double* values, char** headers, int* p_nvalues, int* p_nmiss );

// Output routines
void emmax_error( const char* format, ... );
void emmax_log( const char* format, ... );
void print_help(void);

// Kinship routines
double update_kinship_IBS_mean( double* kins, double* snps, int n );
double update_kinship_IBS_rand( double* kins, double* snps, int n );
double update_kinship_Balding_Nichols( double* kins, double* snps, int n );
void symmetrize_and_normalize_kinship( double* kins, double sum, int n );

// REMLE routines
void ensure_symmetry_and_relax_diag( double* mat, int n, double tol, double eps );
int eigen_L_wo_Z(int n, double* kins, double* eLvals, double* eLvecs);
int eigen_R_wo_Z(int n, int q, double* kins, double* xs, double* eRvals, double* eRvecs);
double REMLE_grid_wo_Z(int n, int q, double* ys, double* xs, double* kins, int ngrids, double llim, double ulim, double* evals, double* evecs, double* optdelta, double* optvg, double* optve, double* pREML0);
double centered_trace(double* kin, int n);
int matrix_invert(int n, double *X, double *Y);
int eigen_decomposition(int n, double* X, double *eigvec, double *eigval);


// GLS routines
void fill_XDX_X0 ( double* X0, double* eLvals, double delta, int n, int q0, int, double* XDX);
void fill_XDX_X1 ( double* X0, double* x1, double* eLvals, double delta, int n, int q0, int, double* XDX );
void fill_XDy_X0 ( double* X0, double* y, double* eLvals, double delta, int n, int q0, double* XDy );
void fill_XDy_X1 ( double* x1, double* y, double* eLvals, double delta, int n, int q0, int, double* XDy );
double compute_yDy ( double* y, double *eLvals, double delta, int n );
double mult_vec_mat_vec ( double* vec, double* mat, int n );

// stat
double tcdf(double t, double nu);

static void compute_phenotype_summary_stats(double *rval_mean,
					    double *rval_var,
					    double *phenotypes,
					    int n) {
  int i;
  double m, msq;

  m = 0.;
  msq = 0.;
  for(i=0; i < n; i++) {
    m += phenotypes[i];
    msq += phenotypes[i]*phenotypes[i];
  }
  m /= n;
  
  *rval_mean = m;
  *rval_var = msq/n - m*m;
}

static void standardize_phenotype(double *phenotypes, int n, double m, double s) {
  int i;
  
  for(i=0; i < n; i++) 
    phenotypes[i] = (phenotypes[i] - m)/s;
}

static double compute_regressed_phenotype_variance(double *phenotypes,
						   int n,
						   double *betas,
						   double *covariates,
						   int n_covariates) {
  int i, j;
  double m, msq;

  m = 0.; msq = 0.;
  for(i=0; i < n; i++) {
    double p = phenotypes[i];

    for(j=0; j < n_covariates; j++)
      /* covariates (X0 below) is column major format in this version, but I 
         see commented out code below for row major format */
      p -= betas[j]*covariates[i + j*n];
    
    m += p;
    msq += p*p;
  }

  m /= n;
  return msq/n - m*m;
}

enum { GT_AA = 0, GT_AB, GT_BB, GT_MISSING, N_GT };
static double genotype_values[N_GT] = { -1.0, 0.0, 1.0, DBL_MISSING };

static int encode_genotype(char *q1, char *q2) {
  if (q1[0] == '0' || q2[0] == '0') {
    return GT_MISSING;
  } else if ((q1[0] == '1' && q2[0] == '2') ||
		 (q1[0] == '2' && q2[0] == '1')) {
    return GT_AB;
  } else if (q1[0] == '1' && q2[0] == '1') {
    return GT_AA;
  } else if (q1[0] == '2' && q2[0] == '2') {
    return GT_BB;
  } else {
    fprintf(stderr,"Covariate genotype %s/%s is not understood\n",
	    q1, q2);
    return GT_MISSING;
  }
}

static double compute_allele_frequency(int *gt, int n) {
  int i, n_obs;
  double b_freq;

  b_freq = 0.;
  n_obs = 0;
  for(i=0; i < n; i++) {
    if (gt[i] != GT_MISSING) {
      n_obs++;
    
      if (gt[i] == GT_AB) {
	b_freq += 0.5;
      } else if (gt[i] == GT_BB) {
	b_freq += 1.0;
      }
    }
  }

  b_freq /= n_obs;

  return b_freq;
}

static void set_genotype_values(double *gv, int *g, int n) {
  int i, nm;
  double m;

  nm = 0;
  m = 0.;
  for(i=0; i < n; i++) {
    if (g[i] != GT_MISSING) {
      gv[i] = genotype_values[ g[i] ];
      m += gv[i];
      nm++;
    }
  }

  m /= nm;
  for(i=0; i < n; i++)
    if (g[i] == GT_MISSING) gv[i] = m;
}

static void set_interaction_values(double *gi, int *g1, int *g2, int n) {
  int i;

  for(i=0; i < n; i++) gi[i] = 0;
  
  for(i=0; i < n; i++) {
    if (g1[i] != GT_MISSING || g2[i] != GT_MISSING) {
      if (g1[i] == GT_AA && g2[i] == GT_AA) {
	gi[i] = -1.;
      } else if (g1[i] == GT_AA && g2[i] == GT_BB) {
	gi[i] = 1.;
      } else if (g1[i] == GT_BB && g2[i] == GT_AA) {
	gi[i] = 1.;
      } else if (g1[i] == GT_BB && g2[i] == GT_BB) {
	gi[i] = -1.;
      }
    }
  }
}

static int extract_covariate_snpids(char ***r_ids, char *snpid_list) {
  char *p, *q, **snp_ids;
  int n_snps;
  
  if (snpid_list == NULL) {
    *r_ids = NULL;
    return 0;
  }

  snp_ids = NULL;
  n_snps = 0;
  p = snpid_list;
  while((q = strsep(&p, ",")) != NULL) {
    if (n_snps % 8 == 0)
      snp_ids = realloc(snp_ids, sizeof(char *)*(n_snps + 8));

    snp_ids[n_snps] = strdup(q);
    n_snps++;
  }

  *r_ids = snp_ids;
  return n_snps;
}


/* KONI - 2015-07-10
   10 megabytes max input line. Properly, should be dynamically set based on
   actual input size.

   TODO: Implement using gzip'd tped file as the rest of the program does */
#define INPUT_LINEMAX (10*1048576)
static int *load_covariate_marker(int n_indv, char *tped_basename,
				  char *covariate_snpid) {
  char *tped_fname;
  char *inputbuf, *p, *q1, *q2, *snp_id;
  int *covariate_genotypes;
  FILE *f;
  int i, n;
  
  tped_fname = (char *) malloc(sizeof(char)*(strlen(tped_basename) + 10));
  sprintf(tped_fname,"%s.tped", tped_basename);
  
  /* Add code to read gzip'd file */
  f = fopen(tped_fname, "r");
  if (f == NULL) {
    fprintf(stderr,"Can't open input file %s (%s)\n", tped_fname, strerror(errno));
    exit(-1);
  }

  inputbuf = (char *) malloc(sizeof(char)*INPUT_LINEMAX);
  covariate_genotypes = (int *) malloc(sizeof(int)*n_indv);
  n = 0;
  while(fgets(inputbuf, INPUT_LINEMAX, f) != NULL) {
    CHOMP(inputbuf);
    if (inputbuf[0] == 0) continue;

    p = inputbuf;
    strsep(&p, " \t");
    snp_id = strsep(&p, " \t");
    if (strcmp(snp_id, covariate_snpid) != 0) continue;

    strsep(&p, " \t");
    strsep(&p, " \t");
    while((q1 = strsep(&p, " \t")) != NULL) {
      q2 = strsep(&p, " \t");
      if (n >= n_indv) {
	fprintf(stderr,"Error: expecting %d individuals for covariate SNP %s but found more\n",
		n_indv, covariate_snpid);
      }
      covariate_genotypes[n] = encode_genotype(q1, q2);
      n++;
    }
    
    if (n < n_indv) {
      fprintf(stderr,"Error: expecting %d individuals for covariate SNP %s but found %d\n",
	      n_indv, covariate_snpid, n);
      exit(-1);
    }
    break;
  }
  fclose(f);

  if (n == 0) {
    fprintf(stderr,"Error: covariate SNP %s was not found in %s\n", covariate_snpid, tped_fname);
    exit(-1);
  }
  
  return covariate_genotypes;
}

/* KONI - 2015-07-11 - added this routine for aiding debugging */
static void print_matrix(FILE *f, char *label, double *mat, int n, int m) {
  int j,k;

  fprintf(f,"%s\n",label);
  for(j=0; j < n; ++j) {
    int k;
    for(k=0; k < m; k++) 
      fprintf(f," %1.2e\t",mat[j + k*n]);
    fprintf(f,"\n");
  }
}

static double calculate_rsq(int *a, int *b, int nf) {
  double n, p1, p2, x00, x01, x10, x11, D;
  int i;
  
  n = 0.;
  p1 = 0.;
  p2 = 0.;
  x00 = x01 = x10 = x11 = 0.;
  for(i=0; i < nf; i++) {
    if (a[i] == GT_MISSING || b[i] == GT_MISSING) continue;
    if (a[i] == GT_AB && b[i] == GT_AB) continue;
    if (a[i] == GT_AA && b[i] == GT_AA) {
      x00+=2;
    } else if (a[i] == GT_AA && b[i] == GT_AB) {
      x00++;
      x01++;
    } else if (a[i] == GT_AA && b[i] == GT_BB) {
      x01+=2;
    } else if (a[i] == GT_AB && b[i] == GT_AA) {
      x00++;
      x10++;
    } else if (a[i] == GT_AB && b[i] == GT_AB) {
      /* Need EM to resolve this */
    } else if (a[i] == GT_AB && b[i] == GT_BB) {
      x01++;
      x11++;
    } else if (a[i] == GT_BB && b[i] == GT_AA) {
      x10+=2;
    } else if (a[i] == GT_BB && b[i] == GT_AB) {
      x10++;
      x11++;
    } else if (a[i] == GT_BB && b[i] == GT_BB) {
      x11+=2;
    }
    n+=2;
  }
  p1 = (x00 + x01)/n;
  p2 = (x00 + x10)/n;
  x00 /= n;
  x01 /= n;
  x10 /= n;
  x11 /= n;
  D = x00 - p1*p2;
  return (D*D/(p1*(1.-p1)*p2*(1.-p2)));
}

static void print_output_header(FILE *f, char **covariate_snpids, int n_covariate_snps,
				int do_genetic_interaction_tests,
				int n_file_covariates) {
  int k;
  
  /* The simplest part, the tested marker in the scan */
  fprintf(f,"Test marker\tchm\tpos\tPVE(model)\tPVE(genetic)\teffect(scan)\tp-value\tallele_freq\tPVE");

  for(k=0; k < n_covariate_snps; k++) {
    fprintf(f,"\t%s\tp-value\trsq\tPVE", covariate_snpids[k]);
  }
  if (do_genetic_interaction_tests) {
    for(k=0; k < n_covariate_snps; k++) {
      fprintf(f,"\tI(%s)\tp-value", covariate_snpids[k]);
    }
  }

  /* Fixed covariates provided as seperate input file. Would normally contain
     the intercept term as the first column */
  fprintf(f,"\tintercept\tp-value");
  for(k=1; k < n_file_covariates; k++)
    fprintf(f,"\t(cov %d)\tp-value", k);

  fprintf(f,"\n");
}

int main(int argc, char** argv) {
  int i, j, k, l, n, nf, q0, q, ngrids, ndigits, istart, iend, nelems, nmiss, *wids, c;
  char *kinf, *phenof, *tpedf, *covf, *outf, *inf, *delims, *lbuf;
  int mphenoflag, write_eig_flag, gz_flag, tped_nheadercols, tfam_nheadercols, zero_miss_flag, isnpid, gen_kin_method, gls_flag;
  double *phenos, *covs, *kins, *snps;
  double sum, llim, ulim, stat, p;
  clock_t cstart, cend, sum0, sum1, sum2, clapstart, clapend;
  struct HFILE phenosh, covsh, kinsh, tpedh, tfamh, eLvalsh, eLvecsh, remlh, outh;
  char **tped_headers, **tfam_headers, **phenos_indids, **covs_indids;
  double maf_threshold = 0.;
  char *covariate_snplist = NULL;
  char **covariate_snpids;
  int n_covariate_snps;
  int **covariate_genotypes;
  double **covariate_values;
  double *covariate_afreq;
  double *rsq;
  int do_genetic_interaction_tests = 0;
  
  cstart = clock();

  sum0 = sum1 = sum2 = 0;
  llim = DEFAULT_LLIM;
  ulim = DEFAULT_ULIM;
  ngrids = DEFAULT_NGRIDS;
  mphenoflag = write_eig_flag = gen_kin_method = 0;
  delims = DEFAULT_DELIMS;
  tped_nheadercols = DEFAULT_TPED_NUM_HEADER_COLS;
  tfam_nheadercols = DEFAULT_TFAM_NUM_HEADER_COLS;
  zero_miss_flag = 1;
  isnpid = DEFAULT_TPED_SNPID_INDEX;
  //dosage_flag = 0;
  gls_flag = 1;

  // Read optional parameters
  q0 = i = 0;
  kinf = covf = tpedf = inf = outf = phenof = NULL;
  phenos = covs = kins = NULL;
  tped_headers = tfam_headers = phenos_indids = covs_indids = NULL;
  ndigits = DEFAULT_NDIGITS;
  istart = 0;
  iend = MAX_NUM_MARKERS;
  while ((c = getopt(argc, argv, "c:d:k:K:S:E:vi:o:p:t:wzD:P:F:ZONm:e:G")) != -1 ) {
    switch(c) {
    case 'c': // covariates input file
      covf = optarg;
      break;
    case 'd': // precision of digits
      ndigits = atoi(optarg);
      break;
    case 'k': //kinship file
      kinf = optarg;
      break;
    case 'K': // create kinship matrix with various methods
              // 1 : IBS matrix with avg fill-in
              // 2 : IBS matrix with rand fill-in
              // 3 : Balding-Nichols matrix
      gen_kin_method = atoi(optarg);
      break;
    case 'S': // start SNP index
      istart = atoi(optarg);
      break;
    case 'E': // end SNP index
      iend = atoi(optarg);
      break;
    case 'v': // turn on verbose mode
      g_verbose = 1;
      break;
    case 'i': // input file prefix
      inf = optarg;
      break;
    case 'o': // output file prefix
      outf = optarg;
      break;
    case 'p': // phenotype file
      phenof = optarg;
      break;
    case 't': // prefix for tped/tfam file
      tpedf = optarg;
      break;
    case 'w': // flag for writing eigenvalues/eigenvectors - ignored when -i is used
      write_eig_flag = 1; // write eigen files
      break;
    case 'z': // turn on gzip input mode
      gz_flag = 1;
      break;
    case 'D': // change delimiters : default " \t\n\r"
      delims = optarg;
      break;
    case 'P' : // change number of header columns in TPED file 
      tped_nheadercols = atoi(optarg);
      break;
    case 'F' : // change the number of header columns in TFAM file
      tfam_nheadercols = atoi(optarg);
      break;
    case 'Z' : // dosage mode with one item per SNP
      zero_miss_flag = 0;
      break;
      //case 'O' :
      //dosage_flag = 1;
      //break;
    case 'N' :
      gls_flag = 0;
      break;
      /* KONI - 2015-07-09 - minor allele frequency/count threshold. If < 1.0, treated as a
	 frequency. If > 1.0, treated as a minimum minor allele count */
    case 'm' : 
      maf_threshold = strtod(optarg, NULL);
      break;
      /* KONI - 2015-07-09 - option for specifying a specific marker to use as a fixed
         effect term in combination with every other marker that is tested. The marker
	 specified itself will not be tested at all */
    case 'e':
      covariate_snplist = optarg;
      break;
      /* KONI - 2015-07-10 - option for performing GxG interaction tests when a 
         covariate SNP is specified. */
    case 'G':
      do_genetic_interaction_tests = 1;
      break;
    default:
      fprintf(stderr,"Error : Unknown option character %c",c);
      print_help();
      abort();
    }
  }

  // Sanity check for the number of required parameters
  if ( argc > optind ) {
    print_help();
    exit(0);//abort();
  }

  n_covariate_snps = extract_covariate_snpids(&covariate_snpids, covariate_snplist);
  
  if ( phenof == NULL ) {
    print_help();
    emmax_error("Phenotype file must be specified\n");
  }
  if ( outf == NULL ) {
    print_help();
    emmax_error("Output prefix must be specified\n");
  }
  if ( tpedf == NULL ) {
    print_help();
    emmax_error("TPED file is not specified\n");
  }

  g_logh = open_file_with_suffix(outf,"log",0,1);

  emmax_log("\nReading TFAM file %s.tfam ....\n",tpedf);
  
  tfamh = open_file_with_suffix(tpedf, "tfam", 0, 0);

  read_matrix_with_col_headers( &tfamh, tfam_nheadercols, delims, 0, 0, &nmiss, NULL, &tfam_headers);
  n = tfamh.nrows;

  snps = (double*)malloc(sizeof(double)*n);  // snp matrix
  tped_headers = (char**)malloc(sizeof(char*)*n);
  lbuf = (char*) malloc(sizeof(char*) * SZ_LONG_BUF);

  if ( gen_kin_method > 0 ) {
    if ( kinf != NULL ) {
      print_help();
      emmax_error("Kinship file cannot be specified with gen_kin_flag. [tped_prefix].kinf will be generated automatically");
    }
    if ( tpedf == NULL ) {
      print_help();
      emmax_error("TPED file must be specified with gen_kin_flag");
    }

    emmax_log( "\nReading TPED file %s.tped to generate a kinship matrix...\n",tpedf);
    tpedh = open_file_with_suffix(tpedf, "tped", gz_flag, 0);

    if ( tpedh.gzflag == 0 ) {
      for(i=0; i < istart; ++i) {
	if ( fgets(lbuf, SZ_LONG_BUF, tpedh.fp ) == NULL ) {
	  emmax_error("The input file %s ended reading before %d lines\n",tpedf,istart);
	}
      }
    }
    else {
      for(i=0; i < istart; ++i) {
	if ( gzgets( tpedh.gzfp, lbuf, SZ_LONG_BUF ) == NULL ) {
	  emmax_error("The input file %s ended reading before %d lines\n",tpedf,istart);
	}
      }
    }

    sum = 0;
    kins = (double*) calloc(n*n, sizeof(double));
    tpedh.nheadercols = tped_nheadercols; 
    nmiss = 0;
    while ( tokenize_line_with_col_headers( &tpedh, tped_nheadercols, delims, zero_miss_flag, lbuf, snps, tped_headers, &nelems, &nmiss) != NULL ) {
      if ( tpedh.nrows % 1000 == 0 ) {
	emmax_log("Reading %d-th line of tped file for creating kinship matrix\n",tpedh.nrows);
      }
      if ( ( atoi(tped_headers[0]) == 0 ) || ( atoi(tped_headers[0]) > 22 ) ) {
	nmiss = 0;
	continue;
      }

      switch(gen_kin_method) {
      case KINSHIP_IBS_MEAN:
	sum += update_kinship_IBS_mean( kins, snps, n );
	break;
      case KINSHIP_IBS_RAND:
	sum += update_kinship_IBS_rand( kins, snps, n );
	break;
      case KINSHIP_BALDING_NICHOLS:
	sum += update_kinship_Balding_Nichols( kins, snps, n);
	break;
      }
      nmiss = 0;
    }
    symmetrize_and_normalize_kinship( kins, sum, n );
    close_file(&tpedh);

    kinsh = open_file_with_suffix(tpedf, "kinf", 0, 1);
    for(i=0; i < n; ++i) {
      for(j=0; j < n; ++j) {
	if ( j > 0 ) fprintf(kinsh.fp,"\t");
	fprintf(kinsh.fp,"%-.*lg",ndigits,kins[i+j*n]);
      }
      fprintf(kinsh.fp,"\t");
    }
    close_file(&kinsh);
  }
  else {
    // Read the kinship matrix from the input file
    emmax_log( "\nReading kinship file %s...\n",kinf);
    if ( kinf == NULL ) {
      kinsh = open_file_with_suffix(tpedf,"kinf",0,0);
      //print_help();
      //emmax_error("Kinship file must be specified without gen_kin_flag\n");
    }
    else {
      kinsh = open_file(kinf, 0, 0);
    }

    read_matrix_with_col_headers(&kinsh, 0, delims, 0, 1, &nmiss, &kins, NULL );  

    emmax_log("  %d rows and %d columns were observed with %d missing values.\n",kinsh.nrows,kinsh.nvaluecols,nmiss);
    if ( ( kinsh.nrows != n ) || ( kinsh.nvaluecols != n ) ) {
      emmax_error("ERROR : Number of rows %d or columns %d is different from %d\n",kinsh.nrows,kinsh.nvaluecols,n);
    }
  }
  
  ensure_symmetry_and_relax_diag( kins, n, TOL, EIGEN_EPS );

  // Read the phenotype matrix from the input file
  emmax_log("\nReading the phenotype file %s...\n",phenof);

  phenosh = open_file(phenof, 0, 0);
  read_matrix_with_col_headers(&phenosh, DEFAULT_PHENO_NUM_HEADER_COLS, delims, 0, 0, &nmiss, &phenos, &phenos_indids );
  emmax_log("  %d rows and %d columns were observed with %d missing values\n",phenosh.nrows,phenosh.nvaluecols, nmiss);
  if ( nmiss > 0 ) {
    mphenoflag = 1;
  }
  fprintf(stderr,"nmiss = %d , mphenoflag = %d\n",nmiss, mphenoflag);

  if ( phenosh.nvaluecols != 1 ) {
    emmax_error("The number of columns in the phenotype file must be 1\n");
  }

  if ( phenosh.nrows != n ) {
    emmax_error("The number of rows in tfam (%d) does not match with the number of phenotypes (%d)",n,phenosh.nrows,n);
  }

  if ( covf == NULL ) {
    emmax_log( "\nNo covariates were found... using intercept only \n");
    covs = (double*)malloc(sizeof(double)*n);
    for(i=0; i < n; ++i) {
      covs[i] = 1.;
    }
    q0 = 1;
  }
  else {
    emmax_log( "\nReading covariate file %s...\n",covf);
    covsh = open_file(covf, 0, 0);
    read_matrix_with_col_headers(&covsh, DEFAULT_PHENO_NUM_HEADER_COLS, delims, 0, 0, &nmiss, &covs, &covs_indids );
    //covs = read_matrix(covf, delims, &nrows, &ncols, &nmiss, 0);
    emmax_log("  %d rows and %d columns were observed with %d missing values. Make sure that the intercept is included in the covariates\n",covsh.nrows,covsh.nvaluecols,nmiss);
    q0 = covsh.nvaluecols;
    if ( nmiss > 0 ) {
      emmax_error("At this point, we do not allow missng covariates");
    }
    if ( covsh.nrows != n ) {
      emmax_error("Number of rows %d is different from %d\n",covsh.nrows,n);
    }
  }

  // scan for missing phenotype and reorganize the inputs
  wids = (int*)malloc(sizeof(int)*n);

  double *K, *eLvals, *eLvecs;
  double trSKS;
  double optdelta, optvg, optve, REML, REML0, hg;
  double *X0;
  double *y, *yt;
  int n_genetic_effects = 1 + n_covariate_snps + (do_genetic_interaction_tests ? n_covariate_snps : 0);
  int tmp = q0 + n_genetic_effects;
  double *XDX = (double*)malloc(sizeof(double)*tmp*tmp);
  double *iXDX = (double*)malloc(sizeof(double)*tmp*tmp);
  double *XDy = (double*)malloc(sizeof(double)*tmp);
  double *betas = (double*)malloc(sizeof(double)*tmp);
  double yDy;
  
  // memory allocation
  y = (double*)malloc(sizeof(double)*n);
  yt = (double*)malloc(sizeof(double)*n);

  //fprintf(stderr,"foo\n");

  //fprintf(stderr,"foo %d\n",mphenoflag);

  if ( mphenoflag == 0 ) {
    // no missing phenotypes - use the full variables with _b
    memcpy(y,phenos,sizeof(double)*n);
    nf = n;
    for(i=0; i < n; ++i) {
      wids[i] = i;
    }
    X0 = covs;
    K = kins;
  }
  else {
    // when missing phenotype exists
    //  new y, X0, x1s matrix with wids are generated
    //  variables are replaced with _p
    nf = 0;
    
    // set wids, and subselect y
    for(k=0; k < n; ++k) {
      if ( phenos[k] != DBL_MISSING ) {
	wids[nf] = k;
	y[nf] = phenos[wids[nf]];
	++nf;
      }
    }
    
    // subselect X0 from covs
    X0 = (double*) malloc(sizeof(double) * q0 * nf);
    for(k=0; k < nf; ++k) {
      for(l=0; l < q0; ++l) {
	// X0[k*q0+l] = covs[wids[k]*q0+l]; // RowMajor
	X0[k+l*nf] = covs[wids[k]+l*n]; // ColMajor
      }
    }
    
    // subselect K
    K = (double*) malloc(sizeof(double) * nf * nf);
    for(k=0; k < nf; ++k) {
      for(l=0; l < nf; ++l) {
	K[k+l*nf] = kins[wids[k]+wids[l]*n]; // RowMajor & ColMajor
      }
    }
  }

  /* KONI - this code should be moved to load_covariate_marker */
  if (n_covariate_snps > 0) {
    covariate_genotypes = (int **) malloc(sizeof(int *)*n_covariate_snps);
    covariate_values = (double **) malloc(sizeof(double *)*n_covariate_snps);
    rsq = (double *) malloc(sizeof(double)*n_covariate_snps);
    covariate_afreq = (double *) malloc(sizeof(double)*n_covariate_snps);
    for(i=0; i < n_covariate_snps; i++) {
      covariate_genotypes[i] = load_covariate_marker(n, tpedf, covariate_snpids[i]);
      covariate_afreq[i] = compute_allele_frequency(covariate_genotypes[i], n);
      covariate_values[i] = (double *) malloc(sizeof(double)*n);
      set_genotype_values(covariate_values[i], covariate_genotypes[i], n);
    }
  } else {
    covariate_genotypes = NULL;
    covariate_values    = NULL;
  }
  
  /* Koni - 2014-07-02 

     Added this to have the phenotype total variance in order to calculate 
     how much of the variance is explained by a marker during marker by marker
     tests */
  double phenotype_mean, phenotype_var;
  compute_phenotype_summary_stats(&phenotype_mean, &phenotype_var, y, nf);
  fprintf(stderr,"phenotype mean = %1.2f, standard deviation = %1.2f\n",
	  phenotype_mean, sqrt(phenotype_var));
  //standardize_phenotype(y, nf, phenotype_mean, sqrt(phenotype_var));
  
  cend = clock();
  emmax_log("File reading - elapsed CPU time is %.6lf\n",((double)(cend-cstart))/CLOCKS_PER_SEC);
  cstart = cend;

  if ( inf == NULL ) {
    double *eRvals, *eRvecs;

    fprintf(stderr,"Computing eigenvectors and eigenvalues.... \n");
    eLvals = (double*)malloc(sizeof(double)*nf);
    eLvecs = (double*)malloc(sizeof(double)*nf*nf);
    eRvals = (double*)malloc(sizeof(double)*(nf-q0));
    eRvecs = (double*)malloc(sizeof(double)*(nf-q0)*nf);
    trSKS = centered_trace(K,nf);
    //fprintf(stderr,"%d %d\n",nf,q0);
    // compute eigen_R and eigen_L
    eigen_R_wo_Z(nf, q0, K, X0, eRvals, eRvecs);
    emmax_log("eRvals[0] = %lf, eRvals[n-q-1] = %lf, eRvals[n-q] = %lf",eRvals[0],eRvals[nf-q0-1],eRvals[nf-q0]);

    cend = clock();
    emmax_log("eigen_R - elapsed CPU time is %.6lf",((double)(cend-cstart))/CLOCKS_PER_SEC);
    cstart = cend;

    REML = REMLE_grid_wo_Z(nf, q0, y, X0, K, ngrids, llim, ulim, eRvals, eRvecs, &optdelta, &optvg, &optve, &REML0);

    free(eRvecs);
    free(eRvals);

    cend = clock();
    emmax_log("REMLE (delta=%lf, REML=%lf, REML0=%lf)- elapsed CPU time is %.6lf",optdelta, REML, REML0, ((double)(cend-cstart))/CLOCKS_PER_SEC);
    cstart = cend;

    eigen_L_wo_Z(nf, K, eLvals, eLvecs);

    cend = clock();
    emmax_log("eigen_L - elapsed CPU time is %.6lf\n",((double)(cend-cstart))/CLOCKS_PER_SEC);
    cstart = cend;

    hg = trSKS/((nf-1.)*optdelta+trSKS);

    remlh = open_file_with_suffix(outf,"reml",0,1);
    fprintf(remlh.fp,"%-.*lg\n",ndigits,REML);
    fprintf(remlh.fp,"%-.*lg\n",ndigits,REML0);
    fprintf(remlh.fp,"%-.*lg\n",ndigits,optdelta);
    fprintf(remlh.fp,"%-.*lg\n",ndigits,optvg);
    fprintf(remlh.fp,"%-.*lg\n",ndigits,optve);
    fprintf(remlh.fp,"%-.*lg\n",ndigits,hg);
    fclose(remlh.fp);
    
    if ( write_eig_flag ) {
      struct HFILE okinsh = open_file_with_suffix(outf,"kinf",0,1);
      for(i=0; i < nf; ++i) {
	for(j=0; j < nf; ++j) {
	  if ( j > 0 ) fprintf(okinsh.fp,"\t");
	  fprintf(okinsh.fp,"%-.*lg",ndigits,K[i+j*nf]); // ColMajor
	}
	fprintf(okinsh.fp,"\n");
      }
      fclose(okinsh.fp);

      eLvalsh = open_file_with_suffix(outf,"eLvals",0,1);
      for(i=0; i < nf; ++i) {
	fprintf(eLvalsh.fp,"%-.*lg\n",ndigits,eLvals[i]);
      }
      fclose(eLvalsh.fp);
      
      eLvecsh = open_file_with_suffix(outf,"eLvecs",0,1);
      for(i=0; i < nf; ++i) {
	for(j=0; j < nf; ++j) {
	  if ( j > 0 ) fprintf(eLvecsh.fp,"\t");
	  fprintf(eLvecsh.fp,"%-.*lg",ndigits,eLvecs[i+j*nf]); // ColMajor
	}
	fprintf(eLvecsh.fp,"\n");
      }
      fclose(eLvecsh.fp);
    }
  }
  else {
    eLvals = (double*)malloc(sizeof(double)*nf);
    eLvecs = (double*)malloc(sizeof(double)*nf*nf);

    remlh = open_file_with_suffix(inf,"reml",0,0);
    fscanf(remlh.fp,"%lf",&REML);
    fscanf(remlh.fp,"%lf",&REML0);
    fscanf(remlh.fp,"%lf",&optdelta);
    fscanf(remlh.fp,"%lf",&optvg);
    fscanf(remlh.fp,"%lf",&optve);
    fscanf(remlh.fp,"%lf",&hg);
    fclose(remlh.fp);

    eLvalsh = open_file_with_suffix(inf,"eLvals",0,0);
    for(i=0; i < nf; ++i) {
      if ( fscanf(eLvalsh.fp,"%lf",&eLvals[i]) == EOF ) {
	emmax_error("EOF reached while reading the eigenvector file");
      }
    }
    fclose(eLvalsh.fp);
      
    eLvecsh = open_file_with_suffix(inf,"eLvecs",0,0);
    for(i=0; i < nf; ++i) {
      for(j=0; j < nf; ++j) {
	if ( fscanf(eLvecsh.fp,"%lf",&eLvecs[i+j*nf]) == EOF ) { // ColMajor
	  emmax_error("FATAL ERROR : EOF reached while reading the eigenvector file\n");
	}
      }
    }
    fclose(eLvecsh.fp);

    cend = clock();
    emmax_log("Reading SVD & REML input - elapsed CPU time is %.6lf\n",((double)(cend-cstart))/CLOCKS_PER_SEC);
    cstart = cend;
  }

  // check if indids matches between phenos, covs, and tfam
  for(i=0; i < DEFAULT_PHENO_NUM_HEADER_COLS; ++i) {
    for(j=0; j < n; ++j) {
      if ( strcmp( tfam_headers[j+i*n], phenos_indids[j+i*n]) != 0 ) {
	emmax_error("Individual IDs do not match exactly in order between tfam and phenos file, but do not match at %s line %d : %s vs %s\n",(i==0) ? "FAMID" : "INDID", j+1, tfam_headers[j+i*n], phenos_indids[j+i*n]);
      }
      
      if ( covf != NULL ) {
	if ( strcmp( tfam_headers[j+i*n], covs_indids[j+i*n]) != 0 ) {
	  emmax_error("Individual IDs do not match exactly in order between tfam and covs file, but do not match at %s line %d : %s vs %s\n",(i==0) ? "FAMID" : "INDID", j+1, tfam_headers[j+i*n], covs_indids[j+i*n]);
	}
      }
    }
  }

  if ( gls_flag == 1 ) {
    // transform y_p, X0_p, x1s_p to yt_p, X0t_p, x1ts_p
    double *X0t, *x1, *x1t; 

    X0t = (double*) malloc(sizeof(double) * q0 * nf);
    x1 = (double*)malloc(sizeof(double)*nf*n_genetic_effects);  // snp matrix
    x1t = (double*)malloc(sizeof(double)*nf*n_genetic_effects);  // snp matrix

    fprintf(stderr,"Solving t(eLvecs) %%*%% X0...\n");
    // X0t = t(eLvecs) %*% X0
    //cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, nf, q0, nf, 1.0, eLvecs, nf, X0, nf, 0.0, X0t, nf );
    dgemm(&ct,&cn,&nf,&q0,&nf,&onef,eLvecs,&nf,X0,&nf,&zerof,X0t,&nf);

    // yt = t(eLvecs) %*% y
    //cblas_dgemv(CblasColMajor, CblasTrans, nf, nf, 1., eLvecs, nf, y, 1, 0., yt, 1);
    dgemv( &ct,&nf, &nf, &onef, eLvecs, &nf, y, &onen, &zerof, yt, &onen);
    
    outh = open_file_with_suffix(outf,"ps",0,1);
    
    // symmetric - both RowMajor & ColMajor
    fill_XDX_X0 ( X0t, eLvals, optdelta, nf, q0, n_genetic_effects, XDX );
    fill_XDy_X0 ( X0t, yt, eLvals, optdelta, nf, q0, XDy );
    yDy = compute_yDy ( yt, eLvals, optdelta, nf );

    // Read the snp matrix from the input file
    emmax_log( "\nReading TPED file %s...\n",tpedf);
    tpedh = open_file_with_suffix(tpedf, "tped", gz_flag, 0);
    tpedh.nheadercols = tped_nheadercols;

    if ( tpedh.gzflag == 0 ) {
      for(i=0; i < istart; ++i) {
	if ( fgets(lbuf, SZ_LONG_BUF, tpedh.fp ) == NULL ) {
	  emmax_error("The input file %s ended reading before %d lines\n",tpedf,istart);
	}
      }
    }
    else {
      for(i=0; i < istart; ++i) {
	if ( gzgets( tpedh.gzfp, lbuf, SZ_LONG_BUF ) == NULL ) {
	  emmax_error("The input file %s ended reading before %d lines\n",tpedf,istart);
	}
      }
    }

    nmiss = 0;
    clapstart = clock();

    fprintf(stderr,"Starting marker tests...\n");
    print_output_header(outh.fp, covariate_snpids, n_covariate_snps, do_genetic_interaction_tests, q0);
    double *gv = (double *) malloc(sizeof(double)*n);
    int *genotypes = (int *) malloc(sizeof(int)*n);
    double **interactions = (double **) malloc(sizeof(double)*n_covariate_snps);
    for(k=0; k < n_covariate_snps; k++)
      interactions[k] = (double *) malloc(sizeof(double)*n);

    /* get a base variance with no genetic effects. This is crude, we are going to substitue
       random values for the genetic terms of the model, solve it, and recover the variance
       explained. That variance is just that which is explained by the kinship matrix and
       any fixed covariate terms supplied as a file that are not genetic terms as per those
       specified on the command line, a genome scan marker, or any interaction term. */
    srand(time(NULL));
    int trial;
    double null_model_residual_var = 0.;
    for(trial=0; trial < 1000; trial++) {
      for(k=0; k < n_genetic_effects; k++) {
	for(j=0; j < nf; j++) {
	  double r = rand()/(RAND_MAX + 1.0) < 0.5;
	  x1[j + nf*k] = r < 0.5 ? 0. : ( r < 0.25 ? -1. : 1.);
	}
      }
      dgemm(&ct, &cn, &nf, &n_genetic_effects, &nf, &onef, eLvecs, &nf, x1, &nf, &zerof, x1t, &nf);
      fill_XDX_X1 ( X0t, x1t, eLvals, optdelta, nf, q0, n_genetic_effects, XDX );
      fill_XDy_X1 ( x1t, yt, eLvals, optdelta, nf, q0, n_genetic_effects, XDy );
      matrix_invert( q0+n_genetic_effects, XDX, iXDX );
      q = q0+n_genetic_effects;
      //cblas_dgemv(CblasColMajor, CblasNoTrans, q0+1, q0+1, 1., iXDX, q0+1, XDy, 1, 0., betas, 1);
      dgemv(&cn, &q, &q, &onef, iXDX, &q, XDy, &onen, &zerof, betas, &onen);

      null_model_residual_var += (yDy - mult_vec_mat_vec(XDy, iXDX, q0 + n_genetic_effects))/
	(nf - q0 - n_genetic_effects);
    }
    null_model_residual_var /= trial;
    fprintf(stderr,"Residual variance with no fixed genetic effects is %1.3f, total phenotype var %1.3f\n",
	    null_model_residual_var, phenotype_var);
    
    /* genome scan loop */
    for(i=istart; (i < iend) && ( tokenize_line_with_col_headers( &tpedh, tped_nheadercols, delims, zero_miss_flag, lbuf, snps, tped_headers, &nelems, &nmiss) != NULL ); ++i) {

      /* Do not analyze a SNP as a genome scan marker if it is also one of the
	 covariate SNPs */
      for(k=0; k < n_covariate_snps; k++)
	if (strcmp(tped_headers[isnpid], covariate_snpids[k]) == 0) break;
      if (n_covariate_snps > 0 && k < n_covariate_snps) continue;
      
      clapend = clock();
      sum0 += (clapend-clapstart);

      if ( i % 10000 == 0 ) {
	emmax_log("Reading %d-th SNP and testing association....\n", i);
      }

      if ( nelems != n ) {
	emmax_error("The SNP file %s.tped do not have the adequate number of columns at line %d - (%d vs %d)\n", tpedf, i, nelems, n);
      }

      /* This nonsense is trying to stay compatible with the tokenize_line... () function
	 call in the loop gaurd above. It is difficult to modify that function to use
	 the genotype enum codes because it has a dual purpose for reading a matrix of
	 floating point numbers as well. Really, a whole different file reading interface
	 needs to be written to free us from that. What I want though is to have genotypes
	 consistent between the covariate SNP[s] and the scan SNP so I can check for LD,
	 and also so that I can change the decimal encoding of the variable, for whatever
	 reason in ONE PLACE. Fortunately the function as is returns values (but in doubles)
	 as 0, 1, 2, DBL_MISSING as is, so I just convert that to ints and hold my nose */
      for(j=0; j < n; j++) {
	if (snps[j] != DBL_MISSING) {
	  genotypes[j] = (int) snps[j];
	} else {
	  genotypes[j] = GT_MISSING;
	}
      }
      double allele_freq = compute_allele_frequency(genotypes, n);
      set_genotype_values(gv, genotypes, n);
      for(j=0; j < nf; ++j)
      	x1[j] = gv[wids[j]];
      
      /* KONI - 2015-07-09 - Skip if minor allele count is less than an integer value > 1
         specified on command line (-m) */
      if (maf_threshold >= 1.0 && (allele_freq*nf < maf_threshold ||
      				   nf*(1 - allele_freq) < maf_threshold))
	      continue;

      /* KONI - 2015-07-09 - Skip if minor allele frequency is less than a proportion < 1.0 
	 specified on the command line (-m) */
      if (maf_threshold < 1.0 && (allele_freq < maf_threshold || 1.0 - allele_freq < maf_threshold))
      	  continue;

      for(j=0; j < nf; ++j) {
	if ( x1[j] == DBL_MISSING ) {
	  x1[j] = sum/(double)(nf-nmiss);
	  //fprintf(stderr,"%d %.5lf\n",j,x1[j]);
	}
      }

      for(k=0; k < n_covariate_snps; k++) {
	rsq[k] = calculate_rsq(genotypes, covariate_genotypes[k], n);
	for(j=0; j < nf; j++) x1[j+(k+1)*nf] = covariate_values[k][wids[j]];
      }

      if (do_genetic_interaction_tests) {
	for(k=0; k < n_covariate_snps; k++) {
	  set_interaction_values(interactions[k], genotypes, covariate_genotypes[k], n);
	  for(j=0; j < nf; j++)
	    x1[j + (k+n_covariate_snps+1)*nf] = interactions[k][wids[j]];
	}
      }

      clapstart = clock();
      // x1t = t(eLvecs) %*% x1
      //cblas_dgemv(CblasColMajor, CblasTrans, nf, nf, 1., eLvecs, nf, x1, 1, 0., x1t, 1);
      //    Ta   TB   M    N         K    alpha  A            B        Beta    C    
      dgemm(&ct, &cn, &nf, &n_genetic_effects, &nf, &onef, eLvecs, &nf, x1, &nf, &zerof, x1t, &nf);
      fill_XDX_X1 ( X0t, x1t, eLvals, optdelta, nf, q0, n_genetic_effects, XDX );
      fill_XDy_X1 ( x1t, yt, eLvals, optdelta, nf, q0, n_genetic_effects, XDy );

    
      matrix_invert( q0+n_genetic_effects, XDX, iXDX );
      q = q0+n_genetic_effects;
      //cblas_dgemv(CblasColMajor, CblasNoTrans, q0+1, q0+1, 1., iXDX, q0+1, XDy, 1, 0., betas, 1);
      dgemv(&cn, &q, &q, &onef, iXDX, &q, XDy, &onen, &zerof, betas, &onen);

      double residual_var = (yDy - mult_vec_mat_vec(XDy, iXDX, q0 + n_genetic_effects))/
	(nf - q0 - n_genetic_effects);

      stat = betas[q0]/sqrt( iXDX[q0*(q0+n_genetic_effects)+q0] * residual_var );
      clapend = clock();
      sum1 += (clapend-clapstart);

      clapstart = clock();
      p = tcdf(stat,nf-q0-n_genetic_effects);
      clapend = clock();
      sum2 += (clapend-clapstart);
      
      /* Koni - 2014-07-02 

	 Attempt at a cheap way to describe the percent of phenotypic variance
	 explained by this marker. Under the simple additive quantitative
	 genetic model, Vm = 2pq*(a + d*(q - p))^2, where Vm means genetic
	 variance for this marker, p is the '1' allele frequency and q = 1 - p,
	 -a is the deviation from the mean for an '00' genotype and +a is the
	 deviation from the mean for a '11' genotype, d the deviation from the
	 mean for a '01' or '10' genotype (a heterozygote). 

	 EMMAX assumes no dominance deviation (d = 0), and the genotype 
	 encoding I have switched to -1, 0, 1 above. Originally it was 0, 1, 2 */
      double percent_variance_explained = 2.*allele_freq*(1. - allele_freq)*
	(betas[q0]*betas[q0])/phenotype_var*100.;

      /* Output the SNP id, chomosome, and position of the scan marker for 
	 this model */
      fprintf(outh.fp,"%s",tped_headers[isnpid]);
      fprintf(outh.fp,"\t%s",tped_headers[0]);
      fprintf(outh.fp,"\t%s",tped_headers[3]);

      /* Output the percent variance explained of the whole model and PVE of genetic effects */
      fprintf(outh.fp,"\t%1.1f", (1. - residual_var/phenotype_var)*100.);
      fprintf(outh.fp,"\t%1.1f", (1. - residual_var/null_model_residual_var)*100.);
      
      /* Output the beta, p-value, allele frequency, and PVE for the scan marker */
      fprintf(outh.fp,"\t%-.*lg",ndigits,betas[q0]);
      fprintf(outh.fp,"\t%-.*lg",ndigits,p);
      fprintf(outh.fp,"\t%1.3f", allele_freq);
      fprintf(outh.fp,"\t%1.2f", percent_variance_explained);

      /* Output terms for covariate SNPs that were modeled along with the scan marker */
      for(k=0; k < n_covariate_snps; k++) {
	fprintf(outh.fp,"\t%-.*lg",ndigits,betas[q0+k+1]);
	stat = betas[q0+k+1]/
	  sqrt( iXDX[(q0+k+1)*(q0+n_genetic_effects)+q0+k+1]*
		( yDy - mult_vec_mat_vec( XDy, iXDX, q0+n_genetic_effects ) ) /
		( nf - q0 - n_genetic_effects ) );
	p = tcdf(stat, nf-q0-n_genetic_effects);
	
	fprintf(outh.fp,"\t%-.*lg", ndigits, p);
	fprintf(outh.fp,"\t%1.3f", rsq[k]);
	percent_variance_explained = 2*covariate_afreq[k]*(1. - covariate_afreq[k])*
	  (betas[q0+k+1]*betas[q0+k+1])/phenotype_var*100.;
	fprintf(outh.fp,"\t%1.3f", percent_variance_explained);
      }

      /* If genetic interaction terms were modeled, output their betas and p-values */
      if (do_genetic_interaction_tests) {
	for(k=0; k < n_covariate_snps; k++) {
	  int col = q0 + n_covariate_snps +1 + k;
	  int df = nf - q0 - n_genetic_effects;
	    
	  fprintf(outh.fp,"\t%-.*lg",ndigits,betas[col]);
	  stat = betas[col]/
	    sqrt( iXDX[(col)*(q0+n_genetic_effects)+col]*
		  ( yDy - mult_vec_mat_vec( XDy, iXDX, q0+n_genetic_effects ) ) / df );
	  p = tcdf(stat, df);
	  fprintf(outh.fp,"\t%-.*lg", ndigits, p);
	}
      }

      /* Output external covariate terms, including intercept, in columns following the
	 genetic effects above */
      for(k=0; k < q0; k++) {
	fprintf(outh.fp,"\t%-.*lg",ndigits,betas[k]);	
	stat = betas[k]/sqrt( iXDX[k*(q0+n_genetic_effects) + k]*( yDy - mult_vec_mat_vec( XDy, iXDX, q0+n_genetic_effects ) ) / (nf - q0 - n_genetic_effects));
	p = tcdf(stat, nf-q0-n_genetic_effects);
	fprintf(outh.fp,"\t%-.*lg", ndigits, p);
      }
      fprintf(outh.fp,"\n");
      
      //memset(snps, 0, sizeof(double)*n);
      nmiss = 0;
      clapstart = clock();
    }
    close_file(&tpedh);
    close_file(&outh);

    free(X0t);
    free(x1t);
    free(x1);
  }

  cend = clock();
  emmax_log("GLS association - elapsed CPU time is %.6lf\n",((double)(cend-cstart))/CLOCKS_PER_SEC);
  emmax_log("Breakdown for input file parsing : %.6lf\n",((double)sum0)/CLOCKS_PER_SEC);
  emmax_log("              matrix computation : %.6lf\n",((double)sum1)/CLOCKS_PER_SEC);
  emmax_log("             p-value computation : %.6lf\n",((double)sum2)/CLOCKS_PER_SEC);
  cstart = cend;
    
  if ( mphenoflag == 1 ) {
    free(covs);
    free(K);
    free(X0);
  }
  if ( tped_headers != NULL ) free(tped_headers);
  if ( covf != NULL ) free(covs_indids);
  if ( snps != NULL ) free(snps);
  if ( lbuf != NULL ) free(lbuf);
  if ( tfam_headers != NULL ) free(tfam_headers);
  if ( phenos_indids!= NULL ) free(phenos_indids);
  free(kins);
  free(phenos);
  free(XDX);
  free(iXDX);
  free(XDy);
  free(betas);
  free(wids);
  free(y);
  free(yt);
  free(eLvals);
  free(eLvecs);

  return 0;
}

// ** eig.L = eigen(K)
// return values are stored in vp_eigL (eigenvalues) and mp_eigL (eigenvectors)
int eigen_L_wo_Z(int n, double* kins, double* eLvals, double* eLvecs ) {
  //  return eigen_decomposition(n,kins,eLvecs,eLvals);
  int i;
  double* kins_copy = (double*)malloc(sizeof(double)*n*n);
  memcpy(kins_copy,kins,sizeof(double)*n*n);
  for(i=0; i < n; ++i) {
    kins_copy[i*n+i] += EIGEN_ADD;
  }
  int info = eigen_decomposition(n, kins_copy, eLvecs, eLvals);
  for(i=0; i < n; ++i) {
    eLvals[i] -= EIGEN_ADD;
  }
  free(kins_copy);
  return info;
}

// ** S = I - X(X'X)^{-1}X'
// ** eig = eigen(S(K+I)S)
// ** ids = (eig$values-1) > 0
// ** values = eig$values[ids]
// ** vectors = eig$vectors[,ids]
// return values are stored in vp_eigL (size n-q) and mp_eigL (n * n-q)
int eigen_R_wo_Z(int n, int q, double* kins, double* X, double* eRvals, double* eRvecs)
{
  int i, j, k, l, info;

  double* kins_copy = (double*)malloc(sizeof(double)*n*n);
  memcpy(kins_copy,kins,sizeof(double)*n*n);
  for(i=0; i < n; ++i) {
    kins_copy[i*n+i] += EIGEN_ADD;
  }

  // first compute XtX = t(X) %*% X
  double* XtX = (double*)calloc(q*q, sizeof(double));
  //cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,q,q,n,1.,X,n,X,n,0.,XtX,q);
  dgemm(&ct,&cn,&q,&q,&n,&onef,X,&n,X,&n,&zerof,XtX,&q);

  // Compute iXtX = solve(XtX)
  double* iXtX = (double*)calloc(q*q, sizeof(double));
  matrix_invert(q, XtX, iXtX);

  // Compute XtiXtXX X %*% solve(t(X)%*%X) %*%t(X)

  double* XiXtX = (double*)calloc(n*q, sizeof(double));
  //cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,q,q,1.,X,n,iXtX,q,0.,XiXtX,n);
  dgemm(&cn,&cn,&n,&q,&q,&onef,X,&n,iXtX,&q,&zerof,XiXtX,&n);
  
  // S = I
  double* S = (double*)calloc(n*n,sizeof(double));
  for(i=0; i < n; ++i) {
    S[i+i*n] = 1.;
  }
  // S = -1*(XtiXtX %*%t(X) + 1*(S=I)
  //cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,n,n,q,-1.,XiXtX,n,X,n,1.,S,n);
  dgemm(&cn, &ct, &n, &n, &q, &minusonef, XiXtX, &n, X, &n, &onef, S, &n);
  //fprintf(stderr,"S[0] = %lf, S[1] = %lf, S[n*n-1] = %lf\n",S[0],S[1],S[n*n-1]);

  // SKS = S %*% K %*% S
  double* SK = (double*)calloc(n*n,sizeof(double));
  //cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.,S,n,kins_copy,n,0.,SK,n);
  dgemm(&cn, &cn, &n, &n, &n, &onef, S, &n, kins_copy, &n, &zerof, SK, &n);
  //fprintf(stderr,"kins[0] = %lf, kins[1] = %lf, kins[n*n-1] = %lf\n",kins[0],kins[1],kins[n*n-1]);

  double* SKS = (double*)calloc(n*n,sizeof(double));
  //cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.,SK,n,S,n,0.,SKS,n);
  dgemm(&cn, &cn, &n, &n, &n, &onef, SK, &n, S, &n, &zerof, SKS, &n);
  //ensure_symmetry_and_relax_diag( SKS, n, TOL, TOL );

  double* evals = (double*)malloc(sizeof(double)*n);
  double* evecs = (double*)malloc(sizeof(double)*n*n);
  info = eigen_decomposition(n, SKS, evecs, evals);
  fprintf(stderr,"evals[0] = %lf, evals[1] = %lf, evals[n-1] = %lf\n",evals[0],evals[1],evals[n-1]);

  // Assuming that the eigenvalues are unordered, 
  //   relocate the smallest eigenvalues at the end
  //   relocated the corresponding eigenvectors accordingly
  //   if negative eigenvalues < -(TOL) observed, report fatal errors
  //   time complexity should be q*n
  double* minevals = (double*)malloc(sizeof(double)*q);
  int* minids = (int*)malloc(sizeof(int)*q);

  // initalize the minimum eigenvalues
  for(i=0; i < q; ++i) {
    minevals[i] = DBL_MAX;
    minids[i] = -1;
  }

  // scan from the first element, and compare with the q-th smallest eigenvalues
  for(i=0; i < n; ++i) {
    if ( evals[i] < 0-TOL ) {
      fprintf(stderr,"FATAL ERROR : Eigenvalues of SKS is %lf\n",evals[i]);
    }

    // minevals[0] is the largest one
    if ( evals[i] < minevals[0] ) {
      // find minimum j that evals[i] > minevals[j]
      for(j=0; (evals[i] < minevals[j]) && ( j < q ); ++j);

      // shift all minevals before j
      for(k=1; k < j-1; ++k) {
	minevals[k-1] = minevals[k];
	minids[k-1] = minids[k];
      }
      minevals[j-1] = evals[i];
      minids[j-1] = i;
    }
  }

  // scan if all minevals are between -TOL and TOL
  for(i=0; i < q; ++i) {
    if ( ( minevals[i] > TOL ) || ( minevals[i] < 0-TOL ) ) {
      fprintf(stderr,"FATAL ERROR : Minimum q eigenvalues of SKS is supposed to be close to zero, but actually %lf\n",minevals[i]);
    }
  }

  /*  
  // extract q eigenvalues closest to zero
  double* minabs = (double*)malloc(sizeof(double)*q);
  int* minids = (int*)malloc(sizeof(int)*q);
  for(i=0; i < q; ++i) {
    minabs[i] = DBL_MAX;
    minids[i] = -1;
  }

  for(i=0; i < n; ++i) {
    for(j=0; j < q; ++j) {
      if ( evals[i] < minabs[j] ) {
	for(k=q-1; k > j; --k) {
	  minabs[k] = minabs[k-1];
	  minids[k] = minids[k-1];
	}
	minabs[j] = evals[i];
	minids[j] = i;
	break;
      }
    }
  }*/

  // exclude minidx[0..q] to copy the eigenvalues and eigenvectors
  for(i=0,k=0; i < n; ++i) {
    for(j=0; j < q; ++j) {
      if ( minids[j] == i )
	break;
    }
    if ( j == q ) {
      eRvals[k] = evals[i]-EIGEN_ADD;
      for(l=0; l < n; ++l) {
	// trying to do eRvecs(l,k) = evecs(l,i)
	// make sure that k < n-q
	//	eRvecs[l*(n-q)+k] = evecs[l*n+i]; //RowMajor
	eRvecs[l+n*k] = evecs[l+i*n];
      }
      ++k;
    }
    if ( k > n-q ) {
      fprintf(stderr,"FATAL ERROR : k = %d > n-q = %d\n", k, n-q);
      abort();
    }
    }
  fprintf(stderr,"k = %d,  n-q = %d\n", k, n-q);

  free(XtX);
  free(iXtX);
  free(XiXtX);
  free(S);
  free(SK);
  free(SKS);
  free(evals);
  free(evecs);
  free(minids);
  free(minevals);
  free(kins_copy);

  return info;
}

double LL_REML_wo_Z(int nq, double logdelta, double* lambdas, double* etas) {
  int i;
  double delta = exp(logdelta);
  long double sum1 = 0.;
  long double sum2 = 0.;
  for(i=0; i < nq; ++i) {
    sum1 += etas[i]*etas[i]/(lambdas[i]+delta);
    sum2 += log(lambdas[i]+delta);
  }
  return( 0.5*(nq*(log(nq/(2*M_PI))-1.-log(sum1))-sum2) );
}

double LL_REML_wo_Z_inf_delta(int nq, double* lambdas, double* etas) {
  int i;
  long double sum1 = 0.;
  for(i=0; i < nq; ++i) {
    sum1 += etas[i]*etas[i];
  }
  return( 0.5*(nq*(log(nq/(2*M_PI))-1.-log(sum1))) );
}

double dLL_REML_wo_Z(int nq, double logdelta, double* lambdas, double* etas) {
  int i;
  double delta = exp(logdelta);
  long double sum1 = 0.;
  long double sum2 = 0.;
  long double sum3 = 0.;
  double f;
  for(i=0; i < nq; ++i) {
    f = (etas[i]*etas[i])/(lambdas[i]+delta);
    sum1 += f/(lambdas[i]+delta);
    sum2 += f;
    sum3 += 1./(lambdas[i]+delta);
  }
  return (0.5*(nq*sum1/sum2-sum3));
}

double bis_root_REML_wo_Z(double x1, double x2, double xacc, int nq, double* lambdas, double* etas)
{
  int j;
  double dx, xmid, rtb;
  double f = dLL_REML_wo_Z(nq,x1,lambdas,etas);
  double fmid = dLL_REML_wo_Z(nq,x2,lambdas,etas);
  //fprintf(stderr,"%lg\t%lg\n",f,fmid);
  if ( f*fmid > 0 ) {
    fprintf(stderr,"Error in bis_root_REML_wo_Z : function value has the same sign in both ends - (%lg, %lg)\n",f,fmid);
    abort();
  }
  rtb = f < 0.0 ? (dx = x2-x1,x1) : (dx = x1-x2,x2);
  for(j=0; j < ITMAX; ++j) {
    fmid = dLL_REML_wo_Z(nq,xmid=rtb+(dx *= 0.5),lambdas,etas);
    //fprintf(stderr,"%d\t%lg\t%lg\n",j,xmid,fmid);
    if ( fmid <= 0.0 ) rtb = xmid;
    if ( fabs(dx) < xacc || fmid == 0.0 ) return rtb;
  }
  fprintf(stderr,"Too many bisections in rtbis\n");
  abort();
}
/*
double prop_search_root_REML_wo_Z(double x1, double x2, double xacc, int nq, double* lambdas, double* etas)
{
  int j;
  double dx, xmid, rtb, prop, fpos;

  // compute the function values at the end
  double f = dLL_REML_wo_Z(nq,x1,lambdas,etas);
  double fmid = dLL_REML_wo_Z(nq,x2,lambdas,etas);
  double w = 0.01;
  if ( f*fmid > 0 ) {
    fprintf(stderr,"Error in bis_root_REML_wo_Z : function value has the same sign in both ends - (%lg, %lg)\n",f,fmid);
    abort();
  }

  // rtb is the negative side
  // prop is the (negative/positive)
  rtb = (f < 0.0) ? (dx = x2-x1, fpos = fmid, fmid = f, x1) : (dx = x1-x2, fpos = f, x2);
  
  for(j=0; j < ITMAX; ++j) {
    //    fmid = dLL_REML_wo_Z(nq,xmid=rtb+(dx *= 0.5),lambdas,etas);
    prop = (0-fmid)/(fpos-fmid);
    fmid = dLL_REML_wo_Z(nq,xmid=rtb+(dx *= (w*0.5+(1.-w)*prop)),lambdas,etas);
    fprintf(stderr,"%d\t%lg\t%lg\n",j,xmid,fmid);
    if ( fmid <= 0.0 ) { rtb = xmid; }
    else { fpos = fmid; }
    if ( fabs(dx) < xacc || fmid == 0.0 ) return rtb;
  }
  fprintf(stderr,"Too many bisections in rtbis\n");
  abort();
  }*/

double brent_REML_wo_Z(double ax, double bx, double cx, double tol, double *xmin, int nq, double* lambdas, double* etas)
{
  int iter;
  double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double d=0.0;
  double e=0.0; 
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=LL_REML_wo_Z(nq,x,lambdas,etas);
  for (iter=1;iter<=ITMAX;iter++) { 
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) { 
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=LL_REML_wo_Z(nq,u,lambdas,etas);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
      SHFT(fv,fw,fx,fu)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    } 
  } 
  fprintf(stderr,"Too many iterations in brent\n");
  exit(-1);
}

double REMLE_grid_wo_Z(int n, int q, double* ys, double* xs, double* kins, int ngrids, double llim, double ulim, double* evals, double* evecs, double* optdelta, double* optvg, double* optve, double* pREML0) 
{
  int i;
  double tmp, tmp2;

  // etas = t(eigR) %*% ys
  double* etas = (double*)calloc(n-q,sizeof(double));
  double* grids = (double*)malloc(sizeof(double)*(ngrids+1));
  double* dLLs = (double*)malloc(sizeof(double)*(ngrids+1));

  // etas = t(evecs) %*% ys
  cblas_dgemv(CblasColMajor, CblasTrans, n, n-q, 1., evecs, n, ys, 1, 0., etas, 1);

  for(i=0; i <= ngrids; ++i) {
    grids[i] = llim + (ulim-llim)/ngrids*i;
    dLLs[i] = dLL_REML_wo_Z(n-q, grids[i], evals, etas);
  }

  for(i=0; i < n-q; ++i) {
    emmax_log("etas[%d] = %lf, lambdas[%d] = %lf",i,etas[i],i,evals[i]);
  }

  //fprintf(stderr,"%lg\n",evals[0]);
  //fprintf(stderr,"%lg\n",evecs[0]);

  double maxREML = 0-DBL_MAX;
  double maxREMLlogdelta = 0;

  if ( dLLs[0] < EPS ) {
    tmp = LL_REML_wo_Z(n-q,llim,evals,etas);
    if ( tmp > maxREML ) {
      maxREMLlogdelta = llim;
      maxREML = tmp;
    }
  }
  if ( dLLs[1] > 0-EPS ) {
    tmp = LL_REML_wo_Z(n-q,ulim,evals,etas);
    if ( tmp > maxREML ) {
      maxREMLlogdelta = ulim;
      maxREML = tmp;
    }
  }
  for(i=0; i < ngrids; ++i) {
    if ( ( dLLs[i]*dLLs[i+1] < 0 ) && ( dLLs[i] > 0 ) && ( dLLs[i+1] < 0 ) ) {
      tmp2 = bis_root_REML_wo_Z(grids[i],grids[i+1],XACCU,n-q,evals,etas);
      tmp = LL_REML_wo_Z(n-q,tmp2,evals,etas);
      //fprintf(stderr,"%lg\t%lg\t%lg\n",grids[i],grids[i+1],tmp2);
      if ( tmp > maxREML ) {
	maxREMLlogdelta = tmp2;
	maxREML = tmp;
      }
    }
    emmax_log("%d\t%lf\t%lf\t%lf\n",i,dLLs[i],dLLs[i+1],LL_REML_wo_Z(n-q,grids[i],evals,etas));
  }

  (*optdelta) = exp(maxREMLlogdelta);
  (*optvg) = 0;
  for(i=0; i < n-q; ++i) {
    (*optvg) += etas[i]*etas[i]/(evals[i]+(*optdelta));
  }
  (*optvg) /= (n-q);
  (*optve) = (*optvg) * (*optdelta);
  (*pREML0) = LL_REML_wo_Z_inf_delta(n-q,evals,etas);

  free(etas);
  free(dLLs);
  free(grids);

  return(maxREML);
}

void close_file(struct HFILE* fhp) {
  if ( fhp->gzflag == 1 ) {
    gzclose(fhp->gzfp);
    fhp->gzfp = NULL;
  }
  else {
    fclose(fhp->fp);
    fhp->fp = NULL;
  }
}

// open_file()
// - filename : file name to open
// - gzflag : gzip flag (use gzfp if gzflag=1, otherwise use fp)
// - wflag : write flag (1 if write mode otherwise read mode)
struct HFILE open_file(char* filename, int gzflag, int wflag) {
  struct HFILE fh;
  fh.gzflag = gzflag;
  fh.wflag = wflag;
  fh.nheadercols = 0;
  fh.nvaluecols = 0;
  fh.nrows = 0;
  if ( gzflag == 1 ) {
    char* mode = (wflag == 1) ? "wb" : "rb";
    fh.gzfp = gzopen(filename,mode);
    fh.fp = NULL;

    if ( fh.gzfp == NULL ) {
      emmax_error("Cannot open file %s for %s",filename,(wflag == 1) ? "writing" : "reading");
    }
  }
  else {
    char* mode = (wflag == 1) ? "w" : "r";
    fh.gzfp = (gzFile) NULL;
    fh.fp = fopen(filename,mode);

    if ( fh.fp == NULL ) {
      emmax_error("Cannot open file %s for %s",filename,(wflag == 1) ? "writing" : "reading");
    }
  }
  return fh;
}

// open_file_with_suffix()
// - [prefix].[suffix] : file name to open
// - gzflag : gzip flag (use gzfp if gzflag=1, otherwise use fp)
// - wflag : write flag (1 if write mode otherwise read mode
struct HFILE open_file_with_suffix(char* prefix, char* suffix, int gzflag, int wflag) {
  char filename[SZBUF];
  sprintf(filename,"%s.%s",prefix,suffix);
  return open_file(filename,gzflag,wflag);
}

// fill X0-only part t(X0) %*% D %*% X0 onto XDX matrix
// for computing t(X) %*% D^{-1} %*%X where D is a diagnoal matrix
// XDX : (q0+1)*(q0+1) matrix
// X0 : n * q0 matrix (left q0 cols of X)
// eLvals : size n vector
// delta : double
void fill_XDX_X0 ( double* X0, double* eLvals, double delta, int n, int q0, int n_genetic_effects,
		   double* XDX) {
  int i, j, k;
  double tmp;

  for(i=0; i < q0; ++i) {
    for(j=0; j <= i; ++j) {
      tmp = 0;
      for(k=0; k < n; ++k) {
	tmp += ((X0[k+i*n]*X0[k+j*n])/(eLvals[k]+delta));
      }
      XDX[i+j*(q0+n_genetic_effects)] = tmp;
      if ( j < i )       
	XDX[j+i*(q0+n_genetic_effects)] = tmp;
    }
  }
}

static int double_compare(const void *c, const void *d) {
  const double *a = (const double *) c;
  const double *b = (const double *) d;
  if (a < b) return -1;
  if (a > b) return  1;
  return 0;
}

// fill X1-dependent part for computing 
// XDX = t(X) %*% D %*% X
// where X = [X0 x1], dim(X0) = (n,q), dim(x1) = (n,1)
void fill_XDX_X1 ( double* X0, double* x1, double* eLvals, double delta, int n, int q0, int m,
		   double* XDX ) {
  int i, j, k;
  double tmp[n];

  for(k=0; k < m; k++) {
    for(i=0; i < q0; ++i) {
      for(j=0; j < n; ++j) {
	tmp[j] = ((X0[j+i*n]*x1[j + k*n])/(eLvals[j]+delta));
      }
      qsort(tmp, n, sizeof(double), double_compare);
      for(j=1; j < n; j++) tmp[0] += tmp[j];

      XDX[i+(q0+m)*(q0+k)] = tmp[0]; // fill the last column vector - XDX[i,q0]
      XDX[(q0+k)+(q0+m)*i] = tmp[0]; // fill the last row vector - XDX[q0,i]
    }
  }

  for(i=0; i < m; i++) {
    for(j=0; j <= i; j++) {
      for(k=0; k < n; k++) 
	tmp[k] = x1[k + i*n]*x1[k + j*n]/(eLvals[k]+delta);
      qsort(tmp, n, sizeof(double), double_compare);
      for(k=1; k < n; k++)
	tmp[0] += tmp[k];
      
      XDX[i+q0+(q0+m)*(q0+j)] = tmp[0];
      XDX[j+q0+(q0+m)*(q0+i)] = tmp[0];
    }
  }
  
}

// fill X0-dependent part (first q0 elements) for computing 
// XDy = t(X) %*% D %*% y
// where X = [X0 x1], dim(X0) = (n,q), dim(x1) = (n,1), dim(y) = (n,1)
void fill_XDy_X0 ( double* X0, double* y, double* eLvals, double delta, int n, int q0, double* XDy ) {
  int i, j;
  double tmp[n];

  for(i=0; i < q0; ++i) {
    for(j=0; j < n; ++j) {
      tmp[j]= ((X0[j+n*i]*y[j])/(eLvals[j]+delta));
    }
    qsort(tmp, n, sizeof(double), double_compare);
    for(j=1; j < n; j++) tmp[0] += tmp[j];
    XDy[i] = tmp[0];
  }
}

// fill X1-dependent part (the ladt element) for computing 
// XDX = t(X) %*% D %*% X
// where X = [X0 x1], dim(X0) = (n,q), dim(x1) = (n,1)
void fill_XDy_X1 ( double* x1, double* y, double* eLvals, double delta, int n, int q0, int m,
		   double* XDy ) {
  int i, k;
  double tmp[n];

  for(k=0; k < m; k++) {
    for(i=0; i < n; ++i) {
      tmp[i] = ((x1[i+(n*k)]*y[i])/(eLvals[i]+delta));
    }
    qsort(tmp, n, sizeof(double), double_compare);
    for(i=1; i < n; i++) tmp[0] += tmp[i];
    XDy[q0+k] = tmp[0];
  }
}

// compute yDy = t(y) %*% D %*% y and return it
double compute_yDy ( double* y, double *eLvals, double delta, int n ) {
  int i;
  double tmp[n];
  
  for(i=0; i < n; ++i ) {
    tmp[i] = ((y[i]*y[i])/(eLvals[i]+delta));
  }
  qsort(tmp, n, sizeof(double), double_compare);
  for(i=1; i < n; i++)
    tmp[0] += tmp[i];
  
  return tmp[0];
}

// compute t(vec) %*% mat %*% vec 
// for a general bector and symmetric matrix vec and mat
double mult_vec_mat_vec ( double* vec, double* mat, int n ) {
  int i, j;
  double tmp1[n];
  double tmp2[n];
  
  for(i=0; i < n; ++i) {
    
    for(j=0; j < n; ++j) {
      tmp2[j] = (vec[i]*vec[j]*mat[i+j*n]);
    }
    qsort(tmp2, n, sizeof(double), double_compare);
    for(j=1; j < n; j++)
      tmp2[0] += tmp2[j];
    
    tmp1[i] = tmp2[0];
  }
  
  qsort(tmp1, n, sizeof(double), double_compare);
  for(i=1; i < n; i++)
    tmp1[0] += tmp1[i];
  
  return tmp1[0];
}

// read_matrix_with_col_headers()
// - read [filename] by each line, separate by delim, and store as a row vector and return as one-dimensional array with row count and column count
// - after reading the matrix in a row major order, the matrix is transposed into a column major order and returned
void read_matrix_with_col_headers( struct HFILE* fhp, int nheadercols, char* delims, int zero_miss_flag, int symmetric, int* p_nmiss, double** matrix, char*** headers) {
  int nvalues, i, j, nmiss;
  double *fmat;
  char **fheaders;
  char* lbuf = (char*) malloc(sizeof(char*) * SZ_LONG_BUF);
  int szmat = DEFAULT_SIZE_MATRIX;
  int szheader = DEFAULT_SIZE_HEADER;
  double* cmat = (double*) malloc(sizeof(double) * szmat );
  char** cheaders = (char**) malloc(sizeof(char*) * szheader );

  fhp->nheadercols = nheadercols; 
  nmiss = 0;

  while( tokenize_line_with_col_headers(fhp, nheadercols, delims, zero_miss_flag, lbuf, &cmat[fhp->nrows*fhp->nvaluecols], &cheaders[fhp->nrows*fhp->nheadercols], &nvalues, &nmiss) != NULL ) {
    if ( fhp->nrows == 1 ) {
      fhp->nvaluecols = nvalues;
    }
    else if ( fhp->nvaluecols != nvalues ) {
      emmax_error("The column size %d do not match to %d at line %d\n",nvalues,fhp->nvaluecols,fhp->nrows);
    }

    if ( (fhp->nrows+1)*(fhp->nheadercols) > szheader ) {
      szheader *= 2;
      fprintf(stderr,"Header size is doubled to %d\n",szheader);
      cheaders = (char**) realloc( cheaders, sizeof(char*) * szheader );
    }

    if ( (fhp->nrows+1)*(fhp->nvaluecols) > szmat ) {
      szmat *= 2;
      fprintf(stderr,"Matrix size is doubled to %d\n",szmat);
      cmat = (double*) realloc( cmat, sizeof(double) * szmat );
    }
  }
  free(lbuf);

  *p_nmiss = nmiss;
  
  fmat = (double*) malloc(sizeof(double)*fhp->nrows*fhp->nvaluecols);
  fheaders = (char**) malloc(sizeof(char*)*fhp->nrows*fhp->nheadercols);
  for(i=0; i < fhp->nrows; ++i) {
    for(j=0; j < fhp->nvaluecols; ++j) {
      fmat[i+j*fhp->nrows] = cmat[i*fhp->nvaluecols+j];
    }
    for(j=0; j < fhp->nheadercols; ++j) {
      fheaders[i+j*fhp->nrows] = cheaders[i*fhp->nheadercols+j];
    }
  }
  free(cmat);
  free(cheaders);
  
  if ( matrix != NULL ) {
    if ( *matrix != NULL ) {
      free(*matrix);
    }
    *matrix = fmat;
  }
  
  if ( headers != NULL ) {
    if ( *headers != NULL ) {
      free(*headers);
    }
    *headers = fheaders;
  }
}

// tokenize_line_with_col_headers()
// read one line from the filestream fp, and tokenizes using delim, and store the elements to tgt with the number of elements to p_nelems. NULL is returned if EOF is reached
double* tokenize_line_with_col_headers( struct HFILE* fhp, int nheadercols, char* delims, int zero_miss_flag, char* lbuf, double* values, char** headers, int* p_nvalues, int* p_nmiss ) {

  int j;
  char *token, *next;

  char *ret = (fhp->gzflag == 1) ? gzgets(fhp->gzfp, lbuf, SZ_LONG_BUF) : fgets( lbuf, SZ_LONG_BUF, fhp->fp );
  int nmiss = 0;

  if ( ret == NULL ) {
    return NULL;
  }

  if ( fhp->nheadercols != nheadercols ) {
    emmax_error("# of header columns mismatch (%d vs %d) at line %d",fhp->nheadercols,nheadercols,fhp->nrows);
  }

  //fprintf(stderr,"tokenize-line called %s\n",lbuf);

  token = strtok(lbuf, delims);
  for( j=0; token != NULL; ++j ) {
    if ( j < nheadercols ) {
      headers[j] = strdup(token);
    }
    // if zero_miss_flag is set, assume the genotypes are encoded 0,1,2
    // Additively encodes the two genotypes in the following way
    // when (j-nheadercols) is even, 0->MISSING, add 1->0, 2->1
    // when (j-nheadercols) is odd, check 0-0 consistency, and add 1->0, 2->1
    else if ( zero_miss_flag == 1 ) {
      if ( (j-nheadercols) % 2 == 0 ) {
	if ( strcmp(token,"0") == 0 ) {
	  values[(j-nheadercols)/2] = DBL_MISSING;
	  ++nmiss;
	}
	else if ( strcmp(token,"1") == 0 ) {
	  values[(j-nheadercols)/2] = 0;
	}
	else if ( strcmp(token,"2") == 0 ) {
	  values[(j-nheadercols)/2] = 1;
	}
	else {
	  emmax_error("Only 0,1,2 are expected for genotypes, but %s is observed",token);
	}
      }
      else {
	if ( values[(j-nheadercols)/2] == DBL_MISSING ) {
	  if ( strcmp(token,"0") != 0 ) {
	    emmax_error("0 0 is expected but 0 %s is observed",token); 
	  }
	}
	else if ( strcmp(token,"0") == 0 ) {
	  emmax_error("0 0 is expected but %s 0 is observed",token);
	}
	else if ( strcmp(token,"1") == 0 ) {
	  //values[(j-nheadercols)/2] += 0;
	}
	else if ( strcmp(token,"2") == 0 ) {
	  values[(j-nheadercols)/2] += 1;
	}
      }
    }
    else {
      values[j-nheadercols] = strtod(token, &next);

      if ( token == next ) {
	if ( ( strcmp(token,"NA") == 0 ) || ( strcmp(token,"N") == 0 ) || ( strcmp(token,"?") == 0 ) ) {
	  values[j-nheadercols] = DBL_MISSING;
	  ++nmiss;
	}
	else {
	  emmax_error("Numeric value is expected at line %d and column %d, but observed %s",fhp->nrows,j,token);
	}
      }
    }
    token = strtok(NULL, delims);
  }
  //fprintf(stderr,"tokenize-line ended %d %d\n",j,nheadercols);
  *p_nvalues = (zero_miss_flag == 1 ) ? (j-nheadercols)/2 : (j-nheadercols);
  *p_nmiss = nmiss;
  ++(fhp->nrows);

  if ( j < nheadercols ) {
    emmax_error("Number of header columns are %d, but only %d columns were observed\n", nheadercols, j);
  }

  return values;
}

void ensure_symmetry_and_relax_diag( double* mat, int n, double tol, double eps ) {
  int i,j;
  /*
  fprintf(stderr,"%d\n",n);
  abort();

  for(i=0; i < n; ++i) {
    for(j=0; j < n; ++j) {
      fprintf(stderr," %lf",mat[i*n+j]);
    }
    fprintf(stderr,"\n");
  }
  */
  for(i=0; i < n; ++i) {
    for(j=0; j < i; ++j) {
      if ( fabs(mat[i*n+j]-mat[j*n+i]) > tol ) {
	fprintf(stderr,"FATAL ERROR : Kinship matrix (%d,%d) and (%d,%d) must be the same, but actually different as %lf and %lf\n",i,j,j,i,mat[i*n+j],mat[j*n+i]);
	abort();
      }
      else {
	mat[i*n+j] = (mat[j*n+i] = (mat[i*n+j]+mat[j*n+i])/2.);
      }
    }
    mat[i*n+i] += eps;
  }
}

double betacf(double a, double b, double x) {
  int m,m2; 
  double aa,c,d,del,h,qab,qam,qap;

  qab=a+b; 
  qap=a+1.0; 
  qam=a-1.0; 
  c=1.0; 
  d=1.0-qab*x/qap; if (fabs(d) < FPMIN) d=FPMIN; d=1.0/d;
  h=d; 
  for (m=1;m<=MAXIT;m++) {
    m2=2*m; aa=m*(b-m)*x/((qam+m2)*(a+m2)); 
    d=1.0+aa*d; 
    if (fabs(d) < FPMIN) d=FPMIN; 
    c=1.0+aa/c; 
    if (fabs(c) < FPMIN) c=FPMIN; 
    d=1.0/d; 
    h *= d*c; 
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d; 
    if (fabs(d) < FPMIN) d=FPMIN; 
    c=1.0+aa/c; 
    if (fabs(c) < FPMIN) c=FPMIN; 
    d=1.0/d; 
    del=d*c; 
    h *= del; 
    if (fabs(del-1.0) < EPS) break;
  } 
  if (m > MAXIT) {
    fprintf(stderr,"a or b too big, or MAXIT too small in betacf %lf %lf %lf",a,b,x); 
    abort();
  }
  return h;
}

double gammln(double xx) {
  double x,y,tmp,ser; 
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5}; 
  int j;
  y=x=xx; 
  tmp=x+5.5; 
  tmp -= (x+0.5)*log(tmp); 
  ser=1.000000000190015; 
  for (j=0;j<=5;j++) 
    ser += cof[j]/++y; 
  return -tmp+log(2.5066282746310005*ser/x);
}

double betai(double a, double b, double x) {
  double gammln(double xx); 
  void nrerror(char error_text[]); 
  double bt;
  if (x < 0.0 || x > 1.0) {
    fprintf(stderr,"Bad x in routine betai"); 
    abort();
  }
  if (x == 0.0 || x == 1.0) bt=0.0; 
  else bt=exp((gammln(a+b)-gammln(a)-gammln(b))+(a*log(x))+(b*log(1.0-x)));
  if (x < (a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a; 
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double tcdf(double t, double nu) {
  if ( isnan(t) ) return 1.;
  else return betai(nu/2.,0.5,nu/(nu+t*t));
}

void emmax_error( const char* format, ... ) {
  va_list args;
  fprintf(g_logh.fp, "ERROR: ");
  va_start (args, format);
  vfprintf(g_logh.fp, format, args);
  va_end (args);
  fprintf(g_logh.fp,"\n");

  fprintf(stderr, "ERROR: ");
  va_start (args, format);
  vfprintf(stderr, format, args);
  va_end (args);
  fprintf(stderr,"\n");
  abort();
}

void emmax_log( const char* format, ... ) {
  va_list args;
  va_start (args, format);
  vfprintf(g_logh.fp, format, args);
  va_end (args);
  fprintf(g_logh.fp,"\n");

  if ( g_verbose ) {
    va_start (args, format);
    vfprintf(stderr, format, args);
    va_end (args);
    fprintf(stderr,"\n");
  }
}

// for matrix S = I-11'/n, and K
// compute centered trace tr(SKS)
// this is equivalent to sum(diag(K))-sum(K)/nrow(K)
double centered_trace(double* kin, int n) {
  int i,j;
  double dsum,asum;
  dsum = asum = 0;
  for(i=0; i < n; ++i) {
    for(j=0; j < n; ++j) {
      asum += kin[i+n*j];
      if ( i == j ) {
	dsum += kin[i+n*j];
      }
    }
  }
  return (dsum - asum/(double)n);
}

void print_help(void) {
  fprintf(stderr,"Usage: emmax [options]\n");
  fprintf(stderr,"Required parameters\n");
  fprintf(stderr,"\t-t [tpedf_prefix] : prefix for tped/tfam files\n");
  fprintf(stderr,"\t-o [out_prefix]  : output file name prefix\n");
  fprintf(stderr,"Likely essential parameters\n");
  fprintf(stderr,"\t-p [phenof] : 3-column phenotype file with FAMID, INDID at the first two colmns, in the same order of .tfam file. Not required only with -K option");
  fprintf(stderr,"\t-k [kinf] : n * n matrix containing kinship values in the individual order consistent to [tpedf].tfam file. [tpedf].kinf will be used if not specified\n");
  fprintf(stderr,"\t-c [covf] : multi-column covariate file with FAMID, INDID at the first two colmns, in the same order of .tfam file");
  fprintf(stderr,"Optional parameters\n");
  fprintf(stderr,"Optional parameters\n");
  fprintf(stderr,"\t-i [in_prefix] : input file name prefix including eigenvectors\n");
  fprintf(stderr,"\t-d [# digits]  : precision of the output values (default : 5)\n");
  fprintf(stderr,"\t-s [start index of SNP] : start index of SNP (default : 0)\n");
  fprintf(stderr,"\t-e [end index of SNP] : end index of SNP (default : #snps)\n");
  fprintf(stderr,"\t-w : flag for writing eigenvalues/eigenvector files\n");
  fprintf(stderr,"\t-D [delimiters] : delimter string in quotation marks\n");
  fprintf(stderr,"\t-P [# heaer cols in tped] : # of column headers in tped file\n");
  fprintf(stderr,"\t-F [# heaer cols in tfam] : # of column headers in tfam file\n");
}

void vector_copy(int n, double* X, double* Y) {
  /*
    Copies a vector X of lenght n to vector Y
  */
  int i;
  for (i=0; i<n; i++) Y[i]=X[i];
}

int check_input(double* x, double*y, char* c) {
  int out=0;
  if (x==y) {
    printf("Error in %s: input equals output \n", c);
    out=1;
  }
  return out;
}

// Wrapper for Lapack machine precision function
//static double dlamch (char CMACH)
//{
//  extern  double  dlamch_ (char *CMACHp);
//  return  dlamch_ (&CMACH);
//}

// Wrapper for Lapack eigenvalue function
/*
int dsyevd (char JOBZ, char UPLO, int N,
	    double *A, int LDA, 
	    double *W, 
	    double *WORK, int LWORK, int *IWORK, int LIWORK) 
{
  extern  void  dsyevd_ (char *JOBZp, char *UPLOp, int *Np,
			 double *A, int *LDAp, 
			 double *W, 
			 double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
			 int *INFOp);
  int  INFO;
  dsyevd_ (&JOBZ, &UPLO, &N, A, &LDA, 
	   W, 
           WORK, &LWORK, IWORK, &LIWORK, &INFO);
  
  return  INFO;
  }*/

int eigen_decomposition(int n, double* X, double *eigvec, double *eigval) {
  /*
    This function calculates the eigenvalues and eigenvectors of 
    the n*n symmetric matrix X. 
    The matrices have to be in Fortran vector format.
    The eigenvectors will be put columnwise in the n*n matrix eigvec,
    where the corresponding eigenvalues will be put in the vector 
    eigval (length n of course). Only the lower triangle of the matrix
    X is used. The content of X is not changed.
    
    This function first queries the Lapack routines for optimal workspace 
    sizes. These memoryblocks are then allocated and the decomposition is 
    calculated using the Lapack function "dsyevr". The allocated memory 
    is then freed. 
  */
  
  double *work;
  int *iwork;
  int  info, lwork, liwork;
  
  if (check_input(X, eigvec, "eigen_decomposition")) return 1;
  
  /*  Use a copy of X so we don't need to change its value or use its memoryblock */
  //Xc=malloc(n*n*sizeof(double));
  
  /*  The support of the eigenvectors. We will not use this but the routine needs it  */
  
  /*  Allocate temporarily minimally allowed size for workspace arrays */
  lwork = 2*n*n+6*n+1;
  liwork = 3+5*n;
  iwork = (int*)malloc(sizeof(int)*liwork);
  work = (double*)malloc(sizeof(double)*lwork);
  
  /*  Check for NULL-pointers.  */
  if ((work==NULL)||(iwork==NULL)) {
    printf("malloc failed in eigen_decomposition\n"); 
    return 2;
  }
  
  vector_copy(n*n, X, eigvec);

  dsyevd(&cv,&cl, &n, eigvec, &n, eigval, work, &lwork, iwork, &liwork, &info);
  //sizeWORK = (int)WORK[0]; 
  //sizeIWORK = IWORK[0]; 
  
  /*  Free previous allocation and reallocate preferable workspaces, Check result  */
  //free(WORK);free(IWORK);
  //WORK = malloc (sizeWORK*sizeof(double));
  //IWORK = malloc (sizeIWORK*sizeof(int));

  //if ((WORK==NULL)||(IWORK==NULL)) {
  //  printf("malloc failed in eigen_decomposition\n"); 
  //  return 2;
  //}
  
  /*  Now calculate the eigenvalues and vectors using optimal workspaces  */
  //info = dsyevd('V','L', n, eigvec, n, eigval, WORK, sizeWORK, IWORK, sizeIWORK);
  //eigvec = Xc;

  //info=dsyevr ('V', 'A', 'L', n, Xc, n, 0, 0, 0, 0, dlamch('S'), &numeig, eigval, eigvec, n, ISUPPZ, WORK, sizeWORK, IWORK, sizeIWORK);
  
  /*  Cleanup and exit  */
  //free(WORK); free(IWORK); //free(ISUPPZ); //free(Xc);
  free(work);
  free(iwork);
  return info;
}

int matrix_invert(int n, double *X, double *Y) {
  /*  
      Calculates the inverse of the n*n matrix X: Y = X^-1
      Does not change the value of X, unless Y=X
  */
  
  int info=0;
  
  /*  When X!=Y we want to keep X unchanged. Copy to Y and use this as working variable  */
  if (X!=Y) vector_copy(n*n, X, Y);
  
  /*  We need to store the pivot matrix obtained by the LU factorisation  */
  //int *ipiv;
  int ipiv[n];
  //  ipiv=malloc(n*sizeof(int));
  if (ipiv==NULL) {
    printf("malloc failed in matrix_invert\n"); 
    return 2;
  }
  
  /*  Turn Y into its LU form, store pivot matrix  */
  //info = clapack_dgetrf (CblasColMajor, n, n, Y, n, ipiv);
  dgetrf(&n,&n,Y,&n,ipiv,&info);
  
  /*  Don't bother continuing when illegal argument (info<0) or singularity (info>0) occurs  */
  if (info!=0) return info;
  
  int lwork = n*16;
  double work[lwork];
  //  double* work = malloc(sizeof(double)*lwork);
  /*  Feed this to the lapack inversion routine.  */
  //info = clapack_dgetri (CblasColMajor, n, Y, n, ipiv);
  dgetri(&n,Y,&n,ipiv,work,&lwork,&info);
  
  /*  Cleanup and exit  */
  //  free(ipiv);
  //  free(work);
  return info;
}

// genotype values are encoded from 0 to 2 additively
double update_kinship_IBS_mean( double* kins, double* snps, int n )
{
  int i,j,nv;
  double sum = 0;

  for(i=0, nv=0; i < n; ++i) {
    if ( snps[i] != DBL_MISSING ) {
      sum += snps[i];
      ++nv;
    }
  }

  for(i=0; i < n; ++i) {
    if ( snps[i] == DBL_MISSING ) {
      snps[i] += sum/(double)nv; 
    }
  }

  for(i=0; i < n; ++i) {
    for(j=0; j < i; ++j) {
      kins[i+j*n] += (2.-fabs(snps[i]-snps[j]));
    }
  }
  return 2.;
}

double update_kinship_IBS_rand( double* kins, double* snps, int n ) 
{
  emmax_error("Not implemented yet");
  return 0;
}

double update_kinship_Balding_Nichols( double* kins, double* snps, int n )
{
  emmax_error("Not implemented yet");
  return 0;
}

void symmetrize_and_normalize_kinship( double* kins, double sum, int n )
{
  int i,j;
  for(i=0; i < n; ++i) {
    for(j=0; j < i; ++j) {
      kins[i+j*n] /= sum;
      kins[j+i*n] = kins[i+j*n];
    }
    kins[i+i*n] = 1.;
  }
}
