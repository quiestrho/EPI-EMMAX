#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
//#include <cblas.h>
//#include <clapack.h>
#include <zlib.h>
//#include "lapack_wrapper.h"

#include <set>
#include <string>

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
#define TOL 1e-6
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
// compute C = alpha* [cTransA(A) %*% diag(1/(evals+delta)) %*% cTransB(B)] + beta*C
void compute_XDY(char cTransA, char cTransB, double* A, double* B, int dimAA, int dimAB, int dimBB, double* eVals, double delta, double* C, double alpha, double beta);

double tcdf(double t, double nu);

int main(int argc, char** argv) {
  int i, j, k, l, n, nf, q0, q, ngrids, ndigits, nelems, nmiss, *wids, c, igeno;
  char *kinf, *phenof, *tpedf, *covf, *outf, *snpidf, *holdf, *delims, *lbuf, buf[SZBUF];
  int mphenoflag, indbeta_flag, gz_flag, tped_nheadercols, tfam_nheadercols, zero_miss_flag, isnpid, n_holdout, n_repeat;
  double *phenos, *covs, *kins, *snps, *genos;
  double llim, ulim;
  struct HFILE phenosh, covsh, kinsh, tpedh, tfamh, snpidsh, outh, holdh;
  char **tped_headers, **tfam_headers, **snpid_headers, **hold_headers, **phenos_indids, **covs_indids;
  std::set<std::string> snpidset, holdset;
  struct timeval tv;

  // set default params
  gettimeofday(&tv, NULL);
  srand((unsigned int)tv.tv_usec);

  llim = DEFAULT_LLIM;
  ulim = DEFAULT_ULIM;
  ngrids = DEFAULT_NGRIDS;
  mphenoflag = 0;
  indbeta_flag = 0;
  delims = DEFAULT_DELIMS;
  tped_nheadercols = DEFAULT_TPED_NUM_HEADER_COLS;
  tfam_nheadercols = DEFAULT_TFAM_NUM_HEADER_COLS;
  zero_miss_flag = 1;
  isnpid = DEFAULT_TPED_SNPID_INDEX;

  // Read optional parameters
  q0 = i = 0;
  kinf = covf = tpedf = snpidf = outf = phenof = holdf = NULL;
  phenos = covs = kins = NULL;
  tped_headers = tfam_headers = phenos_indids = covs_indids = snpid_headers = hold_headers = NULL;
  ndigits = DEFAULT_NDIGITS;
  n_holdout = 0;
  n_repeat = 1;
  
  while ((c = getopt(argc, argv, "ih:c:d:k:n:r:vo:p:t:s:zD:P:F:")) != -1 ) {
    switch(c) {
    case 'i':
      indbeta_flag = 1;
      break;
    case 'h':
      holdf = optarg;
      break;
    case 'c': // covariates input file
      covf = optarg;
      break;
    case 'd': // precision of digits
      ndigits = atoi(optarg);
      break;
    case 'k': //kinship file
      kinf = optarg;
      break;
    case 'n':
      n_holdout = atoi(optarg);
      break;
    case 'r':
      n_repeat = atoi(optarg);
      break;
    case 'v': // turn on verbose mode
      g_verbose = 1;
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
    case 's': // prefix for tped/tfam file
      snpidf = optarg;
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
    default:
      fprintf(stderr,"Error : Unknown option character %c",c);
      print_help();
      abort();
    }
  }

  // Sanity check for the number of required parameters
  if ( argc > optind ) {
    print_help();
    abort();
  }

  if ( outf == NULL ) {
    print_help();
    //emmax_error("Output prefix must be specified\n");
    fprintf(stderr,"ERROR: Output prefix must be specified\n");
    abort();
  }

  g_logh = open_file_with_suffix(outf,"log",0,1);

  if ( ( indbeta_flag ) && ( snpidf == NULL ) ) {
    print_help();
    emmax_error("-i option must be set with -s option\n");
  }

  if ( indbeta_flag ) {
    emmax_log("-i option is set");
  }

  if ( phenof == NULL ) {
    print_help();
    emmax_error("Phenotype file must be specified\n");
  }

  if ( tpedf == NULL ) {
    print_help();
    emmax_error("TPED prefix is not specified\n");
  }

  emmax_log("Reading TFAM file %s.tfam ....",tpedf);
  
  tfamh = open_file_with_suffix(tpedf, "tfam", 0, 0);

  read_matrix_with_col_headers( &tfamh, tfam_nheadercols, delims, 0, 0, &nmiss, NULL, &tfam_headers);
  n = tfamh.nrows;

  snps = (double*)malloc(sizeof(double)*n);  // snp matrix
  tped_headers = (char**)malloc(sizeof(char*)*n);
  lbuf = (char*) malloc(sizeof(char*) * SZ_LONG_BUF);

  // Read the kinship matrix from the input file
  emmax_log("Reading kinship file %s...",kinf);
  if ( kinf == NULL ) {
    kinsh = open_file_with_suffix(tpedf,"kinf",0,0);
    //print_help();
    //emmax_error("Kinship file must be specified without gen_kin_flag\n");
  }
  else {
    kinsh = open_file(kinf, 0, 0);
  }
  
  read_matrix_with_col_headers(&kinsh, 0, delims, 0, 1, &nmiss, &kins, NULL );  
  
  emmax_log("\t%d rows and %d columns were observed with %d missing values.",kinsh.nrows,kinsh.nvaluecols,nmiss);
  if ( ( kinsh.nrows != n ) || ( kinsh.nvaluecols != n ) ) {
    emmax_error("ERROR : Number of rows %d or columns %d is different from %d\n",kinsh.nrows,kinsh.nvaluecols,n);
  }
  
  ensure_symmetry_and_relax_diag( kins, n, TOL, EIGEN_EPS );

  // Read the phenotype matrix from the input file
  emmax_log("\nReading the phenotype file %s...",phenof);

  phenosh = open_file(phenof, 0, 0);
  read_matrix_with_col_headers(&phenosh, DEFAULT_PHENO_NUM_HEADER_COLS, delims, 0, 0, &nmiss, &phenos, &phenos_indids );
  emmax_log("  %d rows and %d columns were observed with %d missing values",phenosh.nrows,phenosh.nvaluecols, nmiss);
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
    emmax_log( "No covariates were specified... using intercept only");
    covs = (double*)malloc(sizeof(double)*n);
    for(i=0; i < n; ++i) {
      covs[i] = 1.;
    }
    q0 = 1;
  }
  else {
    emmax_log( "Reading covariate file %s...",covf);
    covsh = open_file(covf, 0, 0);
    read_matrix_with_col_headers(&covsh, DEFAULT_PHENO_NUM_HEADER_COLS, delims, 0, 0, &nmiss, &covs, &covs_indids );
    //covs = read_matrix(covf, delims, &nrows, &ncols, &nmiss, 0);
    emmax_log("  %d rows and %d columns were observed with %d missing values. Make sure that the intercept is included in the covariates",covsh.nrows,covsh.nvaluecols,nmiss);
    q0 = covsh.nvaluecols;
    if ( nmiss > 0 ) {
      emmax_error("At this point, we do not allow missng covariates");
    }
    if ( covsh.nrows != n ) {
      emmax_error("Number of rows %d is different from %d\n",covsh.nrows,n);
    }
  }

  if ( n_holdout == 0 ) {
    // set default # holdout
    n_holdout = n/5;
  }

  if ( holdf != NULL ) {
    emmax_log("bar1");
    holdh = open_file(holdf, 0, 0);
    read_matrix_with_col_headers( &holdh, 1, delims, 0, 0, &nmiss, NULL, &hold_headers);
    emmax_log("holdh.nrows=%d",holdh.nrows);
    for(i=0; i < holdh.nrows; ++i) {
      std::string s(hold_headers[i]);
      holdset.insert(s);
    }
  }

  // read snpid and tped files to add additional fixed effects
  int ngenos = 0;
  if ( snpidf != NULL ) {
    snpidsh = open_file(snpidf, 0, 0);
    read_matrix_with_col_headers( &snpidsh, 1, delims, 0, 0, &nmiss, NULL, &snpid_headers);
    for(i=0; i < snpidsh.nrows; ++i) {
      std::string s(snpid_headers[i]);
      snpidset.insert(s);
    }

    ngenos = snpidset.size();
    emmax_log("ngenos = %d",ngenos);

    genos = (double*)malloc(sizeof(double)*n*ngenos);
    
    tpedh = open_file_with_suffix(tpedf, "tped", gz_flag, 0);
    tpedh.nheadercols = tped_nheadercols;
    
    nmiss = 0;
    for(i=0, igeno=0; ( tokenize_line_with_col_headers( &tpedh, tped_nheadercols, delims, zero_miss_flag, lbuf, snps, tped_headers, &nelems, &nmiss) != NULL ); ++i) {
      if ( i % 10000 == 0 ) {
	emmax_log("Reading %d-th SNP from TPED file....", i);
      }
      
      if ( nelems != n ) {
	emmax_error("The SNP file %s.tped do not have the adequate number of columns at line %d - (%d vs %d)\n", tpedf, i, nelems, n);
      }
      
      // check if the snp id exist in the snpid
      if ( snpidset.find(tped_headers[isnpid]) != snpidset.end() ) {
	double tmpSum = 0, tmpMean;
	int tmpCnt = 0;
	for(j=0; j < n; ++j) {
	  if ( snps[j] != DBL_MISSING ) {
	    tmpSum += snps[j];
	    ++tmpCnt;
	  }
	}
	tmpMean = tmpSum/tmpCnt;
	
	
	for(j=0; j < n; ++j) {
	  if ( snps[j] == DBL_MISSING ) {
	    genos[j+igeno*n] = tmpMean;
	  }
	  else {
	    genos[j+igeno*n] = snps[j];
	  }
	}
	++igeno;
	snpidset.erase(tped_headers[isnpid]);
      }
    }

    if ( igeno != ngenos ) {
      emmax_error("igeno and ngenos are different - check if all snpids are included in the tped file\n");
    }
    
    if ( snpidset.size() > 0 ) {
      emmax_error("snpidset is nonempty - check if all snpids are included in the tped file");
    }
  }
  else {
    emmax_log("WARNING: -s [snpid_file] option has not been specified. Ignoring all fixed SNP effects..");
  }

  double *K;
  double trSKSh;
  double optdelta, optvg, optve, REML, REML0, hg;
  double *X, *Z;
  double *y;
  
  wids = (int*)malloc(sizeof(double)*n);

  if ( mphenoflag == 0 ) {
    // no missing phenotypes - use the full variables with _b
    y = phenos;
    nf = n;
    for(i=0; i < n; ++i) {
      wids[i] = i;
    }
    K = kins;
  }
  else {
    // when missing phenotype exists
    //  new y, X, x1s matrix with wids are generated
    //  variables are replaced with _p
    emmax_log("foo, n=%d, phenos[0]=%lf, phenos[n-1]=%lf\n",n,phenos[0],phenos[n-1]);
    y = (double*)malloc(sizeof(double)*n);
    nf = 0;
    
    // set wids, and subselect y
    for(k=0; k < n; ++k) {
      if ( phenos[k] != DBL_MISSING ) {
	wids[nf] = k;
	y[nf] = phenos[wids[nf]];
	++nf;
      }
    }
    
    emmax_log("nf=%d\n",nf);
    // subselect K
    K = (double*) malloc(sizeof(double) * nf * nf);
    for(k=0; k < nf; ++k) {
      for(l=0; l < nf; ++l) {
	K[k+l*nf] = kins[wids[k]+wids[l]*n]; // RowMajor & ColMajor
      }
    }
  }

  emmax_log("foo1");
  // concatenate covs and genos
  if ( indbeta_flag ) {
    q = q0;
    Z = (double*)malloc(sizeof(double)*nf*ngenos);
  }
  else {
    q = q0 + ngenos;
  }
  X = (double*)malloc(sizeof(double)*nf*q);
  for(k=0; k < nf; ++k) {
    for(l=0; l < q0; ++l) {
      X[k+l*nf] = covs[wids[k]+l*n];
    }
    for(l=0; l < ngenos; ++l) {
      if ( indbeta_flag ) {
	Z[k+l*nf] = genos[wids[k]+l*n];
      }
      else {
	X[k+(l+q0)*nf] = genos[wids[k]+l*n];
      }
    }
  }

  emmax_log("foo2 - q=%d, X[0]=%lf",q,X[0]);

  // iterate over all possible repeats
  int* hflags = (int*)malloc(sizeof(int)*nf);
  int* hids = (int*)malloc(sizeof(int)*nf);
  int* h2ids = (int*)malloc(sizeof(int)*nf);
  int nh, n2;

  if ( holdf != NULL ) {
    for(j=0, nh=0, n2=0; j < nf; ++j) {
      if ( holdset.find(tfam_headers[wids[j]+n]) != holdset.end() ) {
	h2ids[n2++] = j;
      }
      else {
	hids[nh++] = j;
      }
    }
    n_holdout = n2;
  }

  double* yh = (double*)malloc(sizeof(double)*(nf-n_holdout));
  double* yht = (double*)malloc(sizeof(double)*(nf-n_holdout));
  double* Xh = (double*)malloc(sizeof(double)*(nf-n_holdout)*q);
  double* Xht = (double*)malloc(sizeof(double)*(nf-n_holdout)*q);
  double* Kh = (double*)malloc(sizeof(double)*(nf-n_holdout)*(nf-n_holdout));
  double* eLvalsh = (double*)malloc(sizeof(double)*(nf-n_holdout));
  double* eLvecsh = (double*)malloc(sizeof(double)*(nf-n_holdout)*(nf-n_holdout));
  double* eRvalsh = (double*)malloc(sizeof(double)*(nf-q-n_holdout));
  double* eRvecsh = (double*)malloc(sizeof(double)*(nf-q-n_holdout)*nf);

  double* XDX = (double*)malloc(sizeof(double)*q*q);
  double* XDy = (double*)malloc(sizeof(double)*q);
  double* iXDX = (double*)malloc(sizeof(double)*q*q);
  double* betas = (double*)malloc(sizeof(double)*q);
  double* indbetas;
  if ( indbeta_flag ) {
    indbetas = (double*)malloc(sizeof(double)*ngenos);
  }

  double* XX = (double*)malloc(sizeof(double)*q*q);
  double* iXX = (double*)malloc(sizeof(double)*q*q);
  double* Xy = (double*)malloc(sizeof(double)*q);
  double* betasOLS = (double*)malloc(sizeof(double)*q);
  double* y2cOLS = (double*)malloc(sizeof(double)*n_holdout);

  double* X2 = (double*)malloc(sizeof(double)*n_holdout*q);
  double* y2 = (double*)malloc(sizeof(double)*n_holdout);
  double* K12 = (double*)malloc(sizeof(double)*n_holdout*(nf-n_holdout));
  double* K12t = (double*)malloc(sizeof(double)*n_holdout*(nf-n_holdout));
  double* K22 = (double*)malloc(sizeof(double)*n_holdout*n_holdout);
  double* rht = (double*)malloc(sizeof(double)*(nf-n_holdout));
  double* y2c = (double*)malloc(sizeof(double)*n_holdout);
  double* V22c = (double*)malloc(sizeof(double)*n_holdout*n_holdout);
  
  double *Zh, *Z2;
  if ( indbeta_flag == 1 ) {
    Zh = (double*)malloc(sizeof(double)*(nf-n_holdout)*ngenos);
    Z2 = (double*)malloc(sizeof(double)*n_holdout*ngenos);
  }

  emmax_log("foo3");
  for(i=0; i < n_repeat; ++i) {
    sprintf(buf,"%d.out",i+1);
    outh = open_file_with_suffix(outf,buf,0,1);
    emmax_log("---------------\nRepeating Round %d of out of %d",i+1,n_repeat);
    memset(hflags,0,sizeof(int)*nf);
    if ( holdf == NULL ) {
      for(j=0; j < n_holdout; ++j) {
	while( hflags[(k = (int)floor((double)rand()*(double)nf/(RAND_MAX+1.)))] != 0 );
	hflags[k] = 1;
      }
      
      for(j=0, nh=0, n2=0; j < nf; ++j) {
	if ( hflags[j] == 0 ) {
	  hids[nh++] = j;
	}
	else {
	  h2ids[n2++] = j;
	}
      }
    }
    else {
      // do nothing, and use the same labels
    }

    if ( nh != nf-n_holdout ) {
      emmax_error("nh != nf-n_holdout\n");
    }
    emmax_log("# of holdout = %d, # training = %d",n_holdout,nh);

    // copy training set of yh, X0h, Kh
    for(j=0; j < nh; ++j) {
      yh[j] = y[hids[j]];
      for(k=0; k < q; ++k) {
	Xh[j+k*nh] = X[hids[j]+k*nf];
      }
      for(k=0; k < nh; ++k) {
	Kh[j+k*nh] = K[hids[j]+hids[k]*nf];
      }

      if ( indbeta_flag ) {
	for(k=0; k < ngenos; ++k) {
	  Zh[j+k*nh] = Z[hids[j]+k*nf];
	}
      }
    }

    trSKSh = centered_trace(Kh,nf-n_holdout);
    emmax_log("Eigendecomposition of regressed kinship matrix",nh);
    eigen_R_wo_Z(nh, q, Kh, Xh, eRvalsh, eRvecsh);
    emmax_log("eRvals[0] = %lf, eRvals[nh-q-1] = %lf, eRvals[nh-q] = %lf",eRvalsh[0],eRvalsh[nh-q-1],eRvalsh[nh-q]);

    emmax_log("Estimating REML parameters using %d individuals..",nh);
    REML = REMLE_grid_wo_Z(nh, q, yh, X, Kh, ngrids, llim, ulim, eRvalsh, eRvecsh, &optdelta, &optvg, &optve, &REML0);

    emmax_log("Eigendecomposition of original kinship matrix",nh);
    eigen_L_wo_Z(nh, Kh, eLvalsh, eLvecsh);


    /////////////////
    // Routines for computing betas
    ////////////////
    // beta = (t(X) %*% iV %*% X)^{-1} %*% (t(X) %*% iV %*% y)
    // compute Xt = t(eLvecs)%*%X, yt = t(eLvecs) %*%y
    // Xht = t(eLvecs) %*% Xh
    dgemm(&ct,&cn,&nh,&q,&nh,&onef,eLvecsh,&nh,Xh,&nh,&zerof,Xht,&nh);
    // yht = t(eLvecs) %*% yh
    dgemv(&ct,&nh,&nh,&onef,eLvecsh,&nh,yh,&onen,&zerof,yht,&onen);
    // XDX = t(Xht) %*% diag(1/(delta+eLvals)) %*% Xht
    compute_XDY (ct, cn, Xht, Xht, q, nh, q, eLvalsh, optdelta, XDX, 1., 0.);
    // XDy = t(Xht) %*% diag(1/(delta+eLvals)) %*% yht
    compute_XDY (ct, cn, Xht, yht, q, nh, 1,  eLvalsh, optdelta, XDy, 1., 0.);
    
    // iXDX = solve(XDX)
    matrix_invert( q, XDX, iXDX );

    // betas = iXDX %*% XDy = (t(Xh) %*% solve(V) %*% Xh)^{-1} %*% t(Xh) %*% solve(V) %*% yh
    dgemv(&cn, &q, &q, &onef, iXDX, &q, XDy, &onen, &zerof, betas, &onen);

    // compute indbetas
    if ( indbeta_flag ) {
      double sumy, sumz, sumzy, sumsqz;
      sumy = 0;
      for(j=0; j < nh; ++j) {
	sumy += yh[j];
      }
      for(j=0; j < ngenos; ++j) {
	sumz = sumzy = sumsqz = 0;
	for(k=0; k < nh; ++k) {
	  sumz += Zh[k + j*nh];
	  sumsqz += (Zh[k+j*nh]*Zh[k+j*nh]);
	  sumzy += (yh[j]*Zh[k+j*nh]);
	}
	indbetas[j] = (nh*sumzy-sumz*sumy)/(nh*sumsqz-sumz*sumz);
      }
    }

    // compute pseudo-heritability
    hg = trSKSh/((nf-1.)*optdelta+trSKSh);

    emmax_log("Pseudoheritability = %.5lf",hg);
    emmax_log("REML0 = %.5lf",REML0);
    emmax_log("REML1 = %.5lf",REML);
    emmax_log("REML vg = %.5lf",optvg);
    emmax_log("REML ve = %.5lf",optve);

    // copy X2, y2, K12, K22
    emmax_log("n2 = %d",n2);
    for(j=0; j < n2; ++j) {
      y2[j] = y[h2ids[j]];
      for(k=0; k < q; ++k) {
	X2[j+k*n2] = X[h2ids[j]+k*nf];
      }
      if (indbeta_flag) {
	for(k=0; k < ngenos; ++k) {
	  Z2[j+k*n2] = Z[h2ids[j]+k*nf];
	}
      }
      for(k=0; k < n2; ++k) {
	K22[j+k*n2] = K[h2ids[j]+h2ids[k]*nf];
      }
    }
    
    for(j=0; j < nh; ++j) {
      for(k=0; k < n2; ++k) {
	K12[j+k*nh] = K[hids[j]+h2ids[k]*nf];
      }
    }

    emmax_log("Finished copying y2, X2, K12, K22");

    // OLS part start
    // betas_OLS = (t(X) %*% X)^{-1} %*% t(X)%*%y
    dgemm(&ct, &cn, &q, &q, &nh, &onef, Xh, &nh, Xh, &nh, &zerof, XX, &q);
    dgemv(&ct, &nh, &q, &onef, Xh, &nh, yh, &onen, &zerof, Xy, &onen);
    matrix_invert( q, XX, iXX );
    dgemv(&cn, &q, &q, &onef, iXX, &q, Xy, &onen, &zerof, betasOLS, &onen); // betaOLS = solve(t(X)%*%X)%*%t(X)%*%y
    // y2cOLS = X2 %*% betasOLS
    dgemv(&cn, &n2, &q, &onef, X2, &n2, betasOLS, &onen, &zerof, y2cOLS, &onen);
    if ( indbeta_flag ) {
      dgemv(&cn, &n2, &ngenos, &onef, Z2, &n2, indbetas, &onen, &onef, y2cOLS, &onen);
    }
    // OLS part end

    // rht = -1*Xht*beta+1.0*yht = t(eLvecs)%*%(y-X*beta)
    memcpy(rht, yht, sizeof(double)*nh);
    dgemv(&cn, &nh, &q, &minusonef, Xht, &nh, betas, &onen, &onef, rht, &onen);
    emmax_log("Xht computed, betas[0] = %lf, betasOLS[0] = %lf",betas[0],betasOLS[0]);

    // V12 = optvg*K12+optve*0 = optvg*K12
    // V12t = t(eLvecs)%*%V12 = optvg * t(eLvecs) %*% K12 = optvg * K12t
    dgemm(&ct, &cn, &nh, &n2, &nh, &onef, eLvecsh, &nh, K12, &nh, &zerof, K12t, &nh);

    // y2c = X2%*%beta + t(K12)*V11^{-1}*(y-X%*%beta)
    // y2c = t(K12t)%*%diag(1/(delta+eLvals))%*%rht
    compute_XDY (ct, cn, K12t, rht, n2, nh, 1, eLvalsh, optdelta, y2c, 1., 0.);
    // y2c = y2c + X2%*%beta
    dgemv(&cn, &n2, &q, &onef, X2, &n2, betas, &onen, &onef, y2c, &onen);
    if ( indbeta_flag ) {
      dgemv(&cn, &n2, &ngenos, &onef, Z2, &n2, indbetas, &onen, &onef, y2c, &onen);
    }
    compute_XDY (ct, cn, K12t, K12t, n2, nh, n2, eLvalsh, optdelta, V22c, optvg, 0.);
    for(j=0; j < n2; ++j) {
      for(k=0; k < n2; ++k) {
	V22c[j+k*n2] += optvg*K22[j+k*n2];
      }
      V22c[j+j*n2] += optve;
    }


    // compute cor(y,X*betaOLS)
    double hr, hsum1, hsum2, hsumsq1, hsumsq2, hsum12;
    hsum1 = hsum2 = hsumsq1 = hsumsq2 = hsum12 = 0;
    dgemv(&cn, &nh, &q, &onef, Xh, &nh, betasOLS, &onen, &zerof, rht, &onen);
    for(j=0; j < nh; ++j) {
      hsum1 += yh[j];
      hsumsq1 += (yh[j]*yh[j]);
      hsum2 += rht[j];
      hsumsq2 += (rht[j]*rht[j]);
      hsum12 += (rht[j]*yh[j]);
    }
    hr = (hsum12-hsum1*hsum2/nh)/sqrt((hsumsq1-hsum1*hsum1/nh)*(hsumsq2-hsum2*hsum2/nh));

    double sum1, sum2, sumsq1, sumsq2, sum12, r;
    double sum2OLS, sumsq2OLS, sum12OLS, rOLS;
    sum1 = sum2 = sumsq1 = sumsq2 = sum12 = 0;
    sum2OLS = sumsq2OLS = sum12OLS = 0;
    emmax_log("Round %d Summary\n",i+1);
    emmax_log("-------------------------------------");
    emmax_log("Original\tOLS\tGLS\tSTDEV");
    for(j=0; j < n2; ++j) {
      emmax_log("%-.*lf\t%-.*lf\t%-.*lf\t%-.*lf",ndigits,y[h2ids[j]],ndigits,y2cOLS[j],ndigits,y2c[j],ndigits,sqrt(V22c[j+j*n2]));
      fprintf(outh.fp,"%-.*lf\t",ndigits,y[h2ids[j]]);
      fprintf(outh.fp,"%-.*lf\t",ndigits,y2cOLS[j]);
      fprintf(outh.fp,"%-.*lf\t",ndigits,y2c[j]);
      fprintf(outh.fp,"%-.*lf\n",ndigits,sqrt(V22c[j+j*n2]));

      sum1 += y[h2ids[j]];
      sumsq1 += (y[h2ids[j]]*y[h2ids[j]]);
      sum2 += y2c[j];
      sumsq2 += (y2c[j]*y2c[j]);
      sum2OLS += y2cOLS[j];
      sumsq2OLS += (y2cOLS[j]*y2cOLS[j]);
      sum12 += (y[h2ids[j]]*y2c[j]);
      sum12OLS += (y[h2ids[j]]*y2cOLS[j]);
    }
    emmax_log("-------------------------------------");
    r = (sum12-sum1*sum2/n2)/sqrt((sumsq1-sum1*sum1/n2)*(sumsq2-sum2*sum2/n2));
    rOLS = (sum12OLS-sum1*sum2OLS/n2)/sqrt((sumsq1-sum1*sum1/n2)*(sumsq2OLS-sum2OLS*sum2OLS/n2));
    emmax_log("Correlation:\tOLSPRED=%.5lf (%.5lf)\tGLSPRED=%.5lf (%.5lf)\tOLSFIT=%.5lf (%.5lf)",rOLS,rOLS*rOLS,r,r*r,hr,hr*hr);
    emmax_log("Finished Round %d",i+1);
    fclose(outh.fp);
  }
  /*

  if ( tped_headers != NULL ) free(tped_headers);
  if ( snpid_headers != NULL ) free(snpid_headers);
  if ( covf != NULL ) free(covs_indids);
  if ( snps != NULL ) free(snps);
  if ( lbuf != NULL ) free(lbuf);
  if ( tfam_headers != NULL ) free(tfam_headers);
  if ( phenos_indids!= NULL ) free(phenos_indids);

  if ( mphenoflag != 0 ) free(kins);
  free(K);
  free(Kh);
  free(K12);
  free(K12t);
  free(K22);
  free(V22c);

  free(hflags);
  free(hids);
  free(h2ids);
  free(wids);

  if ( mphenoflag != 0 ) free(phenos);
  free(y);
  free(yh);
  free(yht);
  free(y2c);

  free(rht);
  free(XDX);
  free(iXDX);
  free(XDy);
  free(betas);

  free(covs);
  free(genos);
  free(X);
  free(Xh);
  free(Xht);
  free(X2);

  free(eLvalsh);
  free(eLvecsh);
  free(eRvalsh);
  free(eRvecsh);*/
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
      abort();
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
  double sum1 = 0.;
  double sum2 = 0.;
  for(i=0; i < nq; ++i) {
    sum1 += etas[i]*etas[i]/(lambdas[i]+delta);
    sum2 += log(lambdas[i]+delta);
  }
  return( 0.5*(nq*(log(nq/(2*M_PI))-1.-log(sum1))-sum2) );
}

double LL_REML_wo_Z_inf_delta(int nq, double* lambdas, double* etas) {
  int i;
  double sum1 = 0.;
  for(i=0; i < nq; ++i) {
    sum1 += etas[i]*etas[i];
  }
  return( 0.5*(nq*(log(nq/(2*M_PI))-1.-log(sum1))) );
}

double dLL_REML_wo_Z(int nq, double logdelta, double* lambdas, double* etas) {
  int i;
  double delta = exp(logdelta);
  double sum1 = 0.;
  double sum2 = 0.;
  double sum3 = 0.;
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

  /*
  for(i=0; i < n-q; ++i) {
    emmax_log("etas[%d] = %lf, lambdas[%d] = %lf",i,etas[i],i,evals[i]);
    }*/

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
    //emmax_log("%d\t%lf\t%lf\t%lf\n",i,dLLs[i],dLLs[i+1],LL_REML_wo_Z(n-q,grids[i],evals,etas));
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
    const char* mode = (wflag == 1) ? "wb" : "rb";
    fh.gzfp = gzopen(filename,mode);
    fh.fp = NULL;

    if ( fh.gzfp == NULL ) {
      emmax_error("Cannot open file %s for %s",filename,(wflag == 1) ? "writing" : "reading");
    }
  }
  else {
    const char* mode = (wflag == 1) ? "w" : "r";
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
void fill_XDX_X0 ( double* X0, double* eLvals, double delta, int n, int q0, double* XDX) {
  int i, j, k;
  double tmp;

  for(i=0; i < q0; ++i) {
    for(j=0; j <= i; ++j) {
      tmp = 0;
      for(k=0; k < n; ++k) {
	tmp += ((X0[k+i*n]*X0[k+j*n])/(eLvals[k]+delta));
      }
      XDX[i+j*(q0+1)] = tmp;
      if ( j < i )       
	XDX[j+i*(q0+1)] = tmp;
    }
  }
}

// fill X1-dependent part for computing 
// XDX = t(X) %*% D %*% X
// where X = [X0 x1], dim(X0) = (n,q), dim(x1) = (n,1)
void fill_XDX_X1 ( double* X0, double* x1, double* eLvals, double delta, int n, int q0, double* XDX ) {
  int i, j;
  double tmp;

  for(i=0; i < q0; ++i) {
    tmp = 0;
    for(j=0; j < n; ++j) {
      tmp += ((X0[j+i*n]*x1[j])/(eLvals[j]+delta));
    }
    XDX[i+(q0+1)*q0] = tmp; // fill the last column vector - XDX[i,q0]
    XDX[q0+(q0+1)*i] = tmp; // fill the last row vector - XDX[q0,i]
  }

  tmp = 0;
  for(i=0; i < n; ++i) {
    tmp += (x1[i]*x1[i]/(eLvals[i]+delta));
  }
  XDX[q0*(q0+1)+q0] = tmp; // XDX[q0,q0]
}

// fill X0-dependent part (first q0 elements) for computing 
// XDy = t(X) %*% D %*% y
// where X = [X0 x1], dim(X0) = (n,q), dim(x1) = (n,1), dim(y) = (n,1)
void fill_XDy_X0 ( double* X0, double* y, double* eLvals, double delta, int n, int q0, double* XDy ) {
  int i, j;
  double tmp;

  for(i=0; i < q0; ++i) {
    tmp = 0;
    for(j=0; j < n; ++j) {
      tmp += ((X0[j+n*i]*y[j])/(eLvals[j]+delta));
    }
    XDy[i] = tmp;
  }
}

// fill X1-dependent part (the ladt element) for computing 
// XDX = t(X) %*% D %*% X
// where X = [X0 x1], dim(X0) = (n,q), dim(x1) = (n,1)
void fill_XDy_X1 ( double* x1, double* y, double* eLvals, double delta, int n, int q0, double* XDy ) {
  int i;
  double tmp;
  tmp = 0;
  for(i=0; i < n; ++i) {
    tmp += ((x1[i]*y[i])/(eLvals[i]+delta));
  }
  XDy[q0] = tmp;
}

// C = alpha * cTransA(A) %*% diag(1/(evals+delta)) %*% cTrans(B)
void compute_XDY(char cTransA, char cTransB, double* A, double* B, int dimAA, int dimAB, int dimBB, double* eVals, double delta, double* C, double alpha, double beta) {
  int i, j;

  // computes AD = cTrans(A)%*%diag(1/(evals+delta))
  double* AD = (double*)malloc(sizeof(double) * dimAA * dimAB);
  for(i=0; i < dimAA; ++i) {
    for(j=0; j < dimAB; ++j) {
      AD[i+j*dimAA] = ((cTransA == 'T') ? A[j+i*dimAB] : A[i+j*dimAA]) / (delta + eVals[j]);
    }
  }
  
  // computes C = alpha * AD %*% cTrans(B) + beta * C
  dgemm(&cn,&cTransB,&dimAA,&dimBB,&dimAB,&alpha,AD,&dimAA,B,&dimAB,&beta,C,&dimAA);
  free(AD);
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
  *p_nmiss += nmiss;
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
  fprintf(stderr,"\t-o [out_prefix]  : output file name prefix\n");
  fprintf(stderr,"\t-p [phenof] : 3-column phenotype file with FAMID, INDID at the first two colmns, in the same order of .tfam file. Not required only with -K option\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Likely essential parameters\n");
  fprintf(stderr,"\t-t [tpedf_prefix] : prefix for tped/tfam files\n");
  fprintf(stderr,"\t-s [snp_ids] : list of fixed effect SNP IDs to include\n");
  fprintf(stderr,"\t-n [int] : number to hold out (default - (#inds)/5)\n");
  fprintf(stderr,"\t-r [int] : number of repeats (default : 1)\n");
  fprintf(stderr,"\t-k [kinf] : n * n matrix containing kinship values in the individual order consistent to [tpedf].tfam file. [tpedf].kinf will be used if not specified\n");
  fprintf(stderr,"\t-c [covf] : multi-column covariate file with FAMID, INDID at the first two colmns, in the same order of .tfam file");
  fprintf(stderr,"Optional parameters\n");
  fprintf(stderr,"\t-d [# digits]  : precision of the output values (default : 5)\n");
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
  int *ipiv;
  ipiv=(int*)malloc(n*sizeof(int));
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
  double* work = (double*)malloc(sizeof(double)*lwork);
  /*  Feed this to the lapack inversion routine.  */
  //info = clapack_dgetri (CblasColMajor, n, Y, n, ipiv);
  dgetri(&n,Y,&n,ipiv,work,&lwork,&info);
  
  /*  Cleanup and exit  */
  free(ipiv);
  free(work);
  return info;
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
