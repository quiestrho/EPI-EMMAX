#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <unistd.h>
#include <limits.h>
#include <sys/time.h>
#include <zlib.h>

#include "mkl.h"

//#define FIBS_UNIT 1000
#define DEFAULT_NDIGITS 10
#define SZBUF 1024
#define SZBYTE 256
#define N_GENOTYPES 4
#define NA_GENO_CHAR (N_GENOTYPES-1)
#define DEFAULT_ROW_SIZE 100000
#define DEFAULT_SIZE_MATRIX 1000000
#define DEFAULT_SIZE_HEADER 100000
#define DEFAULT_DELIMS " \t\r\n"
#define SZ_LONG_BUF 1000000
#define DEFAULT_TPED_NUM_HEADER_COLS 4
#define DEFAULT_TFAM_NUM_HEADER_COLS 6
#define DEFAULT_TPED_SNPID_INDEX 1
#define DEFAULT_PHENO_NUM_HEADER_COLS 2

struct HFILE {
  int gzflag;       // 1 if gz if used
  int wflag;        // r(0)/w(1) for plain, rb(0)/wb(1) for gz
  int nheadercols;  // # of header columns (0 if nrows=0)
  int nvaluecols;   // # of value cols (0 if nrows=0)
  int nrows;        // # of rows
  FILE* fp;         // plain file handle
  gzFile gzfp;      // gzip file handle
};

// Input routines
void close_file (struct HFILE* fhp);
struct HFILE open_file(char* filename, int gzflag, int wflag);
struct HFILE open_file_with_suffix(char* prefix, char* suffix, int gzflag, int wflag);
//void read_matrix_with_col_headers( struct HFILE* fhp, int nheadercols, char* delims, int symmetric, int* p_nmiss, unsigned char** matrix, char*** headers);
void read_matrix_with_col_headers( struct HFILE* fhp, int nheadercols, char* delims, int* p_nmiss, unsigned char** matrix, char*** headers);
unsigned char* tokenize_tped_line_with_col_headers( struct HFILE* fhp, int nheadercols, char* delims, char* lbuf, unsigned char* values, char** headers, int* p_nvalues, int* p_nmiss );

void emmax_error( const char* format, ... );
void print_help(void);
FILE* readfile(char* filename);

void print_help(void) {
  fprintf(stderr,"Usage: emmax_kin [tpedf]\n");
  fprintf(stderr,"Required parameters\n");
  fprintf(stderr,"\t[tpedf]     : tped file\n");
  fprintf(stderr,"Optional parameters\n");
  fprintf(stderr,"\t-d [# digits]  : precision of the kinship values (default : 10)\n");
  fprintf(stderr,"\t-M [float] : maximum memory in GB (default: 4.0)\n");
  fprintf(stderr,"\t-s : compute IBS kinship matrix (default is Balding-Nicholas)\n");
  fprintf(stderr,"\t-v : turn on verbose mode\n");
  fprintf(stderr,"\t-r : randomly fill missing genotypes (default is imputation by average)\n");
  fprintf(stderr,"\t-x : include non-autosomal chromosomes in computing kinship matrices\n");
  fprintf(stderr,"\t-S [int] : set random seed\n");
  fprintf(stderr,"\t-m [float] : MAF threshold (default is 0)\n");
  fprintf(stderr,"\t-c [float] : Call rate threshold (default is 0)\n");
}

//void read_matrix_with_col_headers( struct HFILE* fhp, int nheadercols, char* delims, int symmetric, int* p_nmiss, unsigned char** matrix, char*** headers) {
void read_matrix_with_col_headers( struct HFILE* fhp, int nheadercols, char* delims, int* p_nmiss, unsigned char** matrix, char*** headers) {
  char* lbuf = (char*) malloc(sizeof(char*) * SZ_LONG_BUF);
  int szmat = DEFAULT_SIZE_MATRIX;
  int szheader = DEFAULT_SIZE_HEADER;
  unsigned char* cmat = (unsigned char*) malloc(sizeof(unsigned char) * szmat );
  char** cheaders = (char**) malloc(sizeof(char*) * szheader );
  int nvalues, i, j, nmiss;

  fhp->nheadercols = nheadercols; 
  nmiss = 0;

  while( tokenize_tped_line_with_col_headers(fhp, nheadercols, delims, lbuf, &cmat[fhp->nrows*fhp->nvaluecols], &cheaders[fhp->nrows*fhp->nheadercols], &nvalues, &nmiss) != NULL ) {
    if ( fhp->nrows == 1 ) {
      fhp->nvaluecols = nvalues;
    }
    else if ( fhp->nvaluecols != nvalues ) {
      emmax_error("The column size %d do not match to %d at line %d\n",nvalues,fhp->nvaluecols,fhp->nrows);
    }

    if ( (fhp->nrows+1)*(fhp->nvaluecols) > szheader ) {
      szheader *= 2;
      fprintf(stderr,"Header size is doubled to %d\n",szheader);
      cheaders = (char**) realloc( cheaders, sizeof(char*) * szheader );
    }

    if ( (fhp->nrows+1)*(fhp->nvaluecols) > szheader ) {
      szmat *= 2;
      fprintf(stderr,"Matrix size is doubled to %d\n",szmat);
      cmat = (unsigned char*) realloc( cmat, sizeof(unsigned char) * szmat );
    }
  }
  free(lbuf);

  *p_nmiss = nmiss;
  
  unsigned char* fmat = (unsigned char*) malloc(sizeof(unsigned char)*fhp->nrows*fhp->nvaluecols);
  char** fheaders = (char**) malloc(sizeof(char*)*fhp->nrows*fhp->nheadercols);
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

unsigned char* tokenize_tped_line_with_col_headers( struct HFILE* fhp, int nheadercols, char* delims, char* lbuf, unsigned char* values, char** headers, int* p_nvalues, int* p_nmiss ) {
  int j;
  char *token;
  unsigned char ctoken;

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
    else {
      ctoken = (unsigned char)(token[0]-'0');
      
      if ( ctoken > 2 ) {
	fprintf(stderr,"Unrecognized token %s\n",token);
	abort();
      }
      
      if ( (j-nheadercols) % 2 == 0 ) {
	values[(j-nheadercols)/2] = ctoken;
      }
      else {
	if ( ( ctoken > 0 ) && ( values[(j-nheadercols)/2] == 0 ) ) {
	  fprintf(stderr,"Unmatched token pair 0 %s\n",token);
	  abort();
	}
	else if ( ( ctoken == 0 ) && ( values[(j-nheadercols)/2] > 0 ) ) {
	  fprintf(stderr,"Unmatched token pair - %d 0\n",(int)values[(j-nheadercols)/2]);
	  abort();
	}
	values[(j-nheadercols)/2] = (unsigned char)(values[(j-nheadercols)/2]+ctoken);
      }
    }
    token = strtok(NULL, delims);
  }
  //fprintf(stderr,"tokenize-line ended %d %d\n",j,nheadercols);
  if ( (j-nheadercols) % 2 != 0 ) {
    fprintf(stderr,"Number of value tokens are not even %d\n",j-nheadercols);
    abort();
  }

  *p_nvalues = (j-nheadercols)/2;
  *p_nmiss += nmiss;
  ++(fhp->nrows);

  if ( j < nheadercols ) {
    fprintf(stderr,"Number of header columns are %d, but only %d columns were observed\n", nheadercols, j);
    abort();
  }

  return values;
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
      emmax_error("Cannot open file %s for reading",filename);
    }
  }
  else {
    char* mode = (wflag == 1) ? "w" : "r";
    fh.gzfp = (gzFile) NULL;
    fh.fp = fopen(filename,mode);

    if ( fh.fp == NULL ) {
      emmax_error("Cannot open file %s for writing",filename);
    }
  }
  return fh;
}

void emmax_error( const char* format, ... ) {
  va_list args;
  fprintf(stderr, "ERROR: ");
  va_start (args, format);
  vfprintf(stderr, format, args);
  va_end (args);
  fprintf(stderr,"\n");
  abort();
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

int main(int argc, char** argv) {
  int i, j, n, c, ac0, ac1, ac2, nmiss, nelems, nex, n_sum_nin, nin, nex_maf, nex_call_rate, nex_autosomal;
  int verbose, ndigits, tped_nheadercols, tfam_nheadercols, flag_autosomal, rand_fill_flag, ibs_flag, n_unit_lines;
  unsigned char *snprow;
  char *suffix, buf[SZBUF];
  double f, max_memory_GB, maf_thres, call_rate_thres, aaf, call_rate;
  //long *fibs_sums, *scores, mean_score;
  double *kin, *snpunit;
  char *tpedf, *delims, *lbuf;
  char **tfam_headers, **tped_headers;
  struct HFILE tpedh, tfamh, kinsh;
  struct timeval tv;

  // set default params
  gettimeofday(&tv, NULL);
  srand((unsigned int)tv.tv_usec);
  delims = DEFAULT_DELIMS;
  tped_nheadercols = DEFAULT_TPED_NUM_HEADER_COLS;
  tfam_nheadercols = DEFAULT_TFAM_NUM_HEADER_COLS;
  tped_headers = tfam_headers = NULL;
  tpedf = lbuf = 0;
  flag_autosomal = 1;
  rand_fill_flag = 0;
  ibs_flag = 0;
  verbose = 0;
  ndigits = DEFAULT_NDIGITS;
  max_memory_GB = 4.0; // 4.0GB
  maf_thres = 0.0;
  call_rate_thres = 0.0;

  // read arguments and update params
  while ((c = getopt(argc, argv, "d:rsc:vxS:M:m:c:")) != -1 ) {
    switch(c) {
    case 'd': // precision of digits
      ndigits = atoi(optarg);
      break;
    case 'r':
      rand_fill_flag = 1;
      break;
    case 's':
      ibs_flag = 1;
      break;
    case 'v':
      verbose = 1;
      break;
    case 'x':
      flag_autosomal = 0;
      break;
    case 'S':
      srand(atoi(optarg));
      break;
    case 'M':
      max_memory_GB = atof(optarg);
      break;
    case 'm':
      maf_thres = atof(optarg);
      break;
    case 'c':
      call_rate_thres = atof(optarg);
      break;
    default:
      fprintf(stderr,"Error : Unknown option unsigned character %c\n",c);
      abort();
    }
  }

  // Sanity check for the number of required parameters
  if ( argc != optind + 1 ) {
    print_help();
    abort();
  }

  // Read required parameters
  tpedf = argv[optind++];

  if ( verbose) fprintf(stderr,"\nReading TFAM file %s.tfam ....\n",tpedf);

  tfamh = open_file_with_suffix(tpedf, "tfam", 0, 0);
  //read_matrix_with_col_headers( &tfamh, tfam_nheadercols, delims, 0, &nmiss, NULL, &tfam_headers);
  read_matrix_with_col_headers( &tfamh, tfam_nheadercols, delims, &nmiss, NULL, &tfam_headers);
  n = tfamh.nrows; // n is # of individuals
  if ( verbose ) fprintf(stderr,"Identified %d individuals from TFAM file\n",n);

  // compute the # of lines to read together
  // we would need two n*n matrix, and n*m matrix
  // (n*n + n*m)*sizeof(double) < M*1e9 
  // m < M*1e9/sizeof(double)/n - n
  n_unit_lines = (int)floor((max_memory_GB * 1.0e9 / sizeof(double) / n - n)/2)*2;
  n_sum_nin = 0;
  char cn = 'N', ct = 'T';
  double one = 1.;

  if ( verbose ) fprintf(stderr,"Setting # unit lines = %d to fit the memory requirement\n",n_unit_lines);

  snprow = (unsigned char*)malloc(sizeof(unsigned char)*n);
  tped_headers = (char**)malloc(sizeof(char*)*n);
  lbuf = (char*) malloc(sizeof(char*) * SZ_LONG_BUF);

  kin = (double*)calloc(n*n, sizeof(double));
  snpunit = (double*)malloc(n*n_unit_lines * sizeof(double));

  if ( verbose) fprintf(stderr,"Reading TPED file %s.tped ....\n",tpedf);

  tpedh = open_file_with_suffix(tpedf, "tped", 0, 0);
  tpedh.nheadercols = tped_nheadercols;

  nex_autosomal = nex_maf = nex_call_rate = 0;

  for ( i=0, nin=0, nex = 0; tokenize_tped_line_with_col_headers( &tpedh, tped_nheadercols, delims, lbuf, snprow, tped_headers, &nelems, &nmiss) != NULL; ++i) {
    if ( ( verbose ) && ( i % 10000 ) == 0 ) fprintf(stderr,"Reading %d SNPs\n",i);

    if ( ( flag_autosomal == 1 ) && ( ( atoi(tped_headers[0]) == 0 ) || ( atoi(tped_headers[0]) > 22 ) ) ) // if SNP is not in autosomal chromosomes
      //if ( ( flag_autosomal == 1 ) && ( atoi(tped_headers[0]) > 22 ) ) // if SNP is not in autosomal chromosomes
    {
      ++nex; // # excluded snps from the last unit
      ++nex_autosomal;
      continue;
    }

    if ( nelems != n ) {
      emmax_error("Number of values %d in line %d do not match to %d, the number of columns\n", nelems, tpedh.nvaluecols, n);
    }

  /*
    Perform rapid kinship generate IBS, BN, NCOR matrix
--------------------
IBS pairwise matrix
--------------------
i/j  0  1  2  3
0    NA NA NA NA
1    NA 2  1  0 
2    NA 1  2  1
3    NA 0  1  2

In fact is what it does is
X = (m*n) genotype matrix (1,2,3 coded)
Xn = X-2
K = (t(Xn) %*% Xn)/(2*m)
------------------
* NA column may be just averaged or predicted based on r2 with previous SNP
  with a certain window size

---------------------
pairwise BN matrix
--------------------
X : (m*n) genotype matrix
Xn : (m*n) matrix each row standardized, missing assigned to 0
K = t(Xn) %*% Xn / L
--------------------
*/

    ac0 = ac1 = ac2 = 0;
    for(j=0; j < n; ++j) {
      if ( snprow[j] > 0 ) {
	if ( snprow[j] < 2 ) {
	  fprintf(stderr,"Invalid snprow[%d] value %d at line %d, individual %d\n",j,(int)snprow[j],i,j);
	}
	snprow[j] -= 2; // from 2,3,4 to 0,1,2 coding

	switch(snprow[j]) {
	case 0:
	  ++ac0;
	  break;
	case 1:
	  ++ac1;
	  break;
	case 2:
	  ++ac2; 
	  break;
	default:
	  emmax_error("Unknown allele %s, converted to %d\n",buf,(int)snprow[j]);
	  break;
	}
      }
      else {
	snprow[j] = (unsigned char)NA_GENO_CHAR;
      }
    }

    call_rate = (double)(ac0+ac1+ac2)/(double)n;
    //fprintf(stderr,"CallRate = %lf\n", call_rate);
    if ( call_rate <= call_rate_thres ) {
      ++nex;
      ++nex_call_rate;
      continue;
    }
    aaf = (double)(ac1+2*ac2)/(double)(2*(ac0+ac1+ac2));
    if ( ( aaf <= maf_thres ) || ( 1.-aaf <= maf_thres ) ) {
      ++nex;
      ++nex_maf;
      continue;
    }
    
    if ( rand_fill_flag == 1 ) {
      for(j=0; j < n; ++j) {
	if ( snprow[j] == (unsigned char)NA_GENO_CHAR ) {
	  if ( (rand() / (double) RAND_MAX) > aaf ) {
	    if ( (rand() / (double) RAND_MAX) > aaf ) {
	      snprow[j] = (unsigned char)2;
	      //ac1 += 2;
	      ++ac2;
	    }
	    else {
	      snprow[j] = (unsigned char)1;
	      ++ac1;
	      //++ac0;
	      //++ac1;
	    }
	  }
	  else {
	    if ( (rand() / (double) RAND_MAX) > aaf ) {
	      snprow[j] = (unsigned char)1;
	      ++ac1;
	      //++ac0;
	      //++ac1;
	    }
	    else {
	      snprow[j] = (unsigned char)0;
	      ++ac0;
	      //ac0 += 2;
	    }
	  }
	}
      }
      aaf = (double)(ac1+2*ac2)/(double)(2*(ac0+ac1+ac2));
      //aaf = (double)ac1/(double)(ac0+ac1);
    }

    // copy current values to arrays
    // Xn = [-1,1] - two rows per SNP : IBS matrix
    // Xn ~ N(0,1) : BN matrix
    if ( ibs_flag == 1 ) {
      for(j=0; j < n; ++j) {
	if ( snprow[j] == (unsigned char)NA_GENO_CHAR ) {
	  snpunit[nin+  j*n_unit_lines] = 2.*aaf-1.;
	  snpunit[nin+1+j*n_unit_lines] = 2.*aaf-1.;
	}
	else {
	  if ( snprow[j] == 0 ) {
	    snpunit[nin+  j*n_unit_lines] = -1.;
	    snpunit[nin+1+j*n_unit_lines] = -1.;
	  }
	  else if ( snprow[j] == 1 ) {
	    snpunit[nin+  j*n_unit_lines] = 1.;
	    snpunit[nin+1+j*n_unit_lines] = -1.;
	  }
	  else if ( snprow[j] == 2 ) {
	    snpunit[nin+  j*n_unit_lines] = 1.;
	    snpunit[nin+1+j*n_unit_lines] = 1.;
	  }
	  else {
	    emmax_error("Invalid genotype %d\n",snprow[j]);
	  }
	}
      }
      nin += 2;
    }
    else {
      for(j=0; j < n; ++j) {
	if ( snprow[j] == (unsigned char)NA_GENO_CHAR ) {
	    snpunit[nin+j*n_unit_lines] = 0.;
	}
	else {
	  snpunit[nin+j*n_unit_lines] = ((double)snprow[j]-(aaf*2.))/sqrt(4*aaf*(1-aaf));
	}
      }
      ++nin;
    }
    
    // check if nin == n_unit_lines 
    if ( nin >= n_unit_lines ) {
      if ( verbose ) {
	fprintf(stderr,"At SNP %d, intermediately computing kinship matrix with %d SNPs..\n",i,n_unit_lines);
      }
      // kin = kin + t(snpunit)%*%(snpunit)
      dgemm(&ct,&cn,&n,&n,&n_unit_lines,&one,snpunit,&n_unit_lines,snpunit,&n_unit_lines,&one,kin,&n);
      memset(snpunit, 0, sizeof(double)*n_unit_lines*n);

      n_sum_nin += nin;
      nin = 0;
    }
  }
  close_file(&tpedh);

  if ( verbose ) {
    fprintf(stderr,"Succesfully finished reading TPED file\n");
  }
  if ( nin > 0 ) {
    if ( verbose ) {
      fprintf(stderr,"Computing kinship matrix with the remaining %d SNPs..\n",nin);
    }
    dgemm(&ct,&cn,&n,&n,&n_unit_lines,&one,snpunit,&n_unit_lines,snpunit,&n_unit_lines,&one,kin,&n);
    n_sum_nin += nin;
  }

  if ( ibs_flag == 1 ) {
    if ( rand_fill_flag == 1 ) {
      suffix = "rIBS.kinf";
    }
    else {
      suffix = "aIBS.kinf";
    }
  }
  else {
    if ( rand_fill_flag == 1 ) {
      suffix = "rBN.kinf";
    }
    else {
      suffix = "aBN.kinf";
    }
  }
  kinsh = open_file_with_suffix( tpedf, suffix, 0, 1 );


  if ( verbose ) fprintf(stderr,"Printing the kinship matrix to file %s.%s\n",tpedf,suffix);

  for(i=0; i < n; ++i) {
    for(j=0; j < n; ++j) {
      if ( j > 0 ) fprintf(kinsh.fp,"\t");
      f = (double)kin[i+j*n]/(double)(n_sum_nin);
      if ( ibs_flag == 1 ) {
	fprintf(kinsh.fp,"%-.*lf",ndigits,0.5*f+0.5);
      }
      else {
	fprintf(kinsh.fp,"%-.*lf",ndigits,f);  
      }
    }
    fprintf(kinsh.fp,"\n");
  }
  close_file(&kinsh);
  free(snprow);
  free(kin);
  free(snpunit);
  free(lbuf);
  free(tped_headers);
  free(tfam_headers);
  return 0;
}
