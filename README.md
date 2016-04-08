# EPI-EMMAX
Modified version of EMMAX for GxG interactions and other improvements

# Building EPI-EMMAX (non-intel compiler and MKL version)

if cloning from the git repository on github, you will need GNU automake/autoconf installed, and need to do

autoreconf --force --install

first before anything else. This will generate the configure script and necessary other files. Then, to compile, just do

./configure
make

you can do "sudo make install" if you want, that will probably copy the executibles into /usr/local/bin.

# Building EPI-EMMAX for Intel's compiler and MKL library

I don't know how yet to get the configure script to auto-detect the presence of the intel compiler and MKL libraries in your path.
Please try building using the prepared Makefile.intel, but you may need to modify some of the variables in the file so that the
paths to icc and the MKL libraries are correct for your system. To use Makefile.intel, just do

make -f Makefile.intel

Because of site specific variations on how and where the intel compiler and libraries are installed, and my inability to use them myself,
I can not support this feature. If the compiler seems to work but has compilation errors (not warnings), you can send me the output
for help and I will try to figure it out with you. If it just can't seem to find icc or libmkl, please see your system administrator
for help.

# Differences from EMMAX

This program is very similar to the original EMMAX program by [author], published in [link]. I have modified it to extend the
detail of the output, compute additional statistics, and allow the specification of a specific marker to be included as a
fixed effect in addition to testing every other marker in the genome. More than one such marker can be specified. The linear
model coefficients and -log10(p) will be output for them on every test as well as what I call the "scan marker", the individual
markers tested one by one in the GWAS genome-wide scan. Such markers are independent and additive terms, and please note that
this is a convenience feature mainly, this was always possible with the original EMMAX by manually specifying the marker's
genotype codes as a covariate in the separate covariates file (where you might give PCs as fixed effects, or stratification terms). 
By specifying the marker by identifier on the command line, this saves you the trouble of preparing such files.

Additionally, one may request a GxG interaction term for the covariate marker(s) specified with another command line option. In 
this case, the modified EMMAX program will automatically compute the GxG term for each scan marker and add the term to the model,
allowing the explicit test of quantitative epistatsis. This was not possible with the orignal EMMAX. In this case, the GxG
term's coefficient and -log10(p) will be output as well as the independent additive terms for the covariate marker(s) and the scan
marker. If you specify more than one covariate marker, GxG terms will be added for each of them, but not for each other.

New in this derived version, the program will output the "percent variance explained" (PVE) for the scan marker and any covariate
markers, as well as the PVE of the whole model (including random effects) and the the fixed-effect genetic terms,
(PVE of the model without the random effects part). Interpret this at your own risk - PVE is dependent on the population studied
and not the same as the PVE of the same marker or model in a QTL cross. Recall from basic additive quantitative genetics that the
genetic variance attributable to a genetic term in the model is dependent on the product of the absolute allele effect (assuming
additive effects and indepedence) and the **allele frequency**. Study a different population where the allele is at a different
frequency, and a different PVE will result and neither is incorrect.

This version of EMMAX outputs p-values as -log10(p) rather than the p-value itself. Values will be positive and higher ones indicate
stronger statistical significance. The output contains a column header to indicate what means what. Unlike original EMMAX, regression
coefficients will be given for any covariates specified externally as a covariate file.

The EMMAX predict program/source was never modified or changed. It will not work with the output of this derived version of EMMAX.

# Intel library vs. open source alternatives

The program should work via GNU autoconf/automake for either lapack & atlas libraries, or gslcblas. ./configure will detect
(should) either and automatically add flags to the compilation commands to link to these libraries. The atlas libraries
are supposedly architecture optimized. Having the optimized versions installed is way beyond the scope of this document. The
GSL CBLAS libraries will work but are not specifically optimized for architecture. Please contact your system adiminstrator or
Linux savvy friend for help making sure these libraries are installed. ./configure and building will fail if they are not detected.

If you have access to the Intel compiler and Intel's MKL libraries, this is strongly recommended. Matrix computations will be
multithreaded and optimized for your specific CPU. The program will run much faster. Do not distribute statically linked against
the intel MKL library or compiled by Intel icc without checking to make sure your license with Intel permits this. I am not
responsible for the consequences of any user not respecting Intel's licensing of their compilers and libraries.

# Disclaimer

NO WARRANTY
THE PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, BUT WITHOUT ANY WARRANTY. IT IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW THE AUTHOR WILL BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.


