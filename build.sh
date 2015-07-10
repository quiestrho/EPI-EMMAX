export INTELROOT=/opt/intel
export MKLROOT=${INTELROOT}/mkl
export MKLINCLUDE=${MKLROOT}/include

. $INTELROOT/bin/iccvars.sh intel64
. $MKLROOT/bin/mklvars.sh intel64
export MKLLIBPATH=${MKLROOT}/lib/intel64
export INTELLIBPATH=${INTELROOT}/lib/intel64
export BINPATH=${INTELROOT}/bin/intel64

icc -Wall -wd1419 -O2 -I $MKLINCLUDE emmax-kin.c -Wl,--start-group -L $MKLLIBPATH -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -Wl,--end-group -lpthread -lz -o emmax-kin-intel64

icc -Wall -wd1419 -wd1418 -wd869 -wd920 -wd981 -wd1572 -O2 -I $MKLINCLUDE emmax.c -Wl,--start-group $MKLLIBPATH/libmkl_intel_lp64.a $MKLLIBPATH/libmkl_intel_thread.a $MKLLIBPATH/libmkl_core.a -liomp5 -Wl,--end-group -lpthread -lz -o emmax-intel64

icc -Wall -wd1419 -wd1418 -wd869 -wd920 -wd981 -wd1572 -O2 -I $MKLINCLUDE emmax-predict.cc -Wl,--start-group $MKLLIBPATH/libmkl_intel_lp64.a $MKLLIBPATH/libmkl_intel_thread.a $MKLLIBPATH/libmkl_core.a -liomp5 -Wl,--end-group -lpthread -lz -o emmax-predict-intel64
