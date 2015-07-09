export INTELROOT=/home/hmkang/intel/Compiler/11.1/069
export MKLROOT=${INTELROOT}/mkl
export MKLINCLUDE=${MKLROOT}/include

sh $INTELROOT/bin/iccvars.sh intel64
sh $MKLROOT/tools/environment/mklvars64.sh
export MKLLIBPATH=${MKLROOT}/lib/em64t
export INTELLIBPATH=${INTELROOT}/lib/intel64
export BINPATH=${INTELROOT}/bin/intel64

$BINPATH/icc -Wall -wd1419 -O2 -I $MKLINCLUDE emmax-kin.c -Wl,--start-group $MKLLIBPATH/libmkl_intel_lp64.a $MKLLIBPATH/libmkl_intel_thread.a $MKLLIBPATH/libmkl_core.a $INTELLIBPATH/libguide.a -Wl,--end-group -lpthread -lz -o emmax-kin-intel64

$BINPATH/icc -Wall -wd1419 -wd1418 -wd869 -wd920 -wd981 -wd1572 -O2 -I $MKLINCLUDE emmax.c -Wl,--start-group $MKLLIBPATH/libmkl_intel_lp64.a $MKLLIBPATH/libmkl_intel_thread.a $MKLLIBPATH/libmkl_core.a $INTELLIBPATH/libguide.a -Wl,--end-group -lpthread -lz -o emmax-intel64

$BINPATH/icc -Wall -wd1419 -wd1418 -wd869 -wd920 -wd981 -wd1572 -O2 -I $MKLINCLUDE emmax-predict.cc -Wl,--start-group $MKLLIBPATH/libmkl_intel_lp64.a $MKLLIBPATH/libmkl_intel_thread.a $MKLLIBPATH/libmkl_core.a $INTELLIBPATH/libguide.a -Wl,--end-group -lpthread -lz -o emmax-predict-intel64
