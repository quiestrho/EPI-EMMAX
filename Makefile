
export INTELROOT=/opt/intel
export MKLROOT=$(INTELROOT)/mkl
export MKLINCLUDE=$(MKLROOT)/include

export MKLLIBPATH=$(MKLROOT)/lib/intel64
export INTELLIBPATH=$(INTELROOT)/lib/intel64
export BINPATH=$(INTELROOT)/bin/intel64

ICC=/opt/intel/bin/icc #. $(INTELROOT)/bin/iccvars.sh intel64 && . $(MKLROOT)/bin/mklvars.sh intel64 && icc

all: emmax-kin-intel64 emmax-intel64 emmax-predict-intel64

emmax-kin-intel64: emmax-kin.c
	$(ICC) -Wall -wd1419 -O2 -I $(MKLINCLUDE) $< -Wl,--start-group -L $(MKLLIBPATH) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -Wl,--end-group -lpthread -lz -o $@

emmax-intel64: emmax.c
	$(ICC) -Wall -wd1419 -wd1418 -wd869 -wd920 -wd981 -wd1572 -O2 -I $(MKLINCLUDE) $< -Wl,--start-group $(MKLLIBPATH)/libmkl_intel_lp64.a $(MKLLIBPATH)/libmkl_intel_thread.a $(MKLLIBPATH)/libmkl_core.a -liomp5 -Wl,--end-group -lpthread -lz -o $@

emmax-predict-intel64: emmax-predict.cc
	$(ICC) -Wall -wd1419 -wd1418 -wd869 -wd920 -wd981 -wd1572 -O2 -I $(MKLINCLUDE) $< -Wl,--start-group $(MKLLIBPATH)/libmkl_intel_lp64.a $(MKLLIBPATH)/libmkl_intel_thread.a $(MKLLIBPATH)/libmkl_core.a -liomp5 -Wl,--end-group -lpthread -lz -o $@


