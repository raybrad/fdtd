#.SILENT:
SHELL   = /bin/sh
F90     = ifort 
#FOPTS   = -static

#program name
main=fdtd

#module name
mod1=variables
mod2=parameters
mod3=allocmemo
mod4=inputobj
mod5=output
mod6=timeupdate
mod7=source
EXEC=fdtd
VPATH=./variables
##
FFLAGS  = -fpp -openmp -O2  -I./ -fpic #release
FOPTS   = -static
LIBS    =  -lmkl_blas95_lp64 -lmkl_lapack95_lp64  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -openmp -lpthread
##
	##$(cmplr) -O2 -o -fpp  $(main) $(objects)
	#$(cmplr) $(FFLAGS) $(FPPFLAGS) $(FOPTS) \
	#	$(main) $(OBJS) $(LIBS)  -o $(EXEC) 
objects=$(mod1).o $(mod2).o $(mod3).o $(mod4).o $(mod5).o $(mod6).o  $(mod7).o $(main).o 
files= $(mod1).f90 $(mod2).f90 $(mod3).f90 $(mod4).f90 $(mod5).f90 $(mod6).f90 $(mod7).f90 $(main).f90 
cmplr=$(F90) 
$(mod1).o:$(mod1).f90
	$(cmplr) -c $(mod1).f90
$(mod2).o:$(mod2).f90 $(mod1).f90 
	$(cmplr) -c $(mod2).f90
$(mod3).o:$(mod3).f90 $(mod2).f90
	$(cmplr) -c $(mod3).f90
$(mod4).o:$(mod4).f90 $(mod2).f90
	$(cmplr) -c $(mod4).f90
$(mod5).o:$(mod5).f90 $(mod2).f90   
	$(cmplr) -c $(mod5).f90
$(mod6).o:$(mod6).f90 $(mod2).f90
	$(cmplr) -c $(mod6).f90
$(mod7).o:$(mod7).f90 $(mod2).f90
	$(cmplr) -c $(mod7).f90
$(main).o:$(files)
	$(cmplr) -c $(main).f90
$(main):$(objects)
	$(cmplr) $(FFLAGS) $(FOPTS)  -o   $(main) $(objects)$(LIBS)
	rm  *.o  *.mod
clean:
	rm fdtd fdembound fdqmbound 00* sourceE* moleE* job.qsub.o*  Embound.dat

