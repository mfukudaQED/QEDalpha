FC      = ifort
F90     = ifort
#FFLAGS = -O2
#FFLAGS = -O2 -ffree-line-length-huge
# for contracted basis set
OBJS   = module.o simQED.o sub_exvalphys.o sub_eigsym.o sub_diffsolver.o sub_derivative.o sub_makeQmat.o sub_int.o sub_density2.o sub_gauss_int.o sub_arrangeEig.o sub_readBasis.o sub_readVectors.o sub_readDiracOutput.o relvnorm.o sub_Arad.o sub_twoele_divided.o 

sim: $(OBJS)
	$(F90) $(FFLAGS) $(OBJS) -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $*.f

%.o: %.F
	$(FC) $(FFLAGS) -c $*.F

%.o: %.F90
	$(F90) $(FFLAGS) $(INCS) -c -o $@ $<

%.o: %.f90
	$(F90) $(FFLAGS) $(INCS) -c -o $@ $<

clean:
	$(RM) *.o *.mod *.log *~



