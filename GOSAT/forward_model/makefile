FILES = commonData.f90 absc_kg_gen.f90 forward_calc.f90 create_emission.f90
#FILES = forward_calc_abscErr.f90 inverse_abscErr_1.f90
#FILES = commonData_toddtest.f90 absc_kg_gen_toddtest.f90 getabsco.f90

EXE = emission

F95 = ifort 
F95FLAGS = -FR -w95
F95FLAGS1 = -FR -w95 
F95FLAGS5 = -FR -w -cm -CB -traceback -g -assume byterecl

#FFLAGS = -r8 -01 -Vaxlib

LDFLAGS =-r8 -assume byterecl

OBJS  = commonData.o absc_kg_gen.o forward_calc.o create_emission.o
#OBJS  = forward_calc_abscErr.o inverse_abscErr_1.o
#OBJS  = commonData_toddtest.o absc_kg_gen_toddtest.o getabsco.o


all: $(EXE)
#$(EXE): $(OBJS) 
#	$(F95) $(F95FLAGS5) $(LDFLAGS) -o $(EXE) $(OBJS)

$(EXE): $(OBJS)
	$(F95) $(LDFLAGS) -o $(EXE) $(OBJS) -L /public/home/lcc-rt/opt/fftw/.libs/ -lfftw3 -lm

%.o: %.f90
	$(F95) $(F95FLAGS5) -c $*.f90

clean:
	-rm -f $(OBJS) *.mod *.o

#inverse_abscErr_1.o: forward_calc_abscErr.o
#testsinc1.o: absc_kg_gen_toddtest.o commonData_toddtest.o quadlib.o
absc_kg_gen_toddtest.o: commonData.o
forward_calc_abscErr.o: commonData.o
createabscoCO2.o: commonData.o absc_kg_gen.o forward_calc.o
#getabsco.o: absc_kg_gen_toddtest.o
