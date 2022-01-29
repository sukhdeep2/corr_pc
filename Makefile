compiler= mpic++

CFLAGS=-c #-I/opt/gsl/
LDFLAGS=-lgsl -lgslcblas -lgomp -fopenmp
debug= -g2 -gdwarf-2
optimize=  -O2

cc=$(compiler) $(optimize)  $(debug)

Targets=main.o corels.o corr.o do_corr.o wl_corr.o bins_calcs.o read_dat.o initialization.o outp.o calcs.o jk.o sky_calcs.o PB_calcs.o

corr_pc	: $(Targets)
		$(cc) $(Targets) $(LDFLAGS) -o corr_pc
		chmod +X corr_pc
clean:
	rm $(Targets)
	@echo Done

main.o: main.cpp data_def.h corels.h read_dat.h initialization.h outp.h bins_calcs.h do_corr.h corr.h
	$(cc) $(CFLAGS) $(LDFLAGS) main.cpp

do_corr.o: do_corr.cpp data_def.h corels.h read_dat.h initialization.h outp.h bins_calcs.h corr.h do_corr.h
	$(cc) $(CFLAGS) $(LDFLAGS) do_corr.cpp

wl_corr.o: wl_corr.cpp data_def.h corels.h read_dat.h initialization.h outp.h bins_calcs.h corr.h do_corr.h
	$(cc) $(CFLAGS) $(LDFLAGS) wl_corr.cpp

bins_calcs.o: bins_calcs.cpp data_def.h calcs.h bins_calcs.h
	$(cc) $(CFLAGS) $(LDFLAGS) bins_calcs.cpp

corr.o:	corr.cpp data_def.h corels.h calcs.h corr.h
	$(cc) $(CFLAGS) $(LDFLAGS) corr.cpp

corels.o: corels.cpp data_def.h corels.h calcs.h
	$(cc) $(CFLAGS) $(LDFLAGS) corels.cpp

read_dat.o: read_dat.cpp read_dat.h data_def.h calcs.h
	$(cc) $(CFLAGS) $(LDFLAGS) read_dat.cpp

initialization.o: initialization.cpp initialization.h data_def.h calcs.h
	$(cc) $(CFLAGS) $(LDFLAGS) initialization.cpp

calcs.o: calcs.cpp data_def.h calcs.h
	$(cc) $(CFLAGS) $(LDFLAGS) calcs.cpp

sky_calcs.o:	sky_calcs.cpp data_def.h sky_calcs.h
	$(cc) $(CFLAGS) $(LDFLAGS) sky_calcs.cpp

PB_calcs.o:	PB_calcs.cpp data_def.h PB_calcs.h
	$(cc) $(CFLAGS) $(LDFLAGS) PB_calcs.cpp

jk.o:	jk.cpp data_def.h jk.h
		$(cc) $(CFLAGS) $(LDFLAGS) jk.cpp

outp.o:	outp.cpp data_def.h outp.h
	$(cc) $(CFLAGS) $(LDFLAGS) outp.cpp
