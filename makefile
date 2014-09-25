# makefile for hinet_decon

#FC    = ifort
#FFLAGS = -fast

FC = gfortran
FFLAGS = -O3

all: hinet_decon.x

.SUFFIXES:
.SUFFIXES:.F90 .o .f

.F90.o:
	$(FC) $(FFLAGS) -c $< -o $@

m_std.o   : m_std.F90
m_daytim.o: m_daytim.F90 m_std.o
m_endian.o: m_endian.F90 m_std.o
m_system.o: m_system.F90 m_std.o
m_getopt.o: m_getopt.F90 m_system.o m_std.o
m_rtrend.o: m_rtrend.F90 m_std.o
m_tdecon.o: m_tdecon.F90 m_std.o
m_sac.o   : m_sac.F90    m_daytim.o m_endian.o m_std.o

hinet_decon.x: m_std.o m_daytim.o m_endian.o m_system.o m_getopt.o m_rtrend.o m_tdecon.o m_sac.o hinet_decon.F90
	$(FC) $(FFLAGS) $^ -o $@

clean: 
	rm -f *.o *.mod *~ hinet_decon.x
