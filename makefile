#include /soft/apps/phg/intel-2017u4/impi-2017u3/phg-0.9.6/share/phg/Makefile.inc
#include ./Makefile.inc
include /home/yangche/github/phg-0.9.7/Makefile.inc
FC =    mpif90
obj =     main.o
default: all
all: $(obj) main
main: $(obj)
		${LINKER} ${LDFLAGS} -o main $(obj) ${LIBS}
.c.o:
		${CC} ${CFLAGS} ${CPPFLAGS} ${USER_CFLAGS} -c $*.c
%.o : %.f90
		$(FC) $(FLAGS) -c $*.f90

clean:
		rm -rf *.o *.out *.m main main                                