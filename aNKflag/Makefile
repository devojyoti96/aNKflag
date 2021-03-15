DESTDIR = $(PWD)/
CC=gcc
GSL_INCLUDE_DIR=/home/devojyoti/aNKflag/gsl/include/
GSL_LIBRARIES=-L/home/devojyoti/aNKflag/gsl/lib -Wl,"-R /home/devojyoti/aNKflag/gsl/lib"
CLIBS= $(GSL_LIBRARIES) -lgsl -lgslcblas -lm 
CFLAGS=-I. -I./include -I$(GSL_INCLUDE_DIR) -g -fopenmp
c = $(CC) $(CFLAGS) 
o = -o$(DESTDIR)$@
.c.o: ; $c -c  $< 
.f.o: ; $f -c $< 
.c: ; $c $o $< $(CLIBS)
.f: ; $f $o $< 
make: ankflag
ANKOBJ=ankflg.c exc_fns.c flg_fns.c basic_fns.c aux_fns.c init_pars.c
ankflag:$(ANKOBJ)
	$c $(ANKOBJ) $o $(CLIBS)
clean:
	rm ankflag



