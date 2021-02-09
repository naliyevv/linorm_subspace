include make.inc

DSLSRC = \
    AB13HD.o

all: double

double: $(DSLSRC)

$(DSLSRC): $(FRC)

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $<
