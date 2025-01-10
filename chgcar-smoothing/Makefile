CC:=g++
INCLUDE:=-I.
LIBS:=-L.
LIBM:=-lm
LIBGSL:=-lgsl -lgslcblas

INCLUDES = *.H 

TARGET = chgcar_smoothing

all: $(TARGET)

chgcar_smoothing.o: chgcar_smoothing.C $(INCLUDES)
	$(CC) -Wall -c chgcar_smoothing.C

chgcar_smoothing: chgcar_smoothing.o
	$(CC) -o chgcar_smoothing chgcar_smoothing.o $(LIBGSL) $(LIBM)

clean:
	rm -f *.o chgcar_smoothing