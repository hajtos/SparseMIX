
CC=g++
CFLAGS=-I.
DEPS=cluster.h algorithm.h
OBJ=main.o cluster.o algorithm.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

algo: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o *~ algo