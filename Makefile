CFLAGS=-Wall -fopenmp
CC=gcc
LDFLAGS=-fPIC -Wl,-z,noexecstack -I/usr/local/include -L/usr/local/lib -Wl,-rpath,/usr/local/lib
source=$(wildcard *.c)
target=matrix
example=$(basename $(source))
$(target):main.o matric1.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lmpi -lm
main.o:matrix.h
.PHONY: clean
clean:
	$(RM) *.o $(target)
