CFLAGS=-Wall -fopenmp
CC=gcc
LDFLAGS=-fPIC -Wl,-z,noexecstack -I/usr/local/include -L/usr/local/lib -Wl,-rpath,/usr/local/lib
source=$(wildcard *.c)
example=$(basename $(source))
$(example):$(example).o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lpthread -lmpi
.PHONY: clean
clean:
	$(RM) *.o $(example)
