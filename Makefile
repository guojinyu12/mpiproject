CFLAGS=-Wall -fopenmp
CC=gcc
LDFLAGS=-fPIC -Wl,-z,noexecstack -I/usr/local/include -L/usr/local/lib -Wl,-rpath,/usr/local/lib
source=$(wildcard *.c)
obj=$(source:.c=.o)
target=matrix
example=$(basename $(source))
$(target):$(obj)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lmpi -lm
main.o:matrix.h
.PHONY: clean
clean:
	$(RM) *.o $(target)
