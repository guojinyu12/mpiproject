CFLAGS=-Wall -fopenmp
CC=gcc
LDFLAGS=-fPIC -Wl,-z,noexecstack -I/usr/local/include -L/usr/local/lib -Wl,-rpath,/usr/local/lib
source=$(wildcard *.c)
obj=$(source:.c=.o)
target=matrix
example=$(basename $(source))
all:$(target) series
$(target):main.o matrix1.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lmpi -lm
main.o:matrix.h
series: series.o matrix1.o
	$(CC) $(CFLAGS) -o $@ $^
.PHONY: clean all
clean:
	$(RM) *.o $(target) series
