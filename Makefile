CC = clang
CFLAGS = -g -O1 -Wall -Werror -std=gnu11  -fsanitize=address 
LDFLAGS = -pthread

.PHONY: all clean

all: matrix

matrix: main.c matrix.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

clean:
	-rm -f *.o
	-rm -f matrix
	-rm -rf *.dSYM
