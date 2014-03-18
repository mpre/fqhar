CFLAGS	= -g -Wall -DDEBUG -O2 -march=native -Wno-deprecated
LDFLAGS = -lz -lm -std=gnu99

all:
	gcc fqhar.c $(LDFLAGS) $(CFLAGS) -o fqhar
