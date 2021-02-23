EXTRA_CXXFLAGS=
EXTRA_CFLAGS=
CXXFLAGS=-O3 -std=c++17 $(EXTRA_CXXFLAGS)
CFLAGS=-O3 -std=c11 $(EXTRA_CFLAGS)

all: uvgz

clean:
	rm -f uvgz *.o
