all: bds

bds: bds.cpp
	g++ -O3 -I $$(brew --prefix eigen)/include/eigen3/ $< -o $@

clean:
	rm bds
