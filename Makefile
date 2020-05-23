all: bds

bds: bds.cpp
	g++ -O3 -I $$(brew --prefix eigen)/include/eigen3/ $< -o $@

clean:
	rm bds

test: bds.cpp
	g++ -I $$(brew --prefix eigen)/include/eigen3/ $< -o bds_test
	./bds_test
