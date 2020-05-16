all: bds

bds: bds.cpp
	g++ -I $$(brew --prefix eigen)/include/eigen3/ $< -o $@

