

solve: solve.cpp bds.cpp 
	g++ -std=c++14 -O3 -I $$(brew --prefix eigen)/include/eigen3/ $^ -o $@

test: test.cpp bds.cpp
	g++ -std=c++14 -O3 -I $$(brew --prefix eigen)/include/eigen3/ $^ -o $@
	./test


.PHONY: clean test

clean:
	rm -f solve test

