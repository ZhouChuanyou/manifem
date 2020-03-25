

%.o: %.cpp
	g++ -Wshadow -Wall -c -I $(HOME)/include/ -std=c++17 $^
#	g++ -DNDEBUG -c -I $(HOME)/include/ -std=c++17 $^

manifem-exe-%: iterator.o field.o function.o global.o manifold.o mesh.o main-%.o progressive.o
	g++ -o $@ -std=c++17 $^

manifem-exe-9.15: main-9.15.o
	g++ -o $@ -std=c++17 $^

run-%: manifem-exe-%
	./$<

run-test: test
	./$<

test: test.o function.o mesh.o field.o iterator.o manifold.o
	g++ -o $@ -std=c++17 $^

clean:
	rm *.o manifem-exe-*

clean-all: clean
	rm *.eps *.msh

.SECONDARY:

.PHONY: clean clean-all

