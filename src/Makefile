

%.o: %.cpp
	g++ -Wshadow -Wall -c -I $(HOME)/include/ -std=c++17 $^
#	g++ -DNDEBUG -c -I $(HOME)/include/ -std=c++17 $^

manifem-exe-%: iterator.o field.o function.o global.o manifold.o mesh.o main-%.o progressive.o
	g++ -o $@ -std=c++17 $^

manifem-exe-8.15: main-8.15.o
	g++ -o $@ -std=c++17 $^

run-%: manifem-exe-%
	./$<

clean:
	rm *.o manifem-exe-*

clean-all: clean
	rm *.eps *.msh

.SECONDARY:

.PHONY: clean clean-all run-%

