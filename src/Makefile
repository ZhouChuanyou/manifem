

%.o: %.cpp
	g++ -Wshadow -Wall -c -I $(HOME)/include/ -std=c++17 $^
#	g++ -DNDEBUG -c -I $(HOME)/include/ -std=c++17 $^

exe-%: iterator.o field.o function.o global.o manifold.o mesh.o main-%.o progressive.o
	g++ -o $@ -std=c++17 $^

clean:
	rm *.o

run-%: exe-%
	./$<

test: metric-tree.o
	g++ -o $@ -std=c++17 $^

run-test: test
	./test

.SECONDARY:

.PHONY: run-% clean
