

%.o: %.cpp
	g++ -Wshadow -Wall -c -I $(HOME)/include/ -std=c++17 $^
#	g++ -DNDEBUG -c -I $(HOME)/include/ -std=c++17 $^

rectangle: Mesh.o Field.o global.o main.rectangle.o
	g++ -o $@ -std=c++17 $^

3d: Mesh.o Field.o global.o main.3d.o
	g++ -o $@ -std=c++17 $^

exe: iterator.o field.o function.o global.o main.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-cut: Mesh.o Field.o global.o main-cut.o
	g++ -o $@ -std=c++17 $^

exe-test: function.o test.o
	g++ -o $@ -std=c++17 $^

exe-sphere-tri: Mesh.o Field.o global.o main-sphere-tri.o
	g++ -o $@ -std=c++17 $^

exe-sphere-quad: Mesh.o Field.o global.o main-sphere-quad.o
	g++ -o $@ -std=c++17 $^

exe-cartesian: Mesh.o Field.o global.o main-cartesian.o
	g++ -o $@ -std=c++17 $^

exe-1.1: iterator.o field.o function.o global.o main-1.1.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-1.3: iterator.o field.o function.o global.o main-1.3.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-1.4: iterator.o field.o function.o global.o main-1.4.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-1.5: iterator.o field.o function.o global.o main-1.5.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.1: iterator.o field.o function.o global.o main-2.1.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.2: iterator.o field.o function.o global.o main-2.2.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.3: iterator.o field.o function.o global.o main-2.3.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.4: iterator.o field.o function.o global.o main-2.4.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.5: iterator.o field.o function.o global.o main-2.5.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.6: iterator.o field.o function.o global.o main-2.6.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.7: iterator.o field.o function.o global.o main-2.7.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.8: iterator.o field.o function.o global.o main-2.8.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.9: iterator.o field.o function.o global.o main-2.9.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.10: iterator.o field.o function.o global.o main-2.10.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^


clean:
	rm *.o

run: exe
	./exe

run-cut: exe-cut
	./exe-cut

run-test: exe-test
	./exe-test

run-sphere-tri: exe-sphere-tri
	./exe-sphere-tri

run-sphere-quad: exe-sphere-quad
	./exe-sphere-quad

run-cartesian: exe-cartesian
	./exe-cartesian

run-1.1: exe-1.1
	./exe-1.1

run-1.3: exe-1.3
	./exe-1.3

run-1.4: exe-1.4
	./exe-1.4

run-1.5: exe-1.5
	./exe-1.5

run-2.1: exe-2.1
	./exe-2.1

run-2.2: exe-2.2
	./exe-2.2

run-2.3: exe-2.3
	./exe-2.3

run-2.4: exe-2.4
	./exe-2.4

run-2.5: exe-2.5
	./exe-2.5

run-2.6: exe-2.6
	./exe-2.6

run-2.7: exe-2.7
	./exe-2.7

run-2.8: exe-2.8
	./exe-2.7

run-2.9: exe-2.9
	./exe-2.9

run-2.10: exe-2.10
	./exe-2.10

run-3.2: exe-3.2
	./exe-3.2

PHONY: run run-cut run-test run-sphere-tri run-sphere-quad run-cartesian run-1.1 run-1.3 run-1.4 run-1.5 run-2.1 run-2.2 run-2.3 run-2.4 run-2.5 run-2.6 run-3.2 clean
