

%.o: %.cpp
	g++ -c -I $(HOME)/include/ -std=c++17 $^

rectangle: Mesh.o Field.o global.o main.rectangle.o
	g++ -o $@ -std=c++17 $^

3d: Mesh.o Field.o global.o main.3d.o
	g++ -o $@ -std=c++17 $^

exe: Mesh.o Field.o global.o main.o
	g++ -o $@ -std=c++17 $^

exe-cut: Mesh.o Field.o global.o main-cut.o
	g++ -o $@ -std=c++17 $^

exe-tri: Mesh.o Field.o global.o main-tri.o
	g++ -o $@ -std=c++17 $^

exe-sphere-tri: Mesh.o Field.o global.o main-sphere-tri.o
	g++ -o $@ -std=c++17 $^

exe-sphere-quad: Mesh.o Field.o global.o main-sphere-quad.o
	g++ -o $@ -std=c++17 $^

exe-cartesian: Mesh.o Field.o global.o main-cartesian.o
	g++ -o $@ -std=c++17 $^

exe-1.1: Mesh.o Field.o global.o main-1.1.o
	g++ -o $@ -std=c++17 $^

exe-1.3: Mesh.o Field.o global.o main-1.3.o
	g++ -o $@ -std=c++17 $^

exe-1.4: Mesh.o Field.o global.o main-1.4.o
	g++ -o $@ -std=c++17 $^

exe-1.5: Mesh.o Field.o global.o main-1.5.o
	g++ -o $@ -std=c++17 $^

exe-2.2: Mesh.o Field.o global.o main-2.2.o
	g++ -o $@ -std=c++17 $^

exe-2.1: Mesh.o Field.o global.o main-2.1.o
	g++ -o $@ -std=c++17 $^

exe-3.2: Mesh.o Field.o global.o main-3.2.o
	g++ -o $@ -std=c++17 $^

clean:
	rm *.o

run: exe
	./exe

run-cut: exe-cut
	./exe-cut

run-tri: exe-tri
	./exe-tri

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

run-3.2: exe-3.2
	./exe-3.2

PHONY: run run-cut run-tri run-sphere-tri run-sphere-quad run-cartesian run-1.1 run-1.3 run-1.4 run-1.5 run-2.1 run-2.2 run-3.2 clean
