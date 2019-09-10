

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

exe-sphere: Mesh.o Field.o global.o main-sphere.o
	g++ -o $@ -std=c++17 $^

exe-cartesian: Mesh.o Field.o global.o main-cartesian.o
	g++ -o $@ -std=c++17 $^

exe-1.1: Mesh.o Field.o global.o main-1.1.o
	g++ -o $@ -std=c++17 $^

exe-1.3: Mesh.o Field.o global.o main-1.3.o
	g++ -o $@ -std=c++17 $^

exe-1.4: Mesh.o Field.o global.o main-1.4.o
	g++ -o $@ -std=c++17 $^

exe-2.2: Mesh.o Field.o global.o main-2.2.o
	g++ -o $@ -std=c++17 $^

clean:
	rm *.o

run: exe
	./exe

run-cut: exe-cut
	./exe-cut

run-sphere: exe-sphere
	./exe-sphere

run-cartesian: exe-cartesian
	./exe-cartesian

run-1.1: exe-1.1
	./exe-1.1

run-1.3: exe-1.3
	./exe-1.3

run-1.4: exe-1.4
	./exe-1.4

run-2.2: exe-2.2
	./exe-2.2

PHONY: run run-cut run-sphere run-cartesian run-1.1 run-1.3 run-1.4 run-2.2 clean
