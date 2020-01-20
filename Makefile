

%.o: %.cpp
	g++ -Wshadow -Wall -c -I $(HOME)/include/ -std=c++17 $^
#	g++ -DNDEBUG -c -I $(HOME)/include/ -std=c++17 $^


exe-test: iterator.o field.o function.o global.o manifold.o mesh.o progressive.o test.o
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

exe-2.11: iterator.o field.o function.o global.o main-2.11.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.12: iterator.o field.o function.o global.o main-2.12.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.13: iterator.o field.o function.o global.o main-2.13.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.14: iterator.o field.o function.o global.o main-2.14.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-2.15: iterator.o field.o function.o global.o main-2.15.o manifold.o mesh.o
	g++ -o $@ -std=c++17 $^

exe-3.1: iterator.o field.o function.o global.o main-3.1.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.2: iterator.o field.o function.o global.o main-3.2.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.3: iterator.o field.o function.o global.o main-3.3.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.4: iterator.o field.o function.o global.o main-3.4.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.5: iterator.o field.o function.o global.o main-3.5.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.6: iterator.o field.o function.o global.o main-3.6.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.7: iterator.o field.o function.o global.o main-3.7.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.8: iterator.o field.o function.o global.o main-3.8.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.14: iterator.o field.o function.o global.o main-3.14.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.15: iterator.o field.o function.o global.o main-3.15.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.16: iterator.o field.o function.o global.o main-3.16.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.17: iterator.o field.o function.o global.o main-3.17.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.18: iterator.o field.o function.o global.o main-3.18.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

exe-3.19: iterator.o field.o function.o global.o main-3.19.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^


clean:
	rm *.o

run-test: exe-test
	./exe-test

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
	./exe-2.8

run-2.9: exe-2.9
	./exe-2.9

run-2.10: exe-2.10
	./exe-2.10

run-2.11: exe-2.11
	./exe-2.11

run-2.12: exe-2.12
	./exe-2.12

run-2.13: exe-2.13
	./exe-2.13

run-2.14: exe-2.14
	./exe-2.14

run-2.15: exe-2.15
	./exe-2.15

run-3.1: exe-3.1
	./exe-3.1

run-3.2: exe-3.2
	./exe-3.2

run-3.3: exe-3.3
	./exe-3.3

run-3.4: exe-3.4
	./exe-3.4

run-3.5: exe-3.5
	./exe-3.5

run-3.6: exe-3.6
	./exe-3.6

run-3.7: exe-3.7
	./exe-3.7

run-3.8: exe-3.8
	./exe-3.8

run-3.14: exe-3.14
	./exe-3.14

run-3.15: exe-3.15
	./exe-3.15

run-3.16: exe-3.16
	./exe-3.16

run-3.17: exe-3.17
	./exe-3.17

run-3.18: exe-3.18
	./exe-3.18

run-3.19: exe-3.19
	./exe-3.19


PHONY: run-test run-1.1 run-1.3 run-1.4 run-1.5 run-2.1 run-2.2 run-2.3 run-2.4 run-2.5 run-2.6 run-2.7 run-2.8 run-2.9 run-2.10 run-2.11 run-2.12 run-2.13 run-2.14 run-2.15 run-3.1 run-3.2 run-3.3 run-3.4 run-3.5 run-3.6 run-3.7 run-3.8 run-3.14 exe-3.15 exe-3.16 run-3.17 run-3.18 run-3.19 run-progressive clean
