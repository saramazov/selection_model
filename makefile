main:main.o function.o
	g++ main.o function.o -o main
main.o: main.cpp function.h
	g++ -c main.cpp -o main.o -std=c++0x
function.o: function.cpp function.h
	g++ -c function.cpp -o function.o -std=c++0x
