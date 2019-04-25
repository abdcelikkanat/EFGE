# compiler
CC=g++
CFLAGS=-c -Wall -std=c++11

all: hello

hello: main.o Graph.o Model.o Unigram.o
	$(CC) main.o Graph.o Model.o Unigram.o -o hello

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

Graph.o: Graph.cpp
	$(CC) $(CFLAGS) Graph.cpp

Model.o: Model.cpp
	$(CC) $(CFLAGS) Model.cpp

Unigram.o: Unigram.cpp
	$(CC) $(CFLAGS) Unigram.cpp

clean:
	rm *.o hello

