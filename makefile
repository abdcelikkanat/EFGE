# compiler
CC=g++
CFLAGS=-c -Wall -std=c++11

all: run

run: main.o Graph.o Node.o Model.o Unigram.o Vocabulary.o
	$(CC) main.o Graph.o Node.o Model.o Unigram.o Vocabulary.o -o run

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

Graph.o: Graph.cpp
	$(CC) $(CFLAGS) Graph.cpp

Node.o: Node.cpp
	$(CC) $(CFLAGS) Node.cpp

Model.o: Model.cpp
	$(CC) $(CFLAGS) Model.cpp

Unigram.o: Unigram.cpp
	$(CC) $(CFLAGS) Unigram.cpp

Vocabulary.o: Vocabulary.cpp
	$(CC) $(CFLAGS) Vocabulary.cpp


clean:
	rm *.o run

