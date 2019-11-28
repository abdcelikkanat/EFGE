src=./src
build=./build


# compiler
CC=g++
CFLAGS=-c -Wall -std=c++11 -Ilib

.PHONY: all
all: efge

efge: main.o Graph.o Node.o Model.o Unigram.o Vocabulary.o Utilities.o
	$(CC) ${build}/main.o ${build}/Graph.o ${build}/Node.o ${build}/Model.o ${build}/Unigram.o ${build}/Vocabulary.o ${build}/Utilities.o -o efge

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp -o ${build}/main.o

Graph.o: ${src}/Graph.cpp
	$(CC) $(CFLAGS) ${src}/Graph.cpp -o ${build}/Graph.o

Node.o: ${src}/Node.cpp
	$(CC) $(CFLAGS) ${src}/Node.cpp -o ${build}/Node.o

Model.o: ${src}/Model.cpp
	$(CC) $(CFLAGS) ${src}/Model.cpp -o ${build}/Model.o

Unigram.o: ${src}/Unigram.cpp
	$(CC) $(CFLAGS) ${src}/Unigram.cpp -o ${build}/Unigram.o

Vocabulary.o: ${src}/Vocabulary.cpp
	$(CC) $(CFLAGS) ${src}/Vocabulary.cpp -o ${build}/Vocabulary.o

Utilities.o: ${src}/Utilities.cpp
	$(CC) $(CFLAGS) ${src}/Utilities.cpp -o ${build}/Utilities.o

.PHONY: clean
clean:
	rm -r ./build/*.o efge

$(shell   mkdir -p $(build))