################################
# Makefile
#
# author: Yan-Ru Ju
# edited by: 10/2022
################################

CC=g++
DEPS=main.h
CFLAGS=-std=c++11 -O3
.PHONY : clean hw2
objects=../bin/hw2
main: main.cpp $(DEPS)
		$(CC) main.cpp $(CFLAGS) -o $(objects)
clean:
	-rm $(objects)
