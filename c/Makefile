#!/bin/bash

OBJ = test_kde.o
 
all : $(OBJ)
	gcc -std=c99 -g $(OBJ) -o test_kde -lm -lgsl -lgslcblas
 
%.o: %.c
	gcc -std=c99 -c $^
