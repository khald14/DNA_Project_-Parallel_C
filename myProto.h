#pragma once
#include <omp.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>



char* mutant(char* sq1, char* sq2, int index, double w1, double w2, double w3, double w4, int maxMin);
int checkInString(char* str, char x, char y);
int checkInConservativeGroup(char x, char y);
int checkInSemiConservativeGroup(char x, char y);
double aligmentScore(int index, char* sq1, char* sq2, double w1, double w2, double w3, double w4);
char changeToPoint(char c1, char c2);
char changeToColon(char c1, char c2);

