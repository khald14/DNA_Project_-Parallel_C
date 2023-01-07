#include "myProto.h"
#include <stdio.h>

char* mutant(char* sq1, char* sq2, int index, double w1, double w2, double w3, double w4, int maxMin)
{
	w2 *= -1; // multiply with -1 to mach alignment_score.
	w3 *= -1;
	w4 *= -1;

	int n = strlen(sq2);
	char* optSq2;
	int i = 0;
	if ((optSq2 = (char*)malloc(n*sizeof(char)))== NULL){
		printf("error allocationg opt");
	}

	if (maxMin == 1) {//if the input is maximum.
#pragma omp parallel for num_threads(4)
		for (i = 0; i < n ; i++) {//run this for with 4 threads.
			if (sq1[i+index] != sq2[i] && !(checkInConservativeGroup(sq1[i+index], sq2[i])))
				optSq2[i] = sq1[i+index];// we can change every point and space to star, according to the rules.
			else {
				if (sq1[i+index] == sq2[i])//if star then we do not make any changes.
					optSq2[i] = sq2[i];
				else {
					if (w2 < w3 || w2 < w4) {//if colon, then we can't change it to star, according to the rules. 
						if (w3 < w4)//so we check if it is more optimal to change to point or space.
							optSq2[i] = 'B'; // B is not in any group in Conservative Groups, change colon to point 
						if (w4 <= w3)
							optSq2[i] = changeToPoint(sq1[i+index], sq2[i]);//change colon to point.
					}
					else
						optSq2[i] = sq2[i];//do not make changes if the colon is more opt to the maximum than point and space.
				}
			}
		}
		optSq2[n] = '\0';//end the mutant with \0 to make it string and return it witout errors.
		return optSq2;
	}
	else {//if the input is minimum.

#pragma omp parallel for num_threads(4)
		for (i = index; i < n + index; i++) {
			if ((w2 <= w3 && w3 <= w4) || (w2 <= w4 && w4 <= w3)) { //change to colon
				if (!checkInConservativeGroup(sq1[i], sq2[i - index]))
					optSq2[i - index] = changeToColon(sq1[i], sq2[i - index]);
				else
					optSq2[i - index] = sq2[i - index];
			}
			else {
				if ((w3 <= w2 && w2 <= w4) || (w3 <= w4 && w4 <= w2)) {//change to point 
					if (!checkInSemiConservativeGroup(sq1[i], sq2[i - index]))
						optSq2[i - index] = changeToPoint(sq1[i], sq2[i - index]);
					else
						optSq2[i - index] = sq2[i - index];
				}
				else {
					if (sq1[i] != 'B')
						optSq2[i - index] = 'B';
					else
						optSq2[i - index] = 'Z';
				}
			}
		}
	}
	optSq2[n] = '\0';
	return optSq2;
}


int checkInString(char* str, char x, char y)//check if to chars in string.
{
	if (strchr(str, x))
		if (strchr(str, y))
			return 1;

	return 0;
}

int checkInConservativeGroup(char x, char y)//check if two chars in a word in the Conservative Group
{
	char conservativeGroup[9][5] = { "NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF" };
	int i = 0, n = 9;
	for (i = 0; i < n; i++)
	{
		if (x != y && checkInString(conservativeGroup[i], x, y))
			return 1;
	}
	return 0;
}

int checkInSemiConservativeGroup(char x, char y)//check if two chars in a word in the Semi-Conservative Group
{
	char SemiConservativeGroup[11][7] = { "SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM" };

	int i = 0, n = 11;

	for (i = 0; i < n; i++)
	{
		if (x != y && checkInString(SemiConservativeGroup[i], x, y))
			return 1;
	}
	return 0;
}

double aligmentScore(int index, char* sq1, char* sq2, double w1, double w2, double w3, double w4)//calculate the alignment score.
{
	double s = 0;
	int m = strlen(sq2), i = 0;
	int sumOfstars = 0;
	int sumOfCons = 0;
	int sumOfSemiCons = 0;
	int sumOfSpace = 0;
	int countCountiner[16] = { 0 };//we make an array with length (num_of_threads* 4 ( there are 4 ws)).

#pragma omp parallel for num_threads(4)
	for (i = index; i < m + index; i++)// for that runs on sq2 and sq1 according to the offset(m = size of sq2).
	{//every thread writes on his side in the array to prevent that two threads write on the same position in the array.
		if (sq1[i] == sq2[i - index])
			countCountiner[4 * omp_get_thread_num() + 0]++;
		else if (checkInConservativeGroup(sq1[i], sq2[i - index]))
			countCountiner[4 * omp_get_thread_num() + 1]++;
		else if (checkInSemiConservativeGroup(sq1[i], sq2[i - index]))
			countCountiner[4 * omp_get_thread_num() + 2]++;
		else
			countCountiner[4 * omp_get_thread_num() + 3]++;
	}

	//we take the total counts from each process and sum it to new variable that presents a sign.
	for (i = 0; i < 16; i++)
	{
		if (i % 4 == 0)
			sumOfstars += countCountiner[i];
		else if (i % 4 == 1)
			sumOfCons += countCountiner[i];
		else if (i % 4 == 2)
			sumOfSemiCons += countCountiner[i];
		else
			sumOfSpace += countCountiner[i];
	}
	//Calculate the alignment score. 
	s = w1 * (double)sumOfstars - w2 * (double)sumOfCons - w3 * (double)sumOfSemiCons - w4 * (double)sumOfSpace;

	return s;
}

char changeToPoint(char c1, char c2) {// in this func we check if there any chars that can made a point withvthe specific char  
	char lettersInSCG[19] = { 'N','D','E','Q','K','S','T','A','M','I','L','V','H','R','F','Y','P','G','C' };
	int i = 0;// we do not check all the alphabitic chars, just the chars that appear in Semi Conservative Groups,
	char c;
	for (i = 0; i < 19; i++) {
		c = lettersInSCG[i];
		if (checkInSemiConservativeGroup(c1, c) && !checkInConservativeGroup(c1, c) && !checkInConservativeGroup(c2, c))
			return c;
	}
	return c2;
}

char changeToColon(char c1, char c2) {// in this func we check if there any chars that can made a colon withvthe specific char  
	char lettersInCG[17] = { 'N','D','E','Q','K','S','T','A','M','I','L','V','H','R','F','Y','W' };
	int i;// we do not check all the alphabitic chars, just the chars that appear in Conservative Groups.
	char c;
	for (i = 0; i < 17; i++) {
		c = lettersInCG[i];
		if (checkInConservativeGroup(c1, c) && !checkInConservativeGroup(c2, c))
			return c;
	}
	return c2;
}