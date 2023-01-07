#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "myProto.h"


#define SIZE_SQ1 10000
#define SIZE_SQ2 5000

int main(int argc, char* argv[])
{
	int size, rank, i, maxMin = 0;//maxMin is a filter that will be 1 if we saw maximum in input, and 0 for minimum.
	char sq1[SIZE_SQ1], sq2[SIZE_SQ2], maxminfile[10];
	double w1 = 0, w2 = 0, w3 = 0, w4 = 0;
	MPI_Status status;

	double t1, t2;//to messure the time of the program.

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (size != 2)//the number of the processes is just 2.
	{
		printf("Run the example with only two processes\n");
		MPI_Abort(MPI_COMM_WORLD, __LINE__);
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	// Divide the tasks between th two processes. 
	if (rank == 0)
	{
		t1 = MPI_Wtime();//start time.

		FILE* f = fopen("input.txt", "r");
		if (f)
		{
			//read input from the file, proccess 0.
			fscanf(f, "%lf %lf %lf %lf", &w1, &w2, &w3, &w4);
			fscanf(f, "%s", sq1);
			fscanf(f, "%s", sq2);
			fscanf(f, "%s", maxminfile);
			if (strcmp(maxminfile, "maximum") == 0)//check if minimum or maximum.
				maxMin = 1;//maximum.
			else
				maxMin = 0;//minimum.
		}
		else
		{
			printf("Something wrong with file!");
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		}

		int m = strlen(sq1);//size of seq1 
		int n = strlen(sq2);//size of seq2
		int offset = m - n;//how many times we want to slide the seq2 under sq1.
		int offsetStart = 0;//from where process 0 starts.
		int secondOffsetStart = 0;//from where process 1 starts.
		double optScore = 0;//the optimal score for process 0.
		int optOffset = 0; //the optimal offset for process 0.
		int p2OptOffset = 0; //the optimal offset for process 1.
		char* optMutant; //the optimal mutant to process 0;

		optMutant = (char*)malloc(n * sizeof(char));
		if(optMutant == NULL){
			printf("failed in allocation");
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		}

		if (offset % 2 == 1)//devide the offset between the processes.
			secondOffsetStart = offset / 2 + 1;//the process 1 take from offset/2 to the offset.
		else
			secondOffsetStart = offset / 2;

		MPI_Send(&offset, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Send(&secondOffsetStart, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Send(&w1, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		MPI_Send(&w2, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		MPI_Send(&w3, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		MPI_Send(&w4, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		MPI_Send(sq1, SIZE_SQ1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
		MPI_Send(sq2, SIZE_SQ2, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
		MPI_Send(&maxMin, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);


		offset /= 2;//the process 0 take from 0 to the offset/2 .
		int checkScore = 0;//to compare the pre optscore with the new score.
		char* checkOptMutant;//to compare the pre optmutant with the new score.
		for (i = 0; i <= offset; i++) {
			checkOptMutant = mutant(sq1, sq2, i, w1, w2, w3, w4, maxMin);//get the opt mutant for this offset.
			checkScore = aligmentScore(i, sq1, checkOptMutant, w1, w2, w3, w4);//get the optimal alignScore for this offset.
			if (maxMin == 1) {//maximum code.
				if (i > 0) {//do not compare if we were in the first iteration. 
					if (checkScore > optScore) {//compare the score and mutant and offset with the optimal ones from before.
						optMutant = checkOptMutant;
						optScore = checkScore;
						optOffset = i;
					}
				}
				else {
					optMutant = checkOptMutant;//put the optimal score and mutant and offset.
					optScore = checkScore;
					optOffset = i;
				}
			}
			else {
				if (i > 0) {//minimum code.
					if (checkScore < optScore) {//compare the score and mutant and offset with the optimal ones from before.
						optMutant = checkOptMutant;
						optScore = checkScore;
						optOffset = i;
					}
				}
				else {
					optMutant = checkOptMutant;//put the optimal score and mutant and offset.
					optScore = checkScore;
					optOffset = i;
				}
			}
		}
		MPI_Recv(checkOptMutant, n, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);//recive the results from process 1.
		MPI_Recv(&checkScore, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&p2OptOffset, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
		if (maxMin == 1) {//compare process 1 and process 0 optimal mutant score and offset, and save the optimal mutant in the file.
			if (checkScore > optScore) {
				optMutant = checkOptMutant;
				optScore = checkScore;
				optOffset = p2OptOffset;
			}
		}
		else {
			if (checkScore < optScore) {
				optMutant = checkOptMutant;
				optScore = checkScore;
				optOffset = p2OptOffset;
			}
		}

		FILE* fout = fopen("output.txt", "w");//write the results to file.

		if (fout)
			fprintf(fout, "best mutant: %s \nAlignment Score: %lf  and offset : %d ", optMutant, optScore, optOffset);
		else
			printf("Can't open the output file!!");

		t2 = MPI_Wtime();

		printf("The calculation took %lf seconds\n", t2-t1);

	}

	//The process 1.
	else
	{
		//The second process recive the data from the first process and calculate the sum.
		//the two processes calculate the same way.
		//and after all the second process send to the first process its optimal results.
		int offset = 0;
		int offsetStart = 0;
		double optScore = 0;
		int optOffset = 0;
		char* optMutant;
		
		MPI_Recv(&offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&offsetStart, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&w1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&w2, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&w3, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&w4, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(sq1, SIZE_SQ1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(sq2, SIZE_SQ2, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&maxMin, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);


		int n = strlen(sq2);
		optMutant = (char*)malloc(n * sizeof(char));


		int checkScore = 0;
		char* checkOptMutant;
		for (i = offsetStart; i <= offset; i++) {
			checkOptMutant = mutant(sq1, sq2, i, w1, w2, w3, w4, maxMin);
			checkScore = aligmentScore(i, sq1, checkOptMutant, w1, w2, w3, w4);
			if (maxMin == 1) {
				if (i > offsetStart) {
					if (checkScore > optScore) {
						optMutant = checkOptMutant;
						optScore = checkScore;
						optOffset = i;
					}
				}
				else {
					optMutant = checkOptMutant;
					optScore = checkScore;
					optOffset = i;
				}
			}
			else {
				if (i > offsetStart) {
					if (checkScore < optScore) {
						optMutant = checkOptMutant;
						optScore = checkScore;
						optOffset = i;
					}
				}
				else {
					optMutant = checkOptMutant;
					optScore = checkScore;
					optOffset = i;
				}
			}
		}
		MPI_Send(optMutant, n, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&optScore, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&optOffset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	MPI_Finalize();

	return 0;
}
