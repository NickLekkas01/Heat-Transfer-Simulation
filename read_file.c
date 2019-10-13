#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define NXPROB 80                 /* x dimension of problem grid */
//#define NYPROB 64                 /* y dimension of problem grid */

struct Parms { 
	float cx;
	float cy;
} parms = {0.1, 0.1};

int main (int argc, char *argv[]) {
	int i,j;
	float num;
	
	if (argc != 7 ) {
		printf("usage: ./read_file -f <filename> -x <NXPROB> -y <NYPROB>\n");
		return 1;
	}

	char filename[150];
	FILE *ptr;
	int num_threads, ny_prob, nx_prob;
	for (i=0; i<argc; i++) {
		if (strcmp(argv[i],"-f") == 0){
			strcpy(filename,argv[2]);
		}
		if (strcmp(argv[i],"-x") == 0){
			nx_prob = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i],"-y") == 0){
			ny_prob = atoi(argv[i + 1]);
		}
	}
	float u[nx_prob][ny_prob];
	ptr = fopen(filename,"rb");

	fread(&u,sizeof(float),nx_prob*ny_prob,ptr);
	for (i=0; i<nx_prob; i++) {
		for(j=0; j<ny_prob; j++) {
			printf("%6.1f",u[i][j]);
			if (j != ny_prob-1) {
				printf(" ");
			}
		}
		printf("\n");
	}
	return 0;
}