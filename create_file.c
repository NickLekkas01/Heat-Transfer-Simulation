#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//!! anapoda
//#define NXPROB 80                 /* y dimension of problem grid */
//#define NYPROB 64                /* x dimension of problem grid */

void inidat(int nx, int ny, float *u) {
	int ix, iy;
	for (ix = 0; ix <= nx-1; ix++) {
		for (iy = 0; iy <= ny-1; iy++) {
			 *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
		 }
	}
}

int main (int argc, char *argv[]) {
	int ix, iy, i;
	FILE *a;

	if (argc != 7 ) {
		printf("usage: ./create_file -f <filename> -x <NXPROB> -y <NYPROB>\n");
		return 1;
	}

	char filename[150];
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

	//initialize with "random" values
	inidat(nx_prob,ny_prob,(float*)u);

	a = fopen(filename, "w" );
	fwrite(u,nx_prob*ny_prob,sizeof(float),a);
	printf("created outt file with %d %d\n", nx_prob, ny_prob);
	return 0;
}
