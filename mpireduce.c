#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define NXPROB 80                /* x dimension of problem grid */
//#define NYPROB 64                 /* y dimension of problem grid */
#define STEPS 2000             /* number of time steps */
#define MIN_CHANGE 0.15

#define NONE MPI_PROC_NULL      /* indicates no neighbor */

#define LTAG 2                  /* message tag */
#define RTAG 3                  /* message tag */
#define UTAG 4                  /* message tag */
#define DTAG 5                  /* message tag */

struct Parms {
	float cx;
	float cy;
} parms = {0.1, 0.1};

int main (int argc, char *argv[]) {
	int	taskid,numtasks,i,j,ix,iy,it,iz;

	if (argc != 7 ) {
		printf("example: ./mpi -t 2 -x 64 -y 80\n");
		exit(1);
	}

	int num_threads, ny_prob, nx_prob;
	for (i=0; i<argc; i++) {
		if (strcmp(argv[i],"-t") == 0){
			num_threads = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i],"-x") == 0){
			nx_prob = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i],"-y") == 0){
			ny_prob = atoi(argv[i + 1]);
		}
	}

	/* First, find out my taskid and how many tasks are running */
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	int numworkers = numtasks-1;

	//create cartecian topology
	int dims[2]={0,0};
	MPI_Dims_create(numtasks,2,dims);
	MPI_Comm mpi_com_cart;
	int periods[2]={0,0};
	MPI_Cart_create (MPI_COMM_WORLD, 2, dims, periods, 1, &mpi_com_cart);

	// find possision of worker
	int cords[2]={0,0};
  	MPI_Comm_rank(mpi_com_cart,&taskid);
	MPI_Cart_coords(mpi_com_cart, taskid, 2, cords);

	//determin rows and cols per worker
	int rows = nx_prob / dims[0];
	int cols = ny_prob / dims[1];

	//determin and distributing extra rows and cols
	int extra_rows = nx_prob % dims[0];
	int extra_cols = ny_prob % dims[1];
	int offset_row = taskid / dims[0];
	int offset_col = taskid % dims[1];
	rows = (offset_row < extra_rows) ? rows+1 : rows;
	cols = (offset_col < extra_cols) ? cols+1 : cols;

	//determin neighbors
	int up,down,left,right;
	cords[0] -= 1;
	if (cords[0] < 0)
		up = NONE;
	else
		MPI_Cart_rank(mpi_com_cart, cords, &up);
	cords[0] += 2;
	if (cords[0] >= dims[0])
		down = NONE;
	else
		MPI_Cart_rank(mpi_com_cart, cords, &down);
	cords[0] -= 1;
	cords[1] -= 1;
	if (cords[1] < 0)
		left = NONE;
	else
		MPI_Cart_rank(mpi_com_cart, cords, &left);
	cords[1] += 2;
	if (cords[1] >= dims[1])
		right = NONE;
	else
		MPI_Cart_rank(mpi_com_cart, cords, &right);
	cords[1] -= 1;

	//add extra rows and cols for halo points
	rows += 2;
	cols += 2;
	float  u[2][rows][cols];

	//parallel reading
	MPI_File fh;
	MPI_Status status;
	int file_offset;
	char string[100];
	int tmp_len;
	int ret_val;

	for (iz=0; iz<2; iz++)
		for (ix=0; ix<rows; ix++)
			for (iy=0; iy<cols; iy++)
				u[iz][ix][iy] = 0.0;

	double read_time = MPI_Wtime();
	double final_read_time;

	MPI_File_open(mpi_com_cart,"outt",MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	file_offset = dims[1]*cords[0]*(rows-2)*(cols-2)*sizeof(float) + cords[1]*(cols-2)*sizeof(float);
	for (i=1; i <= rows-2 ; i++){
		MPI_File_seek(fh,file_offset,MPI_SEEK_SET);
		MPI_File_read(fh,&u[0][i][1],(cols-2),MPI_FLOAT,&status);
		file_offset += dims[1]*(cols-2)*sizeof(float);
	}
	MPI_File_close(&fh);

	read_time = MPI_Wtime() - read_time;
	MPI_Reduce(&read_time, &final_read_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (taskid == 0) {
		printf("parallel read time: %lf\n", final_read_time);		
	}

	//creating cols datatype
	MPI_Datatype my_col_type;
	MPI_Type_vector(rows-2,1,cols,MPI_FLOAT,&my_col_type);
	MPI_Type_commit( &my_col_type );

	//creating rows datatype
	MPI_Datatype my_row_type;
	MPI_Type_contiguous(cols-2, MPI_FLOAT, &my_row_type);
	MPI_Type_commit( &my_row_type );

	//initialize send and rcv messages
	MPI_Request	up_sreq[2],up_rreq[2],down_sreq[2],down_rreq[2],left_sreq[2],left_rreq[2],right_sreq[2],right_rreq[2];
	for(i=0;i<2;i++){
		MPI_Send_init(&u[i][1][cols-2], 1, my_col_type, right, RTAG, mpi_com_cart, &right_sreq[i]);//send to right
		MPI_Send_init(&u[i][1][1], 1, my_col_type, left, LTAG, mpi_com_cart, &left_sreq[i]);//send to left
		MPI_Send_init(&u[i][1][1], 1, my_row_type, up, UTAG, mpi_com_cart, &up_sreq[i]);//send to up
		MPI_Send_init(&u[i][rows-2][1], 1, my_row_type, down, DTAG, mpi_com_cart, &down_sreq[i]);//send to down

		MPI_Recv_init(&u[i][1][cols-1], 1, my_col_type, right, LTAG, mpi_com_cart, &right_rreq[i]);//rcv from right
		MPI_Recv_init(&u[i][1][0], 1, my_col_type, left, RTAG, mpi_com_cart, &left_rreq[i]);//rcv from left
		MPI_Recv_init(&u[i][0][1], 1, my_row_type, up, DTAG, mpi_com_cart, &up_rreq[i]);//rcv from up
		MPI_Recv_init(&u[i][rows-1][1], 1, my_row_type, down, UTAG, mpi_com_cart, &down_rreq[i]);//rcv from down
	}

	//MPI_Barrier(mpi_com_cart);

	//calculating "useful" points
	int start_x,end_x,start_y,end_y;
	int timeNotPrinted = 1;
	if (up != NONE)
		start_x = 1;
	else
		start_x = 2;

	if (down != NONE)
		end_x = rows-2;
	else
		end_x = rows-3;

	if (left != NONE)
		start_y = 1;
	else
		start_y = 2;

	if (right != NONE)
		end_y = cols-2;
	else
		end_y = cols-3;

	int start_x_plus=start_x+1,start_x_minus=start_x-1;
	int end_x_plus=end_x+1,end_x_minus=end_x-1;
	int start_y_plus=start_y+1,start_y_minus=start_y-1;
	int end_y_plus=end_y+1,end_y_minus=end_y-1;
	iz=0;
	double local_time, final_time;
	int diff_flag,flag_sum;
	int reduce_step;

	//starting heat transfer
	local_time=MPI_Wtime();

	for (it = 1; it <= STEPS; it++) {

		//sigglisi flags
		diff_flag = 0;
		reduce_step = (it%(STEPS/50)==0);

		//send outter points
		MPI_Start(&right_sreq[iz]);//send to right
		MPI_Start(&left_sreq[iz]);//send to left
		MPI_Start(&up_sreq[iz]);//send to up
		MPI_Start(&down_sreq[iz]);//send to down

		//receive outter points
		MPI_Start(&right_rreq[iz]);//receive from right
		MPI_Start(&left_rreq[iz]);//receive from left
		MPI_Start(&up_rreq[iz]);//receive from up
		MPI_Start(&down_rreq[iz]);//receive from down

		//compute inner points
		for (ix = 2; ix < rows-2; ix++) {
			for (iy = 2; iy < cols-2; iy++) {
				u[1-iz][ix][iy] = u[iz][ix][iy]  +
					parms.cx * (u[iz][ix+1][iy] +
					u[iz][ix-1][iy] -
					2.0 * u[iz][ix][iy]) +
					parms.cy * (u[iz][ix][iy+1] +
					u[iz][ix][iy-1] -
					2.0 * u[iz][ix][iy]);
				if(reduce_step)
					if (abs(u[1-iz][ix][iy] - u[iz][ix][iy]) > MIN_CHANGE)
						diff_flag++;
			}
		}
		//wait to recieve halo points
		MPI_Wait(&right_rreq[iz],&status);
		MPI_Wait(&left_rreq[iz],&status);
		MPI_Wait(&up_rreq[iz],&status);
		MPI_Wait(&down_rreq[iz],&status);


		if (up != NONE)
			for (iy = start_y; iy <= end_y; iy++) {
				//up halo points
				u[1-iz][start_x][iy] = u[iz][start_x][iy]  +
					parms.cx * (u[iz][start_x_plus][iy] +
					u[iz][start_x_minus][iy] -
					2.0 * u[iz][start_x][iy]) +
					parms.cy * (u[iz][start_x][iy+1] +
					u[iz][start_x][iy-1] -
					2.0 * u[iz][start_x][iy]);
				if(reduce_step)
					if (abs(u[1-iz][start_x][iy] - u[iz][start_x][iy]) > MIN_CHANGE)
						diff_flag++;
			}

		if (down != NONE)
			for (iy = start_y; iy <= end_y; iy++) {
				//down halo points
				u[1-iz][end_x][iy] = u[iz][end_x][iy]  +
					parms.cx * (u[iz][end_x_plus][iy] +
					u[iz][end_x_minus][iy] -
					2.0 * u[iz][end_x][iy]) +
					parms.cy * (u[iz][end_x][iy+1] +
					u[iz][end_x][iy-1] -
					2.0 * u[iz][end_x][iy]);
				if(reduce_step)
					if (abs(u[1-iz][end_x][iy] - u[iz][end_x][iy]) > MIN_CHANGE)
						diff_flag++;
			}

		if (left != NONE)
			for (ix = start_x; ix <= end_x; ix++) {
				//left halo points
				u[1-iz][ix][start_y] = u[iz][ix][start_y]  +
					parms.cx * (u[iz][ix+1][start_y] +
					u[iz][ix-1][start_y] -
					2.0 * u[iz][ix][start_y]) +
					parms.cy * (u[iz][ix][start_y_plus] +
					u[iz][ix][start_y_minus] -
					2.0 * u[iz][ix][start_y]);
				if(reduce_step)
					if (abs(u[1-iz][ix][start_y] - u[iz][ix][start_y]) > MIN_CHANGE)
						diff_flag++;
			}

		if (right != NONE)
			for (ix = start_x; ix <= end_x; ix++) {
				//right halo points
				u[1-iz][ix][end_y] = u[iz][ix][end_y]  +
					parms.cx * (u[iz][ix+1][end_y] +
					u[iz][ix-1][end_y] -
					2.0 * u[iz][ix][end_y]) +
					parms.cy * (u[iz][ix][end_y_plus] +
					u[iz][ix][end_y_minus] -
					2.0 * u[iz][ix][end_y]);
				if(reduce_step)
					if (abs(u[1-iz][ix][end_y] - u[iz][ix][end_y]) > MIN_CHANGE)
						diff_flag++;
			}

		//wait to send
		MPI_Wait(&right_sreq[iz],&status);
		MPI_Wait(&left_sreq[iz],&status);
		MPI_Wait(&up_sreq[iz],&status);
		MPI_Wait(&down_sreq[iz],&status);

		//checking for sigglisi
		if(reduce_step){
			MPI_Allreduce(&diff_flag,&flag_sum,1,MPI_INT,MPI_SUM,mpi_com_cart);
			if (flag_sum == 0) {
				if(timeNotPrinted)
				{
					printf("Reduced at %d\n",it);
					timeNotPrinted = 0;
					local_time = MPI_Wtime() - local_time;

					MPI_Reduce(&local_time, &final_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
					if(taskid == 0) {
						printf("steps:%d nx_prob:%d ny_prob:%d MIN_CHANGE:%f\n",STEPS,nx_prob,ny_prob,MIN_CHANGE);
						printf("threads:%d numtasks:%d\n",num_threads,numtasks);
						printf("Total time pased %lf\n",final_time);
					}
				}
			}
		}


		iz = 1 - iz;
	}

	if(timeNotPrinted)
	{
		timeNotPrinted = 0;
		local_time = MPI_Wtime() - local_time;

		MPI_Reduce(&local_time, &final_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		if(taskid == 0) {
			printf("steps:%d nx_prob:%d ny_prob:%d MIN_CHANGE:%f\n",STEPS,nx_prob,ny_prob,MIN_CHANGE);
			printf("threads:%d numtasks:%d\n",num_threads,numtasks);
			printf("Total time pased %lf\n",final_time);
		}
	}

	local_time = MPI_Wtime();

	//parallel write
	MPI_File_open(mpi_com_cart,"test.bin", MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	file_offset= dims[1]*cords[0]*(rows-2)*(cols-2)*sizeof(float) + cords[1]*(cols-2)*sizeof(float);
	for (i=1; i <= rows-2 ; i++){
		MPI_File_seek(fh,file_offset,MPI_SEEK_SET);
		MPI_File_write(fh,&u[0][i][1],(cols-2),MPI_FLOAT,&status);
		file_offset += dims[1]*(cols-2)*sizeof(float);
	}

	MPI_File_close(&fh);
	//parallel write

	local_time = MPI_Wtime() - local_time;
	MPI_Reduce(&local_time, &final_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (taskid == 0) {
		printf("parallel write time: %lf\n", final_time);
	}

	MPI_Type_free(&my_col_type);
	MPI_Type_free(&my_row_type);

	MPI_Finalize();
}
