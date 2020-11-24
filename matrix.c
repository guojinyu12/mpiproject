#include <stdio.h>
#include "mpi.h"

//维数
#define DIMENSION 2
#define DIMENSION_SIZE 2
#define DIMENSION_SIZE_POW_2 DIMENSION_SIZE * DIMENSION_SIZE
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int** matrix1(int** restrict a, int** restrict b, int** restrict out, int axlength, int aylength);

int main(int argc, char* argv[]){
	int size;//进程数
	int rank;//MPI进程序号
	int dims[DIMENSION] = {DIMENSION_SIZE, DIMENSION_SIZE};//每个维度的元素数
	int periods[DIMENSION] = {1, 1};//循环
	int dimension_size[2];// 矩阵行列数
	MPI_Comm cartcomm;//笛卡尔通信域
	MPI_Init(&argc, &argv);//MPI环境
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if(size == 4){
		int coords[2] = {0, 0};//坐标
		int neighbors[4];//保存邻居
		// 2*2笛卡尔坐标坐标
		MPI_Cart_create(MPI_COMM_WORLD, DIMENSION, dims, periods, 1, &cartcomm);
		// 获得坐标和笛卡尔进程编号
		MPI_Comm_rank(cartcomm, &rank);
		if(rank <= 1){
			int rank2;//MPI进程序号
			for(int i = 0; i < DIMENSION_SIZE_POW_2; i++){
				/* MPI_Cart_coords(cartcomm, rank, DIMENSION, coords); */
				coords[1]++;
				if(coords[1] == DIMENSION_SIZE){
					coords[1] = 0;
					coords[0]++;
				}
		 		MPI_Cart_rank(cartcomm, coords, &rank2);
				if(rank == 0){
					dimension_size[0] = (40 + DIMENSION_SIZE - 1) / DIMENSION_SIZE;
					dimension_size[1] = (40 + DIMENSION_SIZE - 1) / DIMENSION_SIZE;
					MPI_Bcast(dimension_size, 2, MPI_INT, 0, cartcomm);
				}
				else{// rank 等于1
					// pass
				}
				//获取数据

			}
		}
		if(rank != 0)
			MPI_Bcast(dimension_size, 2, MPI_INT, 0, cartcomm);
		
		
		// 两个维度上的邻居编号
		MPI_Cart_shift(cartcomm, 0, 1, &neighbors[UP], &neighbors[DOWN]);
		MPI_Cart_shift(cartcomm, 1, 1, &neighbors[LEFT], &neighbors[RIGHT]);
		MPI_Request request[4];// 请求对象
		MPI_Status status[4];// 状态对象
		double inbuffer[4];// 输入缓冲
		double outbuffer[4]={ rank, rank };// 输出缓冲
		// 从右和下接收，向左和上发送
		MPI_Barrier(cartcomm);
		// outbuffer[0]++;
		MPI_Isend(&outbuffer[0], 1, MPI_DOUBLE, neighbors[DOWN], 1, cartcomm, &request[DOWN]);
		MPI_Irecv(&inbuffer[UP], 1, MPI_DOUBLE, neighbors[UP], 1, cartcomm, &request[UP]);
		// outbuffer[1]++;
		MPI_Isend(&outbuffer[1], 1, MPI_DOUBLE, neighbors[RIGHT], 1, cartcomm, &request[RIGHT]);
		MPI_Irecv(&inbuffer[LEFT], 1, MPI_DOUBLE, neighbors[LEFT], 1, cartcomm, &request[LEFT]);
		MPI_Waitall(4, request, status);// 等待接收完成
		MPI_Barrier(cartcomm);

		// 结果
		printf("rank = %d coords = (%d, %d) neighbors = %d,%d,%d,%d ", rank,
				coords[0], coords[1], neighbors[0], neighbors[1], neighbors[2], neighbors[3]);
		printf("up = %f, left = %f dimension_size = %d * %d\n", inbuffer[UP], inbuffer[LEFT], 
				dimension_size[0], dimension_size[1]);
	}
	else
		fprintf(stderr, "there is not enough process"); 
	MPI_Finalize();//MPI环境结束
	return 0;
}

//a * b, b 输入前转置
int** matrix1(int** restrict a, int** restrict b, int** restrict out, int axlength, int aylength){
#	pragma omp parallel for
	for(int i = axlength - 1; i >= 0; i--){
		for(int j = axlength - 1; j >= 0; j--){
			for(int k = aylength - 1; k >= 0; k--){
				out[i][j] += a[i][k]*b[j][k];
			}
		}
	}
	return out;
}
//a * b，b不转置
int** matrix2(int** restrict a, int** restrict b, int** restrict out, int axlength, int aylength){
#	pragma omp parallel for
	for(int i = axlength - 1; i >= 0; i--){
		for(int j = axlength - 1; j >= 0; j--){
			for(int k = aylength - 1; k >= 0; k--){
				out[i][j] += a[i][k]*b[k][j];
			}
		}
	}
	return out;
}
