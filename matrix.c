#include <stdio.h>
#include "mpi.h"

#define DIMENSION 2
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int** matrix1(int** restrict a, int** restrict b, int** restrict out, int axlength, int aylength);

int main(int argc, char* argv[]){
	int size;//进程数
	int rank;//MPI进程序号
	int rank2;//MPI进程序号
	int dims[DIMENSION] = {2,2};//每个维度的元素数
	int periods[DIMENSION] = {1, 1};//循环
	MPI_Comm cartcomm;//笛卡尔通信域

	scanf("%d", &size);
	MPI_Init(&argc, &argv);//MPI环境
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if(size == 4){
		int neighbors[4];//保存邻居
		int coords[2];//坐标
		//2*2笛卡尔坐标坐标
		MPI_Cart_create(MPI_COMM_WORLD, DIMENSION, dims, periods, 0, &cartcomm);
		//获得坐标和笛卡尔进程编号
		MPI_Comm_rank(cartcomm, &rank);
		MPI_Cart_coords(cartcomm, rank, 2, coords);
		MPI_Cart_rank(cartcomm, coords, &rank2);
		//两个维度上的邻居编号
		MPI_Cart_shift(cartcomm, 0, 1, &neighbors[UP], &neighbors[DOWN]);
		MPI_Cart_shift(cartcomm, 1, 1, &neighbors[LEFT], &neighbors[RIGHT]);

		int inbuffer[4];//输入缓冲
		int outbuffer[4]={rank};//输出缓冲
		MPI_Request request[4];//请求对象
		MPI_Status status[4];//状态对象
		//从右和下接收，向左和上发送
		MPI_Isend(&outbuffer[0], 1, MPI_INT, neighbors[DOWN], 1, cartcomm, &request[DOWN]);
		MPI_Irecv(&inbuffer[UP], 1, MPI_INT, neighbors[UP], 1, cartcomm, &request[UP]);
		MPI_Isend(&outbuffer[0], 1, MPI_INT, neighbors[RIGHT], 1, cartcomm, &request[RIGHT]);
		MPI_Irecv(&inbuffer[LEFT], 1, MPI_INT, neighbors[LEFT], 1, cartcomm, &request[LEFT]);
		MPI_Waitall(4, request, status);//等待接收完成

		//结果
		printf("rank = %d rank2 = %d coords = (%d, %d) neighbors = %d,%d,%d,%d ", rank, rank2,
				coords[0], coords[1], neighbors[0], neighbors[1], neighbors[2], neighbors[3]);
		printf("up = %d, left = %d\n", inbuffer[UP], inbuffer[LEFT]);
	}
	else
		fprintf(stderr, "error"); 
	MPI_Finalize();//MPI环境结束
	return 0;
}

//a * b, b 输入前倒置
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
