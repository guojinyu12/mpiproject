#include "mpi.h"
#include "matrix.h"

//维数
#define DIMENSION 2
//每个维度的元素个数
#define DIMENSION_SIZE 2
//总的元素个数
#define DIMENSION_SIZE_POW_2 (DIMENSION_SIZE * DIMENSION_SIZE)
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int main(int argc, char* argv[]){
	int size;//进程数
	MPI_Init(&argc, &argv);//MPI环境
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if(size == 4){
		int rank;//MPI进程序号
		int dims[DIMENSION] = {DIMENSION_SIZE, DIMENSION_SIZE};//每个维度的元素数
		int periods[DIMENSION] = {1, 1};//循环
		int matrix_size[4];// 两个矩阵行列数
		int coords[2] = {0, 0};//坐标
		int neighbors[4];//保存邻居
		double *local_a = NULL, *local_b = NULL, *local_c = NULL;//矩阵
		double *a= NULL, *b = NULL, *c = NULL;//矩阵

		MPI_Comm cartcomm;//笛卡尔通信域
		// 2*2笛卡尔坐标坐标
		MPI_Cart_create(MPI_COMM_WORLD, DIMENSION, dims, periods, 1, &cartcomm);
		// 获得坐标和笛卡尔进程编号
		MPI_Comm_rank(cartcomm, &rank);
		//进程 0 读取并发送矩阵数据
		if(rank == 0){
			int* block;
			// 矩阵数据的获取和内存分配
			if(read_matrix("data/a.txt", &a, &matrix_size[0]) == -1 ||
			read_matrix("data/b.txt", &b, &matrix_size[2]) == -1 || 
			matrix_size[1] != matrix_size[2]){
				fprintf(stderr, "there is not file that is fit my need");
				if(a != NULL){
					free(local_a);
				}
				if(b != NULL){
					free(b);
				}
				MPI_Finalize();
			
			}
			int block_num = (int)sqrt(size);
			int x = block(matrix_size[0], block_num), y = block(matrix_size[1], block_num);
			int local_size[2];
			matrix_group(&block, block_num, block_num);
			for(int i = 0; i < size; i++){
				printf("%d,", block[i]);
			}
			printf("\n");
			matrix_split(&a, matrix_size, block, block_num);
			free(block);
			block = NULL;
			/* //矩阵显示 */
			/* for(int i = 0; i < size; i++){ */
			/* 	local_size[0] = (int)a[(x * y + 2) * i]; */
			/* 	local_size[1] = (int)a[(x * y + 2) * i + 1]; */
			/* 	print_matrix(stdout, &a[(x * y + 2) * i + 2], local_size); */
			/* 	printf("\n"); */
			/* } */

			// c = (double*)calloc(matrix_size[0] * matrix_size[3], sizeof(double));

			// 矩阵乘法测试
			/* matrix1(a, b, c, matrix_size); */
			/* for(int i = 0; i < matrix_size[0]; i++){ */
			/* 	for(int j = 0; j < matrix_size[3]; j++){ */
			/* 		printf("%.0f ", c[i * matrix_size[3]+ j]); */
			/* 	} */
			/* 	printf("\n"); */
			/* } */
			//矩阵数据的获取完毕

			/* int rank2;//MPI进程序号 */
			/* for(int i = 0; i < DIMENSION_SIZE_POW_2; i++){ */
			/* 	/1* MPI_Cart_coords(cartcomm, rank, DIMENSION, coords); *1/ */
			/* 	coords[1]++; */
			/* 	if(coords[1] == DIMENSION_SIZE){ */
			/* 		coords[1] = 0; */
			/* 		coords[0]++; */
			/* 	} */
		 		/* MPI_Cart_rank(cartcomm, coords, &rank2); */

			/* } */
			//广播矩阵的行标和列标
			MPI_Bcast(matrix_size, 4, MPI_INT, 0, cartcomm);
		}//进程 0 独有动作结束
		else{
			//广播接收
			MPI_Bcast(matrix_size, 4, MPI_INT, 0, cartcomm);

			/* int length[3]; */
			/* length[0] = (matrix_size[0] + DIMENSION_SIZE - 1) / DIMENSION_SIZE */ 
			/* 	* (matrix_size[1] + DIMENSION_SIZE - 1) / DIMENSION_SIZE * sizeof(double); */
			/* local_a = (double*) malloc(length[0]); */
			/* length[1] = (matrix_size[2] + DIMENSION_SIZE - 1) / DIMENSION_SIZE */ 
			/* 	* (matrix_size[3] + DIMENSION_SIZE - 1) / DIMENSION_SIZE; */
			/* local_b = (double*) malloc(length[1]); */
			/* length[2] = (matrix_size[0] + DIMENSION_SIZE - 1) / DIMENSION_SIZE */ 
			/* 	* (matrix_size[3] + DIMENSION_SIZE - 1) / DIMENSION_SIZE; */
			/* local_c = (double*) malloc(length[2]); */
		}
		
		
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
		printf("up = %f, left = %f matrix_size = ((%d * %d), (%d, %d))\n", inbuffer[UP], inbuffer[LEFT], 
				matrix_size[0], matrix_size[1], matrix_size[2], matrix_size[3]);
		/* if(rank == 0){ */
		/* 	for(int i = 0; i < matrix_size[0]; ++i){ */
		/* 		printf("(%d)", i); */
		/* 		for(int j = 0; j < matrix_size[1]; ++j){ */
		/* 			printf("%.0f, ", local_a[i * matrix_size[0] + j]); */
		/* 		} */
		/* 		printf("\n"); */
		/* 	} */
		/* } */
		// {{{ 释放内存
		if(a != NULL) {
			free(a);
			a = NULL;
		}
		if(b != NULL){
		   	free(b);
			b = NULL;
		}
		if(c != NULL){
		   	free(c);
			c = NULL;
		}
		if(local_a != NULL) {
			free(local_a);
			local_a = NULL;
		}
		if(local_b != NULL){
		   	free(local_b);
			local_b = NULL;
		}
		if(local_c != NULL){
		   	free(local_c);
			local_c = NULL;
		}// }}} 释放内存结束
	}
	else
		fprintf(stderr, "there is not enough process"); 
	MPI_Finalize();//MPI环境结束
	return 0;
}
