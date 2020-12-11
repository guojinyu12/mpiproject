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
	if(size == 4 || size == 9 || size == 16){
		int rank;//MPI进程序号
		int block_num = (int)sqrt(size);
		int dims[DIMENSION] = {block_num, block_num};//每个维度的元素数
		int periods[DIMENSION] = {1, 1};//循环
		
		double *a= NULL, *b = NULL, *c = NULL;//矩阵
		int matrix_size[4];// 两个矩阵行列数

		MPI_Comm cartcomm;//笛卡尔通信域
		// 2*2笛卡尔坐标坐标
		MPI_Cart_create(MPI_COMM_WORLD, DIMENSION, dims, periods, 1, &cartcomm);
		// 获得坐标和笛卡尔进程编号
		MPI_Comm_rank(cartcomm, &rank);
		//进程 0 读取并发送矩阵数据
		if(rank == 0){
			// 矩阵数据的获取和内存分配{{{
			if(read_matrix("data/a.txt", &a, &matrix_size[0]) == -1 ||
			read_matrix("data/b.txt", &b, &matrix_size[2]) == -1 || 
			matrix_size[1] != matrix_size[2]){
				fprintf(stderr, "there is not file that is fit my need");
				if(a != NULL){
					free(a);
				}
				if(b != NULL){
					free(b);
				}
				MPI_Finalize();
			
			}
			//矩阵数据的获取完毕}}}
			//矩阵分块{{{
			int* block1 = malloc(size*sizeof(int));
			//矩阵A生成映射关系并分块
			matrix_group(&block1, block_num, block_num);
			matrix_split(&a, &matrix_size[0], block1, block_num);
			//矩阵B生成映射关系并分块
			matrix_group2(&block1, block_num, block_num);
			matrix_split(&b, &matrix_size[2], block1, block_num);
			free(block1);
			block1 = NULL;
			//矩阵分块结束}}}
		}//进程 0 独有动作结束

		//广播矩阵的行标和列标
		MPI_Bcast(matrix_size, 4, MPI_INT, 0, cartcomm);
		//计算矩阵所需最少空间{{{
		int local_matrix_size[4] = {// 矩阵块的最大行列数
			block(matrix_size[0], block_num), 
			block(matrix_size[1], block_num), 
			block(matrix_size[2], block_num), 
			block(matrix_size[3], block_num), 
		};
		int local_length[3] = {
			local_matrix_size[0] * local_matrix_size[1] + 2,
			local_matrix_size[2] * local_matrix_size[3] + 2,
			local_matrix_size[0] * local_matrix_size[3] + 2,
		};
		double *local_a = malloc(local_length[0] * sizeof(double));
		double *local_b = malloc(local_length[1] * sizeof(double));
		double *local_c = calloc(local_length[2], sizeof(double)); //保证初始化为0
		// 本地内存准备结束}}}
		// 分发数据
		MPI_Scatter(a, local_length[0], MPI_DOUBLE, local_a, local_length[0], MPI_DOUBLE, 0, cartcomm);
		MPI_Scatter(b, local_length[1], MPI_DOUBLE, local_b, local_length[1], MPI_DOUBLE, 0, cartcomm);
		//提前释放内存，节约空间
		if(rank == 0){
			if(a != NULL) {
				free(a);
				a = NULL;
			}
			if(b != NULL){
				free(b);
				b = NULL;
			}
		}
		// 计算并在进程间传输{{{
		// 矩阵块C的大小
		local_c[0] = local_a[0];
		local_c[1] = local_b[1];
		matrix_multiplus(local_a, local_b, local_c);
		// 两个维度上的邻居编号
		int neighbors[4];//保存邻居
		MPI_Cart_shift(cartcomm, 0, 1, &neighbors[UP], &neighbors[DOWN]);
		MPI_Cart_shift(cartcomm, 1, 1, &neighbors[LEFT], &neighbors[RIGHT]);
		MPI_Request request[4];// 请求对象
		MPI_Status status[4];// 状态对象
		// 从右和下接收，向左和上发送
		for(int i = block_num - 1; i > 0; i--){
			MPI_Sendrecv_replace(local_a, local_length[0], MPI_DOUBLE, neighbors[RIGHT], 1, neighbors[LEFT], 1, cartcomm, &status[1]);
			MPI_Sendrecv_replace(local_b, local_length[1], MPI_DOUBLE, neighbors[DOWN], 1, neighbors[UP], 1, cartcomm, &status[0]);
			matrix_multiplus(local_a, local_b, local_c);
		}
		// 所有运算部分结束}}}
		// 数据的回收
		if(rank == 0){
			c =  malloc(size * local_length[2] * sizeof(double));
		}
		MPI_Gather(local_c, local_length[2], MPI_DOUBLE, c, local_length[2], MPI_DOUBLE, 0, cartcomm);
		// 转换与结果显示
		if(rank == 0){
			int* block1 = malloc(size * sizeof(int));
			for(int i = 0; i < size; i++){
				block1[i] = i;
			}
			matrix_size[1] = matrix_size[3];
			c = matrix_merge(&c, matrix_size, block1, block_num);
			print_matrix(c, matrix_size);
		}
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
		fprintf(stderr, "there is not enough process\n"); 
	MPI_Finalize();//MPI环境结束
	return 0;
}
