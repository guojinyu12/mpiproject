#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
//分块方法，待定
#define block(x, y) (((x) + (y) - 1) / (y))
//矩阵读取
extern int read_matrix(char *restrict filename, double *restrict *restrict matrix,
	   	int*restrict matrix_size);
extern int print_matrix(FILE* pfile,  double *restrict matrix, int*restrict matrix_size);
//进程的笛卡尔坐标和矩阵块的映射关系
extern int matrix_group(int** block, int column, int row);
extern int matrix_group2(int** block, int column, int row);
//矩阵分块
extern double* matrix_split(double **a, int* dim, const int* block, const int block_num);
//矩阵乘法，b转置
extern double* matrix1(double* restrict a, double* restrict b, double* restrict out, const int* length);
//矩阵乘法，b不转置
extern int** matrix2(int** restrict a, int** restrict b, int** restrict out, int axlength, int aylength);
