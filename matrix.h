#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
//分块方法，待定
#define block(x, y) (((x) + (y) - 1) / (y))
//矩阵读取
extern int read_matrix(char *filename, double **matrix, int*matrix_size);
extern int read_matrix2(FILE *pfile,  double **matrix, int*matrix_size);
#define print_matrix(matrix, matrix_size) write_matrix(stdout, matrix, matrix_size)
extern int write_matrix(FILE *pfile,  double *matrix, int *matrix_size);
extern int write_matrix2(FILE* pfile,  double *matrix, int matrix_size_x, int matrix_size_y);
//进程的笛卡尔坐标和矩阵块的映射关系
extern int matrix_group(int** block, int column, int row);
extern int matrix_group2(int** block, int column, int row);
//矩阵分块
extern double* matrix_split(double **a, int* dim, const int* block, const int block_num);
extern double* matrix_merge(double **matrix, int*matrix_size, const int* block, int block_num);
//矩阵分块乘法
extern double* matrix1(double* a, double* b, double* out, const int* length);
//矩阵分块乘法包装
#define matrix_multiplus(a, b, out) do{\
	int size[4] = {(int)a[0], (int)a[1], (int)b[0], (int)b[1]};\
	matrix1(a + 2, b + 2, out + 2, size);\
}while(0)
