#include <stdlib.h>
#include "matrix.h"

int main(void){
	double *a= NULL, *b = NULL, *c = NULL;//矩阵
	int matrix_size[4];// 两个矩阵行列数
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
	}
	c = calloc(matrix_size[0] * matrix_size[3], sizeof(double));
	matrix1(a, b, c, matrix_size);
	matrix_size[1] = matrix_size[3];
	print_matrix(c, matrix_size);
	free(a); a = NULL;
	free(b); b = NULL;
	free(c); c = NULL;
	return 0;
}

