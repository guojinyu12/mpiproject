#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define div_ceil(x, y) (((x) + (y) - 1) / (y))
//矩阵 a 的矩阵块和进程的映射关系，注意：实际上， column 必须等于 row, out 总长为row * column
int matrix_group(int** block, int column, int row){
	int* out = *block;
	int k = 0;
	for(int i = 0; i < row; i++){
		k = i * column;
		for(int j = 0; j < column; j++){
			out[k + j] = k + j + i;
			if(i + j >= column){
				out[k + j] -= column;
			}
		}
	}
	return 0;
}
//矩阵 b 的矩阵块和进程的映射关系，注意：实际上， column 必须等于 row, out 总长为row * column
int matrix_group2(int** block, int column, int row){
	int* out = *block;
	int k = 0;
	for(int i = 0; i < row; i++){
		k = i * column;
		for(int j = 0; j < column; j++){
			out[k + j] = k + j + j * column;
			if( i + j >= row){
				out[k + j] -= row * column;
			}
		}
	}
	return 0;
}

/*result 是矩阵， matrix_size 保存 result 的行列个数，block表示分块的顺序，block_num分块的数量的平方根*/
double* matrix_split(double **result, int *matrix_size, const int* block, const int block_num){
	int row = div_ceil(matrix_size[0], block_num); //每一小块的最大行数，可能小于
	int column  = div_ceil(matrix_size[1], block_num);//没一小块的最大列数，可能小于
	//每块两行前两行表示矩阵行列数目，分为num块	
	const int length = row * column + 2;
	int num = block_num * block_num;//分块块数
	double* matrixes = calloc(length * num, sizeof(double));
	double* matrix = *result;

	for(int i = num - 1; i >= 0; i--){
		int tem = i * length + 2;
		int coords[2];//矩阵坐标
		coords[0] = block[i] / block_num * row;
		coords[1] = block[i] % block_num * column;
		//实际矩阵块的行列数
		int row_real = row, column_real = column;
		if(coords[0] + row > matrix_size[0]) 
			row_real = matrix_size[0] - coords[0];
		if(coords[1] + column > matrix_size[1])
			column_real = matrix_size[1] - coords[1];

		for(int j = row_real - 1; j >= 0; --j){
			// 唯一的区别
			memcpy(matrixes + tem + j * column,
					matrix + (coords[0] + j) * matrix_size[1] + coords[1], column_real * sizeof(double));
		}
		matrixes[tem - 2] = row_real;
		matrixes[tem - 1] = column_real;
	}
	
	free(*result);
	*result = matrixes;
	return *result;	
}

//前面函数的逆过程
double* matrix_merge(double **result, int *matrix_size, const int* block, const int block_num){
	int row = div_ceil(matrix_size[0], block_num); //每一小块的最大行数，可能小于
	int column  = div_ceil(matrix_size[1], block_num);//没一小块的最大列数，可能小于
	//每块两行前两行表示矩阵行列数目，分为num块	
	const int length = row * column + 2;
	int num = block_num * block_num;//分块块数
	// 区别
	double* matrixes = *result;
   	double* matrix = calloc(matrix_size[0] * matrix_size[1], sizeof(double));

	for(int i = num - 1; i >= 0; i--){
		int tem = i * length + 2;
		int coords[2];//矩阵坐标
		coords[0] = block[i] / block_num * row;
		coords[1] = block[i] % block_num * column;
		//实际矩阵块的行列数
		int row_real = row, column_real = column;
		if(coords[0] + row > matrix_size[0]) 
			row_real = matrix_size[0] - coords[0];
		if(coords[1] + column > matrix_size[1])
			column_real = matrix_size[1] - coords[1];
		for(int j = row_real - 1; j >= 0; --j){
			// 的区别
			memcpy(matrix + (coords[0] + j) * matrix_size[1] + coords[1],
				   	matrixes + tem + j * column, column_real * sizeof(double));
		}
	}
	// 区别
	free(*result);
	*result = matrix;
	return *result;	
}

//矩阵数据的文件读取
int read_matrix(char *restrict filename, double *restrict *restrict matrix,
	   	int*restrict matrix_size){

	char buffer[512];//文件读取缓冲
	int i = 0;//索引
	FILE *pfile = fopen(filename, "r");
	//读取矩阵的行标和列标
	fgets(buffer, 512, pfile);
	matrix_size[0] = atoi(buffer);
	while(buffer[i++] != ',');
	matrix_size[1] = atoi(&buffer[i]);
	if(matrix_size[0] <= 0 || matrix_size[1] <= 0){
		return -1;//失败
	}
	//矩阵 a 的创建
	int length = matrix_size[0] * matrix_size[1];
	double* ptr = (double*) calloc(length, sizeof(double));//分配内存
	int index = 0;//矩阵存储对象 a 的索引
	//读取、转换并保存数据
	while(fgets(buffer, 512, pfile)){
		char *pstr = buffer;
		char *str_end = buffer;
		while(str_end[0] != '\n'){
			if(index >= length)
				return -1;
			ptr[index++] = strtod(pstr, &str_end);
			pstr = str_end + 1;
		}
		if(index % matrix_size[1] != 0){
			return -1;
		}
	}
	//长度正常
	if(index != length){
		free(ptr);
		return -1;//失败
	}
	*matrix = ptr;
	fclose(pfile);
	return 1;
}

//矩阵数据的文件读取或输入，不做文件结束检查，备用
int read_matrix2(FILE* pfile,  double *restrict *restrict matrix,
	   	int*restrict matrix_size){
	fscanf(pfile, "%d,%d", &matrix_size[0], &matrix_size[1]);
	if(matrix_size[0] <= 0 || matrix_size[1] <= 0){
		return -1;//失败
	}
	*matrix = malloc(matrix_size[0] * matrix_size[1] * sizeof(double));
	for(int i = 0; i < matrix_size[0]; i ++){
		fscanf(pfile, "%lf", &(*matrix)[i * matrix_size[1]]);
		for(int j = 1; j < matrix_size[1]; j++){
			fscanf(pfile, ",%lf", &(*matrix)[i * matrix_size[1] + j]);
		}
	}
	return 1;
}

//矩阵数据的文件写入或显示
int write_matrix2(FILE* pfile,  double *restrict matrix, int matrix_size_x, int matrix_size_y){
	if(pfile == NULL || matrix == NULL)
		return -1;
	// 写入矩阵的行标和列标
	fprintf(pfile, "%d,%d\n", matrix_size_x, matrix_size_y);
	for(int i = 0; i < matrix_size_x; i ++){
		fprintf(pfile, "%.2f", matrix[i * matrix_size_y]);
		for(int j = 1; j < matrix_size_y; j++){
			fprintf(pfile, ",%.2f", matrix[i * matrix_size_y + j]);
		}
		fprintf(pfile, "\n");
	}
	return 1;
}
int write_matrix(FILE* pfile,  double *restrict matrix, int*restrict matrix_size){
	return write_matrix2(pfile, matrix, matrix_size[0], matrix_size[1]);
}

//a * b
double* matrix1(double* restrict a, double* restrict b, double* restrict out,
	   const int* length){
#	pragma omp parallel for
	for(int i = length[0] - 1; i >= 0; i--){
		for(int j = length[3] - 1; j >= 0; j--){
			int x = i * length[3] + j, y = i * length[1];
			for(int k = length[1] - 1; k >= 0; k--){
				out[x] += a[y + k] * b[k * length[3] + j];
			}
		}
	}
	return out;
}
