#include "matrix.h"
//注意：实际上， column 必须等于 row, out 总长为row * column
int matrix_group(int** block, int column, int row){
	int* out = *block;
	int k = 0;
	for(int i = 0; i < row; i++){
		k = i * column;
		for(int j = 0; j < column; j++){
			out[k + j] = k + j + i;
			if(i + j >= column){
				out[k + j] -= i * column;
			}
		}
	}
	return 0;
}
int matrix_group2(int** block, int column, int row){
	int* out = *block;
	int k = 0;
	for(int i = 0; i < row; i++){
		k = i * column;
		for(int j = 0; j < column; j++){
			out[k + j] = k + j + j * column;
			if( i + j >= row){
				out[k + j] -= (row + j - 1) * column;
			}
		}
	}
	return 0;
}
/* int matrix_group_destr */
/*matrix 是矩阵， matrix_size 保存 matrix 的行列个数，block表示分块的顺序，block_num分块的数量的平方根*/
double* matrix_split(double **matrix, int* matrix_size, const int* block, const int block_num){
	int row = block(matrix_size[0], block_num); //每一小块的最大行数，可能小于
	int column  = block(matrix_size[1], block_num);//没一小块的最大列数，可能小于
	int num = block_num * block_num;//分块块数
	//每块两行前两行表示矩阵行列数目，分为num块	
	int length = row * column + 2;
	double* matrixes = calloc(length * num, sizeof(double));
	for(int i = num - 1; i >= 0; i--){
		int tem = block[i] * length + 2;
		int coords[2];//矩阵坐标
		coords[0] = i / block_num * row;
		coords[1] = i % block_num * column;
		//实际矩阵块的行列数
		if(coords[0] + row > matrix_size[0]) 
			row = matrix_size[0] - coords[0];
		if(coords[1] + column > matrix_size[1])
			column = matrix_size[1] - coords[1];
		for(int j = row - 1; j >= 0; --j){
			memcpy(matrixes + tem + j * column, 
					*matrix + (coords[0] + j) * matrix_size[1] + coords[1], column * sizeof(double));
		}
		matrixes[tem - 2] = row;
		matrixes[tem - 1] = column;
	}
	
	free(*matrix);
	*matrix = matrixes;
	return *matrix;	
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

//矩阵数据的文件写入或显示
int print_matrix(FILE* pfile,  double *restrict matrix, int*restrict matrix_size){
	if(pfile == NULL || matrix == NULL)
		return -1;
	// 写入矩阵的行标和列标
	fprintf(pfile, "%d,%d\n", matrix_size[0], matrix_size[1]);
	for(int i = 0; i < matrix_size[0]; i ++){
		fprintf(pfile, "%.2f", matrix[i * matrix_size[1]]);
		for(int j = 1; j < matrix_size[1]; j++){
			fprintf(pfile, ",%.2f", matrix[i * matrix_size[1] + j]);
		}
		fprintf(pfile, "\n");
	}
	return 0;
}

//a * b, b 输入
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
