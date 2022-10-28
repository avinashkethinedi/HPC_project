#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void kronecker_graph(float **mat, int n, long nnz, float *prob)
{
	int i, j, k;
	k = log2(n);
	float c_prob[4];
	c_prob[0] = prob[0];
	for(i=1;i<4;i++)
		c_prob[i] = c_prob[i-1] + prob[i];
	printf("cumulative prob: %f, %f, %f, %f\n", c_prob[0], c_prob[1], c_prob[2], c_prob[3]);
	printf("prob: %f, %f, %f, %f\n", prob[0], prob[1], prob[2], prob[3]);
	/*for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			mat[i][j] = (float)rand()/RAND_MAX;*/
	int nnz_count=0;
	int mat_dim, p_idx, prow, pcol, row, col, collisions=0;
	float p;
	while(nnz_count < nnz)
	{
		i=0;
		mat_dim=n;
		row=0;
		col=0;
		while(i++<=k)
		{
			p = (float)rand()/RAND_MAX;
			p_idx = p<c_prob[0]?0:(p<c_prob[1]?1:(p<c_prob[2]?2:3));
			prow = p_idx/2;
			pcol = p_idx%2;
			mat_dim/=2;
			row = row + prow*mat_dim;
			col = col + pcol*mat_dim;
			//printf("\t%.2f, %d, (%d, %d), (%d, %d), %d\n", p, p_idx, prow, pcol, row, col, mat_dim);
		}
		//printf("%d, %0.2f, (%d, %d)\n", nnz_count, p, row, col);
		if(!mat[row][col])
		{
			mat[row][col] = p;
			nnz_count++;
		}
		else
			collisions++;
		//mat[row][col] = p;
		//nnz_count++;
	}
	printf("collisions: %d\n", collisions);
}
void print_mat(float **mat, int n)
{
	int i, j;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			printf("%0.3f ", mat[i][j]);
		printf("\n");
	}
}
void prob_dist(float **mat, int n)
{
	int count[] = {0, 0, 0, 0};
	for(int x=0;x<2;x++)
		for(int y=0;y<2;y++)
		{
			int dim = n/2;
			int row_s = x*dim;
			int col_s = y*dim;
			int p_idx = x*2+y;
			printf("%d %d\n", row_s, col_s);
			for(int i=row_s;i<row_s+dim;i++)
				for(int j=col_s;j<col_s+dim;j++)
					if(mat[i][j]) count[p_idx]++;
		}
	printf("distribution: %d %d %d %d\n", count[0], count[1], count [2], count[3]);

}
void main(int argc, char **argv)
{
	int i, j, n, edge_factor;
	long nnz;
	float prob[] = {0.45, 0.2, 0.2, 0.15};
	n = 1<<(argc>=2?atoi(argv[1]):3);
	if((n & (n-1)) || n == 1)
	{
		printf("matrix dimensions:(%d, %d) aren't in powers of 2\n", n, n);
		exit(1);
	}
	edge_factor = argc>=3?atoi(argv[2]):7;
	nnz = n*edge_factor;
	printf("%d\n", nnz);
	nnz = nnz>(long)n*n?(long)n*n/2:nnz;
	printf("n: %d nnz: %ld\n", n, nnz);
	printf("size of size_t: %d\n", sizeof(size_t));
	size_t mat_size = n;
	mat_size *= mat_size;
	float **mat = (float**)malloc(sizeof(float*)*n);
    mat[0] = (float*)calloc(mat_size, sizeof(float));
    for(i=1;i<n;i++)
        mat[i] = mat[0]+i*n;
	kronecker_graph(mat, n, nnz, prob);
	//print_mat(mat, n);
	prob_dist(mat, n);
}
