#include<math.h>
typedef struct block_probablity
{
	float a, b, c, d;
	long nnz;
}block;
typedef struct compressed_sparse_row
{
	int *row_ptr, *col_ptr, num_rows;
	float *val_ptr;
	long nnz;
}CSR;
long** block_allocation(int n, int m)
{
	size_t size = (size_t)n*m;
	long **pe_blocks = (long**)malloc(sizeof(long*)*n);
	pe_blocks[0] = (long*)malloc(sizeof(long)*size);
	for(int i=1;i<n;i++)
		pe_blocks[i] = pe_blocks[0]+(size_t)i*m;
	return pe_blocks;
}
long nnz_distribution(long**, int, int);
long** calculate_prob_distribution(int mat_blocks, int rank, int npes, block *mat_prop)
{
	if(mat_blocks & (mat_blocks-1))
	{
 		printf("matrix blocks %ld isn't in powers of 2\n", mat_blocks);
 		exit(1);
 	}
	long **pe_blocks = block_allocation(mat_blocks, mat_blocks);
	int i, j, k, l, m, x, y;
	int n=log2(mat_blocks);
	for(i=0;i<mat_blocks;i++)
		for(j=0;j<mat_blocks;j++)
			pe_blocks[i][j] = mat_prop->nnz;
	int num_blocks, grid_size, block_row, block_col, grid_row, grid_col;
	float prob;
	for(i=0;i<n;i++)
	{
		num_blocks = 1<<i;
		grid_size = 1<<(n-1-i);
		for(j=0;j<num_blocks;j++)
		{
			block_row = j*grid_size*2;
			for(k=0;k<num_blocks;k++)
			{
				block_col = k*grid_size*2;
				for(l=0;l<2;l++)
				{
					grid_row = block_row+l*grid_size*l;
					for(m=0;m<2;m++)
					{
						grid_col = block_col+grid_size*m;
						prob = l?(m?mat_prop->d:mat_prop->c):(m?mat_prop->b:mat_prop->a);
						for(x=0;x<grid_size;x++)
						{
							for(y=0;y<grid_size;y++)
								pe_blocks[grid_row+x][grid_col+y] = (long double)round(pe_blocks[grid_row+x][grid_col+y]*prob);
								//pe_blocks[grid_row+x][grid_col+y] *= prob;
						}
					}
				}
			}
		}
	}
	//if(!rank) nnz_distribution(pe_blocks, mat_blocks, mat_blocks);
	n = mat_blocks/npes;
	long **nnz_dist = block_allocation(n, mat_blocks);
	size_t size = mat_blocks*n*sizeof(long);
	memcpy(nnz_dist[0], pe_blocks[n*rank], size);
	/*for(i=0;i<npes;i++)
	{
		if(rank == i)nnz_distribution(nnz_dist, n, mat_blocks);
		MPI_Barrier(MPI_COMM_WORLD);
	}*/
	free(pe_blocks[0]);
	free(pe_blocks);
	return nnz_dist;
}
long nnz_distribution(long **pe_blocks, int n, int m)
{
	int i, j;
	long nnz=0;
	for(i=0;i<n;i++)
	{
		//printf("PE: %d -> ", i);
		for(j=0;j<m;j++)
		{
			//printf("%ld ", pe_blocks[i][j]);
			nnz += pe_blocks[i][j];
		}
		//printf("\n");
	}
	//printf("nnz: (%ld)\n", nnz);
	return nnz;
}

