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
long calculate_nnz(long*, int);
long* calculate_nnz_distribution(int rank, int npes, block *mat_prop)
{
	long **pe_blocks = block_allocation(npes, npes);
	int i, j, k, l, m, x, y;
	int n=log2(npes);
	//if(!rank) printf("n: %d\n", n);
	for(i=0;i<npes;i++)
		for(j=0;j<npes;j++)
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
								pe_blocks[grid_row+x][grid_col+y] = (long)round(pe_blocks[grid_row+x][grid_col+y]*prob);
								//pe_blocks[grid_row+x][grid_col+y] *= prob;
						}
					}
				}
			}
		}
	}
	long *nnz_dist = (long*)malloc(sizeof(long)*npes);
	memcpy(nnz_dist, pe_blocks[rank], (size_t)npes*sizeof(long));
	free(pe_blocks[0]);
	free(pe_blocks);
	return nnz_dist;
}
long calculate_nnz(long *nnz_dist, int n)
{
	int i;
	long nnz=0;
	for(i=0;i<n;i++)
	{
		//printf("%ld ", pe_blocks[i][j]);
		nnz += nnz_dist[i];
	}
	//printf("nnz: (%ld)\n", nnz);
	return nnz;
}

