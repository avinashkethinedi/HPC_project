#pragma once
#include<math.h>
#include<algorithm>
typedef struct block_probablity
{
	float a, b, c, d;
	long nnz;
}block;
typedef struct EDGE
{
	//edge from v to u with weight w
	int v;
	long u;
	float w;
}edge;
typedef struct CSR_MATRIX
{
	int *row_ptr;
	long *col_ptr;
	float *val_ptr;
	long nnz, rows;
}csr_data;
typedef struct SpMV_DATA
{
	csr_data *csr_mat;
	float *multi_vector;
	float *result_vector;
	long mat_size;
}SpMV_data;
typedef struct TIME
{
	double t, avg_t, min_t, max_t;
}time_stats;
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
		//printf("%ld ", pe_blocks[i][j]);
		nnz += nnz_dist[i];
	//printf("nnz: (%ld)\n", nnz);
	return nnz;
}
void stochastic_Kronecker_grpah(edge *edge_array, long pe_nnz, long dim, long start_idx, float *prob, float *c_prob, int *row_nnz_count)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srand(start_idx+rank);
	long mat_dim, row, col, nnz_count=0;
	int p_idx, prow, pcol, i, k;
	float p;
	k = log2(dim);
	while(nnz_count < pe_nnz)
	{
		i=0;
		mat_dim = dim;
		row=0;
		col=0;
		while (i++<k)
		{
			p = (float)rand()/RAND_MAX;
			p_idx = p<c_prob[0]?0:(p<c_prob[1]?1:(p<c_prob[2]?2:3));
			prow = p_idx/2;
			pcol = p_idx%2;
			mat_dim/=2;
			row = row + mat_dim*prow;//optimize with +=
			col = col + mat_dim*pcol;
		}
		edge_array[nnz_count].v = row;
		edge_array[nnz_count].u = start_idx+col;
		edge_array[nnz_count++].w = p;
		row_nnz_count[row]++;
	}
}
//for validation
void print_csr(csr_data *csr_mat)
{
	printf("rows: %ld, nnz: %ld\n", csr_mat->rows, csr_mat->nnz);
	for(long i=0;i<csr_mat->rows;i++)
	{
		printf("%ld: ",i);
		for(long j=csr_mat->row_ptr[i];j<csr_mat->row_ptr[i+1];j++)
			printf("%ld ", csr_mat->col_ptr[j]);
		printf("\n");
	}
}
//for validation
void validate_csr(csr_data *csr_mat)
{
	for(long i=0;i<csr_mat->rows;i++)
		for(long j=csr_mat->row_ptr[i]+1;j<csr_mat->row_ptr[i+1];j++)
			if(csr_mat->col_ptr[j] < csr_mat->col_ptr[j-1])
			{
				printf("data incorrect r:%ld, c:%ld\n", i, j);
				return;
			}
	printf("data validated\n");
}
csr_data* create_csr_data(edge *edge_array, int *row_nnz_count, long pe_nnz, long rows_per_pe)
{
	size_t n = (size_t)sizeof(csr_data)+(rows_per_pe+1)*sizeof(int) + pe_nnz*(sizeof(long) + sizeof(float));
	csr_data *csr_mat = (csr_data*)malloc(n);
	csr_mat->nnz = pe_nnz;
	csr_mat->rows = rows_per_pe;
	long i, idx;
	//printf("SIZE:%ld B", n);
	if(csr_mat == NULL)
	{
		printf("csr memory allocation failed (%ld B)\n", n);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	n = sizeof(csr_data);
	csr_mat->row_ptr = (int*)((char*)csr_mat + n);
	n+=((rows_per_pe+1)*sizeof(int));
	csr_mat->col_ptr = (long*)((char*)csr_mat + n);
	n+=(pe_nnz*sizeof(long));
	csr_mat->val_ptr = (float*)((char*)csr_mat + n);

	//row_ptr
	csr_mat->row_ptr[0] = 0;
	for(i=1;i<=rows_per_pe;i++)
		csr_mat->row_ptr[i] = csr_mat->row_ptr[i-1] + row_nnz_count[i-1];
	for(i=0;i<pe_nnz;i++)
	{
		idx = edge_array[i].v;
		row_nnz_count[idx]--;
		idx = csr_mat->row_ptr[idx] + row_nnz_count[idx];
		csr_mat->col_ptr[idx] = edge_array[i].u;
		csr_mat->val_ptr[idx] = edge_array[i].w;
	}
	for(i=0;i<rows_per_pe;i++)
		std::sort(&csr_mat->col_ptr[csr_mat->row_ptr[i]], &csr_mat->col_ptr[csr_mat->row_ptr[i+1]]);
	return csr_mat;
}
csr_data* create_matrix_data(long *nnz_dist, long pe_nnz, long rows_per_pe, int npes, block *mat_prob)
{
	edge *edge_array = (edge*)malloc((size_t)pe_nnz*sizeof(edge));
	int i, *row_nnz_count = (int*)calloc(rows_per_pe, sizeof(int));
	csr_data *csr_mat;
	long block_nnz, start_idx, idx = 0;
	float prob[] = {mat_prob->a, mat_prob->b, mat_prob->c, mat_prob->d}, c_prob[4];
	c_prob[0] = prob[0];
	for(i=1;i<4;i++)
		c_prob[i] = prob[i] + c_prob[i-1];
	for(i=0;i<npes;i++)
	{
		block_nnz = nnz_dist[npes-1-i];
		start_idx = rows_per_pe*(npes-1-i);
		stochastic_Kronecker_grpah(&edge_array[idx], block_nnz, rows_per_pe, start_idx, prob, c_prob, row_nnz_count);
		idx += block_nnz;
	}
	csr_mat = create_csr_data(edge_array, row_nnz_count, pe_nnz, rows_per_pe);
	//free data;
	free(edge_array);
	free(row_nnz_count);
	return csr_mat;
}
void SpMV_kernel(SpMV_data *data)
{
	long i, j;
	for(i=0;i<data->csr_mat->rows;i++)
	{
		data->result_vector[i] = 0;
		for(j=data->csr_mat->row_ptr[i];i<data->csr_mat->row_ptr[i+1];i++)
			data->result_vector[i] += data->csr_mat->val_ptr[j]*data->multi_vector[data->csr_mat->col_ptr[j]];
	}
}
void iterative_SpMV(SpMV_data *data, int iterations)
{
	for(int i=0;i<iterations;i++)
	{
		SpMV_kernel(data);
		//TODO: compress and broadcast
		//TODO: check accuracy
		float *temp = data->multi_vector;
		data->multi_vector = data->result_vector;
		data->result_vector = temp;
	}
}