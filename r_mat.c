#include<stdio.h>
#include<stdlib.h>
#include<math.h>
typedef struct block_probablity
{
	float a, b, c, d;
	long nnz;
}block;
long** block_allocation(int n_pes)
{
	long **pe_blocks = (long**)malloc(sizeof(long*)*n_pes);
	pe_blocks[0] = (long*)malloc(sizeof(long)*n_pes*n_pes);
	for(int i=1;i<n_pes;i++)
		pe_blocks[i] = pe_blocks[0]+i*n_pes;
	return pe_blocks;
}
void calculate_prob_distribution(long **pe_blocks, int n_pes, block *mat_prop)
{
	int i, j, k, l, m, x, y;
	int n=log2(n_pes);
	for(i=0;i<n_pes;i++)
		for(j=0;j<n_pes;j++)
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
}
void prob_distribution(long **pe_blocks, int n_pes)
{
	int i, j;
	long nnz=0;
	for(i=0;i<n_pes;i++)
	{
		printf("PE: %d -> ", i);
		for(j=0;j<n_pes;j++)
		{
			printf("%ld ", pe_blocks[i][j]);
			nnz += pe_blocks[i][j];
		}
		printf("\n");
	}
	printf("nnz: (%ld)\n", nnz);
}
void row_distribution(int n, long rows_per_pe)
{
	for(int i=0;i<n;i++)
		printf("PE: %d, rows: %ld -> %ld\n", i, rows_per_pe*i, rows_per_pe*(i+1));	
}
int main(int argc, char **argv)
{
	block mat_prob;
	int scaling_factor, n_pes, edge_factor;
	long mat_size, nnz, rows_per_pe;
	scaling_factor = argc >= 2 ? atoi(argv[1]) : 15;
	mat_size = (long)1<<scaling_factor;
	edge_factor = argc >= 3 ? atoi(argv[2]) : 14;
	mat_prob.nnz = mat_size*edge_factor;
	n_pes = argc >= 4 ? atoi(argv[3]) : 8;
	rows_per_pe = mat_size/n_pes;
	if((n_pes & (n_pes-1)) || n_pes == 1)
	{
		printf("number of PEs:%d aren't in powers of 2\n", n_pes);
		exit(1);
	}
	mat_prob.a = argc >=5 ? atof(argv[4]):0.57;
	mat_prob.b = argc >=6 ? atof(argv[5]):0.19;
	mat_prob.c = argc >=7 ? atof(argv[6]):0.19;
	mat_prob.d = 1 - mat_prob.a - mat_prob.b - mat_prob.c;
	printf("matrix distribution a: %.3f, b: %.3f, c: %.3f, d: %.3f\n", mat_prob.a, mat_prob.b, mat_prob.c, mat_prob.d);
	printf("scaling factor: %d, matric_size: %ld, n_pes: %d, edge_factor: %d, nnz: %ld\n", scaling_factor, mat_size, n_pes, edge_factor, mat_prob.nnz);
	long **pe_blocks = block_allocation(n_pes);
	calculate_prob_distribution(pe_blocks, n_pes, &mat_prob);
	prob_distribution(pe_blocks, n_pes);
	row_distribution(n_pes, rows_per_pe);
}
