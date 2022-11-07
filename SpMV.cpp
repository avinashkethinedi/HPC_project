//#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<mpi.h>
#include"r_mat.h"
void print_help();
void set_parameters(int, char**, int*, int*, block*, long*);
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int scaling_factor = 7, edge_factor = 27;
	long mat_size, rows_per_pe, pe_nnz, nnz=0, *nnz_dist;
	block mat_prob = {0.25, 0.25, 0.25, 0.25, 1};//{0.57, 0.19, 0.19, 0.05, 1};
	set_parameters(argc, argv, &scaling_factor, &edge_factor, &mat_prob, &mat_size);
	int rank, npes, mat_blocks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	rows_per_pe = mat_size/npes;
	mat_blocks = mat_size/rows_per_pe;
	//probability distribution
	nnz_dist = calculate_nnz_distribution(rank, npes, &mat_prob);
	pe_nnz = calculate_nnz(nnz_dist, npes);
	//matrix creation
	csr_data *csr_mat = create_matrix_data(nnz_dist, pe_nnz, rows_per_pe, npes, &mat_prob);
	MPI_Reduce(&pe_nnz, &nnz, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	//free data
	free(nnz_dist);
	free(csr_mat);
	MPI_Barrier(MPI_COMM_WORLD);
	if(!rank)
	{
		printf("scaling factor: %d\n", scaling_factor);
		printf("edge factor: %d (%0.3lf)\n", edge_factor, ((double)nnz/mat_size));
		printf("probabilities: %0.3f, %0.3f, %0.3f, %0.3f\n", mat_prob.a, mat_prob.b, mat_prob.c, mat_prob.d);
		printf("matrix size: %ld\n", mat_size);
		printf("nnz: %ld (%ld)\n", nnz, mat_prob.nnz);
		printf("rows per PE: %ld\n", rows_per_pe);
		printf("mat blocks: %d\n",mat_blocks);
	}
	MPI_Finalize();
	return EXIT_SUCCESS;
}

void print_help()
{
	printf("scaling factor -s (15)\n");
	printf("edge factor -e (20)\n");
	printf("probability values: -a(0.57) -b(0.19) -c(0.19)\n");
	printf("help -h");
	exit(1);
}
void set_parameters(int argc, char **argv, int *scaling_factor, int *edge_factor, block *mat_prob, long *mat_size)
{
	int opt;
	while((opt = getopt(argc, argv, ":s:e:a:b:c:h:")) != -1)
	{
		switch(opt)
		{
			case 's':
				*scaling_factor=atoi(optarg);
				if(*scaling_factor < 1)
				{
					printf("invalid value for -s (it requires positive integer)\n");
					print_help();
				}
				break;
			case 'e':
				*edge_factor=atoi(optarg);
				if(*edge_factor < 1)
				{
					printf("invalid value for -e (it requires positive integer)\n");
					print_help();
				}
				break;
			case 'a':
				mat_prob->a = atof(optarg);
				if(mat_prob->a > 1 || mat_prob->a < 0)
				{
					printf("invalid value for -a (it requires positive value [0,1])\n");
					print_help();
				}
				break;
			case 'b':
				mat_prob->b = atof(optarg);
				if(mat_prob->b > 1 || mat_prob->b < 0)
				{
					printf("invalid value for -b (it requires positive value [0,1])\n");
					print_help();
				}
				break;
			case 'c':
				mat_prob->c = atof(optarg);
				if(mat_prob->c > 1 || mat_prob->c < 0)
				{
					printf("invalid value for -c (it requires positive value [0,1])\n");
					print_help();
				}
				break;
			case 'h':
				print_help();
				break;
			default:
				if(optopt == 's')
					*scaling_factor=15;
				else if(optopt == 'e')
					*edge_factor = 14;
				else if(optopt == 'a')
					mat_prob->a = 0.25;
				else if(optopt == 'b')
					mat_prob->b = 0.25;
				else if(optopt == 'c')
					mat_prob->c = 0.25;
				else
				{
					print_help();
					exit(1);
				}
		}
	}
	mat_prob->d = 1 - (mat_prob->a + mat_prob->b + mat_prob->c);
	if(mat_prob->d > 1 || mat_prob->d < 0)
	{
		printf("sum of probabilities (a: %0.3f, b: %0.3f, c: %0.3f, d: %0.3f) != 1\n", mat_prob->a, mat_prob->b, mat_prob->c, mat_prob->d);
		print_help();
	}
	*mat_size = (long)1<<(*scaling_factor);
	mat_prob->nnz = (*mat_size)*(*edge_factor);
}