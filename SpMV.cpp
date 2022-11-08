//#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<mpi.h>
#include"r_mat.h"
void print_help();
void set_parameters(int, char**, int*, int*, block*, long*, int*);
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int scaling_factor = 7, edge_factor = 27, iterations = 10;
	long mat_size, rows_per_pe, pe_nnz, nnz=0, *nnz_dist;
	time_stats graph_time = {0, 0, 0, 0}, SpMV_time = {0, 0, 0, 0};
	block mat_prob = {0.25, 0.25, 0.25, 0.25, 1};
	set_parameters(argc, argv, &scaling_factor, &edge_factor, &mat_prob, &mat_size, &iterations);
	int rank, npes, mat_blocks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	rows_per_pe = mat_size/npes;
	mat_blocks = mat_size/rows_per_pe;
	//probability distribution
	nnz_dist = calculate_nnz_distribution(rank, npes, &mat_prob);
	pe_nnz = calculate_nnz(nnz_dist, npes);
	//matrix creation
	graph_time.t = MPI_Wtime();
	csr_data *csr_mat = create_matrix_data(nnz_dist, pe_nnz, rows_per_pe, npes, &mat_prob);
	graph_time.t = MPI_Wtime() - graph_time.t;
	//SpMV kernel
	SpMV_data data;
	data.csr_mat = csr_mat;
	data.multi_vector = (float*)malloc(mat_size*2*sizeof(float));
	data.result_vector = &data.multi_vector[mat_size];
	//initilise multiplication vector
	for(long i=0;i<mat_size;i++) data.multi_vector[i] = rand()/RAND_MAX;
	//iterative SpMV
	SpMV_time.t = MPI_Wtime();
	iterative_SpMV(&data, iterations);
	SpMV_time.t = MPI_Wtime() - SpMV_time.t;
	//free data
	free(nnz_dist);
	free(csr_mat);
	free(data.multi_vector);
	MPI_Barrier(MPI_COMM_WORLD);
	//stats
	MPI_Reduce(&pe_nnz, &nnz, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&graph_time.t, &graph_time.min_t, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&graph_time.t, &graph_time.max_t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&graph_time.t, &graph_time.avg_t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&SpMV_time.t, &SpMV_time.min_t, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&SpMV_time.t, &SpMV_time.max_t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&SpMV_time.t, &SpMV_time.avg_t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	graph_time.avg_t /= npes;
	SpMV_time.avg_t /= npes;
	MPI_Barrier(MPI_COMM_WORLD);
	if(!rank)
	{
		printf("scaling factor: %d\n", scaling_factor);
		printf("edge factor: %d (%0.3lf)  iterations: %d\n", edge_factor, ((double)nnz/mat_size), iterations);
		printf("probabilities: %0.3f, %0.3f, %0.3f, %0.3f\n", mat_prob.a, mat_prob.b, mat_prob.c, mat_prob.d);
		printf("matrix size: %ld\n", mat_size);
		printf("nnz: %ld (%ld)\n", nnz, mat_prob.nnz);
		printf("rows per PE: %ld\n", rows_per_pe);
		printf("mat blocks: %d\n",mat_blocks);
		printf("graph time: avg: %0.5lf s,  min: %0.5lf s,  max: %0.5lf s\n", graph_time.avg_t, graph_time.min_t, graph_time.max_t);
		printf("SpMV time: avg: %0.5lf s,  min: %0.5lf s,  max: %0.5lf s\n", SpMV_time.avg_t, SpMV_time.min_t, SpMV_time.max_t);
	}
	MPI_Finalize();
	return EXIT_SUCCESS;
}

void print_help()
{
	printf("scaling factor -s (15)\n");
	printf("edge factor -e (20)\n");
	printf("probability values: -a(0.57) -b(0.19) -c(0.19)\n");
	printf("number of iterations -i (10)\n");
	printf("help -h\n");
	exit(1);
}
void set_parameters(int argc, char **argv, int *scaling_factor, int *edge_factor, block *mat_prob, long *mat_size, int *iterations)
{
	int opt;
	while((opt = getopt(argc, argv, ":s:e:a:b:c:h:i:")) != -1)
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
			case 'i':
				*iterations=atoi(optarg);
				if(*iterations < 5 || *iterations > 50)
				{
					printf("invalid value for -i (range 5 - 50)\n");
					print_help();
				}
				break;
			case 'h':
				print_help();
				break;
			default:
				if(optopt == 's')
					*scaling_factor = 15;
				else if(optopt == 'e')
					*edge_factor = 14;
				else if(optopt == 'a')
					mat_prob->a = 0.25;
				else if(optopt == 'b')
					mat_prob->b = 0.25;
				else if(optopt == 'c')
					mat_prob->c = 0.25;
				else if(optopt == 'i')
					*iterations = 10;
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