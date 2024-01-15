#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>

void Get_input(int my_rank,int *m,int *n)
{
	if(my_rank==0){
		printf("Please enter m,n:\n");
		scanf("%d %d",m,n);
	}
	MPI_Bcast(m,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);
}




//得到矩阵
void Get_matrix(int n, int m, double *local_matrix, int local_m, int my_rank)
{
    double *A;
    if (!my_rank)
    {
        A = (double *)malloc(m * n * sizeof(double));
        printf("Please enter the matrix:\n");
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                scanf("%lf", &A[i * n + j]);
    }
	//MPI_Scatter函数将矩阵分发出去
    MPI_Scatter(A, local_m * n, MPI_DOUBLE, local_matrix, local_m * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

//打印矩阵
//MPI_Gather函数将local_matrix从各个进程聚集到0号进程输出
void Print_matrix(int my_rank,int n,int m,int local_m,double *local_matrix,MPI_Comm comm)
{
	double *matrix = NULL;
	int i,j;
	if(my_rank==0)
	{
		matrix = malloc(m*n*sizeof(double));
		MPI_Gather(local_matrix,local_m*n,MPI_DOUBLE,matrix,local_m*n,MPI_DOUBLE,0,comm);
		printf("The matrix is:\n");
		for(i=0;i<m;++i)
		{
			for(j=0;j<n;++j)
			{
				printf("%f ",matrix[i*n+j]);
			}
			printf("\n");
		}
		free(matrix);
	}
	else{
		MPI_Gather(local_matrix,local_m*n,MPI_DOUBLE,matrix,local_m*n,MPI_DOUBLE,0,comm);
	}
}

//得到向量并分发
void Get_vector(int my_rank,int n,int local_n,double *local_vector,MPI_Comm comm)
{
	double *vector = NULL;
	int i;
	if(my_rank==0)
	{
		vector=(double *)malloc(n*sizeof(double));
		printf("Please enter the vector:\n");
		for(i=0;i<n;i++)
		{
			scanf("%lf",&vector[i]);
		}
	}
	printf("\n");
	MPI_Scatter(vector,local_n,MPI_DOUBLE,local_vector,local_n,MPI_DOUBLE,0,comm);
}



//聚合向量到0号进程并且输出
void Print_vector(int my_rank,int n,int local_n,double *local_vector,MPI_Comm comm)
{
	double *vector = NULL;
	int i,j;
	if(my_rank==0)
	{
		vector = malloc(n*sizeof(double));
		MPI_Gather(local_vector,local_n,MPI_DOUBLE,vector,local_n,MPI_DOUBLE,0,comm);
		printf("The vector is:\n");
		for(i=0;i<n;i++){
			printf("%f ",vector[i]);
		}
		printf("\n");
		free(vector);
	}
	else{
		MPI_Gather(local_vector,local_n,MPI_DOUBLE,vector,local_n,MPI_DOUBLE,0,comm);
	}
}


//实现矩阵乘法

void Mat_vect_mult(double *local_matrix,double *local_vector,double *local_y,int local_m,int n,int local_n,MPI_Comm comm)
{
	int local_i,j;
	double *x;
	
	x=malloc(n*sizeof(double));
	//将向量聚合到所有进程，MPI_Allgather和MPI_Gather的区别就在于Allgather的所
	//有进程都会知道你聚合到的的向量，相当于聚合到0号进程之后又bcast广播了一次
	MPI_Allgather(local_vector,local_n,MPI_DOUBLE,x,local_n,MPI_DOUBLE,comm);
	
	for(local_i=0;local_i<local_m;local_i++)
	{
		local_y[local_i]=0.0;
		for(j=0;j<n;j++)
		{
			local_y[local_i]+=local_matrix[local_i*n+j]*x[j];
		}
	}
	free(x);
}

//打印结果
void Print_y(int my_rank,double *local_y,int m,int local_m,MPI_Comm comm)
{
	double *y=NULL;
	int i;
	if(my_rank==0){
		y=malloc(m*sizeof(double));
		MPI_Gather(local_y,local_m,MPI_DOUBLE,y,local_m,MPI_DOUBLE,0,comm);
		printf("The vector y is:\n");
		for(i=0;i<m;i++)
		{
			printf("%lf ",y[i]);
		}
		printf("\n");
		free(y);
	}
	else{
		MPI_Gather(local_y,local_m,MPI_DOUBLE,y,local_m,MPI_DOUBLE,0,comm);
	}
}

int main()
{
	int comm_sz,my_rank,i;
	int m,n,local_m,local_n;
	double *local_matrix,*local_vector;
	double *local_y;
	
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	Get_input(my_rank,&m,&n);
	local_m=m/comm_sz;
	local_n=n/comm_sz;
	local_matrix=(double *)malloc(local_m*n*sizeof(double));
	local_vector=(double *)malloc(local_n*sizeof(double));
	local_y=(double *)malloc(local_m*sizeof(double));

	Get_matrix(n,m,local_matrix,local_m,my_rank);
	Print_matrix(my_rank,n,m,local_m,local_matrix,MPI_COMM_WORLD);
	Get_vector(my_rank,n,local_n,local_vector,MPI_COMM_WORLD);
	Print_vector(my_rank,n,local_n,local_vector,MPI_COMM_WORLD);
	Mat_vect_mult(local_matrix,local_vector,local_y,local_m,n,local_n,MPI_COMM_WORLD);
	Print_y(my_rank,local_y,m,local_m,MPI_COMM_WORLD);
	
	MPI_Finalize();
	return 0;
}
