/*
********************************************************************************
MULTICORE RADIX CLUSTERING FOR SKEWED DATA USING TASK-BASED PARALLELIZATION:
********************************************************************************
This C program implements a version of radix clustering for skewed data using thread-level parallelism
for handling new passes. This is done in such a way that requires no synchronization between threads.

Since the parallelization alternative studied here arises from the logic of the original program, 
without considering special optimizations, we consider this could be called a "naive" parallelization.
Also, because it is task-based and relying on the concept of adjustable parallelism, we call this
adjustable parallel radix clustering.

This version could serve as a framework, upon which optimization options can be built (i.e.
the use of an active pool of threads, the use of helper threads for pre-fetching, configuration of thread affinity, 
vectorized execution, among others).

It is coded after Boncz et al serial Radix join implementation,* expanding upon it by allowing input 
(a series of unsigned integers), skewness (variable cluster sizes) and the framework for parallelization.

The handling of skewness is achieved by performing an extra pass in which the number of tuples per cluster
is counted by doing a prefix sum over this result. Because of this our implementation would take about 1.5 
the time of the Boncz et al code, for serial cases.


*Available in http://homepages.cwi.nl/~boncz/mimuw/. 

********************************************************************************
CALCULATION OF MAXIMUM POSSIBLE TASK BASED PARALLELIZATION
********************************************************************************
The main novelty introduced in our proposal, consists on a scheme for the calculation of the maximum 
possible task based parallelization given the maximum number of threads which can be used, the number 
of passes and the considered bits per pass. Its worth to point out that this only 
refers to a naive parallelization, to be done at the moment of calling for new
passes over existing clusters.

This calculation is done before the execution of the algorithm, its results are stored
in an arraz that contains the number of new threads that any thread from a given pass can 
use for parallelizing the processing of the new clusters.

By default, the first pass will always be done by one thread alone.

As an example: 
Suppose that the first pass creates 8 clusters, and the function is to be called again over each of those clusters for the second pass.

Here would be the results of this array (called Threads_per_pass) given the available number of threads:
Number of Threads= 1 Threads_per_pass={1,0} 	=>No parallelization possible for calling the second pass
Number of Threads= 2 Threads_per_pass={1,1} 	=>The second pass will be called in a parallel way.
						For the existing 8 clusters, the second pass will be called in the existing thread for 4 clusters
						and in a new thread for the remaining clusters.
Number of Threads= 3 Threads_per_pass={1,1}	=>Same as the previous case.
Number of Threads= 4 Threads_per_pass={1,3}	=>The second pass will be called in a prallel way.
						For the existing 8 clusters, the second pass will be called in the existing thread for 2 clusters
						and in 3 new threads (one for each 2 clusters)
Number of Threads= 5 Threads_per_pass={1,3}	=>Same as the previous case.
Number of Threads= 6 Threads_per_pass={1,3}	=>Same as the previous case.
Number of Threads= 7 Threads_per_pass={1,3}	=>Same as the previous case.
Number of Threads= 8 Threads_per_pass={1,7}	=>Per cluster parallelization possible. The second pass will be called in a different thread over each 
						of the existing 8 clusters.
Number of Threads= 9 or more Threads_per_pass={1,7} =>Same as the previous case.


********************************************************************************
OBSERVATIONS
********************************************************************************
- We considered implementing to allow tuples of type char, but conversion took too much time. 
An alternative implementation for chars and floats could be easily done.

********************************************************************************
INPUT
********************************************************************************
To be passed to the program as a file, with the following format:
On the first line: Number of repeats
On the second line: Maximum number of threads.
On the third line: Number of passes
On the forth line: Number of bits per pass, separated by spaces or commas
On the fifth line: Expected affinity of threads between former passes. 1 if Affinity, 0 if not. 
On the sixth line: Number of tuples
On the seventh line, onwards: Values of the tuples (unsigned integers).

********************************************************************************
OUTPUT
********************************************************************************
This prints the resulting clustering of the tuples.

********************************************************************************
FUNCTIONS INCLUDED
********************************************************************************
- uint timer(void): Returns the current time in miliseconds.
- int load_Data (char [], uint&, uint&, uint*&, uint&, uint&, uint* &): 
	Loads the data from the file in the given format. Returns 1 if successful, 0 if invalid input, -1 if missing tuples, -2 if too many tuples.
- void prefix_sum(uint*&, uint): Utility function for calculating a prefix sum over an array.
- void radix_cluster_for_skewed_data(uint, uint, uint, uint[], uint*, uint*, uint *): 
	Following Boncz, if the number of passes is even, the result will be stored in the second to last parameter, otherwise in the third to last.
- int main(int, char*[])

*******************************************************************************
KNOWN BUGS
********************************************************************************
It has been observed that the load data function sometimes returns missing tuples 
when all that is missing is an endline character at the end of the input file. 
This uncommon bug has not been fixed in the current implementation, but can be 
solved (in the meantime) by adding an extra empty line to the input file. 

Affinity setting was suggested and its implementation is proposed here. However the
results for this have been sub-optimal, which lead us to think that something maybe
from faulty design/implementation, from the pthreads library, or even from the OS
prevent us from seeing the result of this affinity setting. This needs to be evaluated
further.
********************************************************************************
OvGU Magdeburg June 2014.

Authors: 
Gabriel Campero gabrielcampero@gmail.com
Angel Yordanov angellyordanov@gmail.com
Svetoslav Ermakov slav.ermakov@gmail.com

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> //For logarithm calculations related to maximum naive parallelization
#include <sys/time.h> //For testing
#include <pthread.h>
#include <unistd.h> //For determining number of cores
#define MAXLINESIZE 128  /* or other suitable maximum line size for reading from file, and maximum tuple length as strings*/
typedef unsigned int uint;

//GLOBAL VARIABLES.
uint* threads_per_pass; 
bool*  thread_affinity;
uint total_passes;

struct thread_data {
	uint number_tuples;
	uint passes;
	uint total_bits; 
	uint* bits_per_pass;
	uint* result0;
	uint* result1;
	uint* tuples;
};

struct block_thread_data {
	uint* number_tuples;
	uint* passes;
	uint* total_bits; 
	uint** bits_per_pass;
	uint** result0;
	uint** result1;
	uint** tuples;
	uint n;
};



uint timer(void) {
        struct timeval tp;
        gettimeofday(&tp, NULL);
        return tp.tv_sec * 1000 + tp.tv_usec / 1000;
}

/*The following function loads the parameters for the clustering from a file, with the following format:
On the first line: Number of repeats
On the second line: Maximum number of threads.
On the third line: Number of passes
On the forth line: Number of bits per pass, separated by spaces or commas
On the fifth line: Expected affinity of threads between former passes. 1 if Affinity, 0 if not. 
On the sixth line: Number of tuples
On the seventh line, onwards: Values of the tuples (unsigned integers).
*/
int load_Data (char filename [], uint& nthreads, uint& repeats, uint& passes, uint*& bits_per_pass, uint& total_bits, uint& number_tuples, uint* &tuples){
	FILE *file = fopen ( filename, "r" );
	if ( file != NULL )
		{
			char line [ MAXLINESIZE ];


			fgets ( line, sizeof line, file );
			if (atoi(line)<=0)
				return 0;
			repeats=atoi(line);

			fgets ( line, sizeof line, file );
			if (atoi(line)<=0)
				return 0;
			nthreads=atoi(line);			

			fgets ( line, sizeof line, file );
			if (atoi(line)<=0)
				return 0;
			passes=atoi(line);

			bits_per_pass=(uint*)malloc(sizeof(uint)*passes);
			thread_affinity=(bool*)malloc(passes*sizeof(bool));
			fgets ( line, sizeof line, file );
			total_bits=0;
			char* ptr;
			if (strtol (line, &ptr, 10)<=0)
				return 0;
			bits_per_pass[0]=strtol (line, &ptr, 10);
			total_bits+=bits_per_pass[0];
			for (uint i = 1; i < passes; i++) 
			{
				int aux= strtol (ptr, &ptr, 10);
				if (aux<=0)
					return 0;
				bits_per_pass[i]=aux;
				total_bits+=bits_per_pass[i];
			}

			fgets ( line, sizeof line, file );
			if (strtol (line, &ptr, 10)<0)
				return 0;
			thread_affinity[0]=strtol (line, &ptr, 10)==0?false:true;
			for (uint i = 1; i < passes; i++) 
			{
				int aux= strtol (ptr, &ptr, 10);
				if (aux<0)
					return 0;
				thread_affinity[i]=aux==0?false:true;
			}

			fgets ( line, sizeof line, file );
			if (atoi(line)<=0)
				return 0;
			number_tuples=atoi(line);
			tuples=(uint*)malloc(sizeof(uint)*number_tuples);
			uint currentTuple=-1;
			while ( fgets ( line, sizeof line, file ) != NULL ) 
			{
				currentTuple++;
				if (currentTuple<number_tuples){
					tuples[currentTuple]=atoi(line);
				}
				if (currentTuple>number_tuples)
					return -2;

			}
			if (currentTuple!=number_tuples)
				return -1;
			fclose ( file );
		}
	else
		{
			perror ( filename );
			return 0;
		}
	return 1;

}
/*Utility function. It is used to determine the lower borders of each cluster with skewed data.*/
void prefix_sum(uint*& array, uint size){
	uint tmp1=0;	
	uint tmp2=0;
	for (int i=0; i<size; i++){
		tmp2+=array[i];
		array[i]=tmp1;
		tmp1=tmp2;
	}
}

void* _parallel_radix_clustering_naive(void* threadarg); 

/*Recursive version of radix clustering for skewed data (when clusters from a given pass have different number of tuples).*/
void radix_cluster_for_skewed_data(uint number_tuples, uint passes, uint total_bits, uint bits_per_pass[], uint* result0, uint* result1, uint * tuples) {
	uint i, clustersize = number_tuples >> bits_per_pass[0], nclusters = 1 << bits_per_pass[0], mask = nclusters-1, lowbits = total_bits-bits_per_pass[0];
	uint s;

	uint *borders = (uint*) alloca(nclusters*sizeof(uint)); 
	for(i=0; i<nclusters; i++) 
		borders[i] = 0;
	for(i=0; i < number_tuples; i++) 
		borders[( tuples[i] >> lowbits) & mask]++; //A first pass is performed, counting the number of tuples per cluster.
	prefix_sum(borders, nclusters); //By a prefix-sum we calculate the lower borders of clusters.

	/* do the partitioning work of this pass <===Same as in Boncz. */
	for(i=0; i < number_tuples; i++) {
		uint clustnr = ( tuples[i] >> lowbits) & mask;
		result0[borders[clustnr]++]= tuples[i];	
	/*Note that the border array is changed so instead of storing lower borders, it now stores the position of the last tuple in each cluster.*/
	}
	if (passes > 1){ /* more passes to do? call recursively ... Similar as in Boncz, but changed to handle uneven cluster borders for skewed data.*/
		s=0; //s= number of tuples below given cluster	
		for(i=0; i<nclusters; i++){
			radix_cluster_for_skewed_data(borders[i]-s, passes-1, lowbits, bits_per_pass+1, result1+s, result0+s, result0+s);
			s=borders[i];
		}
	}

}


/*Function for calling the serial/recursive version in a per-block way.*/
void* _per_block_parallel_radix_clustering_naive(void* threadarg){
	struct block_thread_data *my_data= (struct block_thread_data *) threadarg;
	for (uint i=0; i<my_data->n; i++){
		radix_cluster_for_skewed_data(my_data->number_tuples[i], my_data->passes[i], my_data->total_bits[i], my_data->bits_per_pass[i], my_data->result0[i], my_data->result1[i], my_data->tuples[i]);
	}
}

/*Function implementing parallel radix clustering, given the pre-computation of threads per pass.*/
void* _parallel_radix_clustering_naive(void* threadarg) {
	struct thread_data *my_data= (struct thread_data *) threadarg;
	uint i, clustersize = my_data->number_tuples >> my_data->bits_per_pass[0], nclusters = 1 << my_data->bits_per_pass[0], mask = nclusters-1, lowbits = my_data->total_bits-my_data->bits_per_pass[0];
	uint s;

	uint *borders = (uint*) alloca(nclusters*sizeof(uint)); 
	for(i=0; i<nclusters; i++) 
		borders[i] = 0;
	for(i=0; i < my_data->number_tuples; i++) 
		borders[( my_data->tuples[i] >> lowbits) & mask]++; //A first pass is performed, counting the number of tuples per cluster.
	prefix_sum(borders, nclusters); //By a prefix-sum we calculate the lower borders of clusters.

	/* do the partitioning work of this pass <===Same as in Boncz. */
	for(i=0; i < my_data->number_tuples; i++) {
		uint clustnr = ( my_data->tuples[i] >> lowbits) & mask;
		my_data->result0[borders[clustnr]++]= my_data->tuples[i];	
	/*Note that the border array is changed so instead of storing lower borders, it now stores the position of the last tuple in each cluster.*/
	}
	if (my_data->passes > 1){ /* more passes to do?... Similar as in Boncz, but changed to handle uneven cluster borders for skewed data. And naive parallelism handling.*/
		s=0; //s= number of tuples below given cluster	
		uint nthreads=threads_per_pass[1+total_passes-my_data->passes];		
		if(nthreads==0){		//call recursively
			for(i=0; i<nclusters; i++){
				radix_cluster_for_skewed_data(borders[i]-s, my_data->passes-1, lowbits, my_data->bits_per_pass+1, my_data->result1+s, my_data->result0+s, my_data->result0+s); 
/*Note, this was called with another function, thus increasing instructions count, but reducing number of operations from avoiding packing and unpacking of data.*/
				s=borders[i];
			}
		}
		else if(nthreads==nclusters-1){ //call with a new thread per cluster
			pthread_t threads[nthreads];
			pthread_attr_t attr;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			if(thread_affinity[total_passes-my_data->passes]){//setaffinity
				cpu_set_t currentCPU;
				CPU_ZERO(&currentCPU);
				CPU_SET(sched_getcpu(), &currentCPU);
				pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &currentCPU);
			}
			int rc; 
			struct thread_data my_data2[nthreads+1];
			for(i=0; i<nthreads; i++){
				my_data2[i].number_tuples=borders[i]-s;
				my_data2[i].passes=my_data->passes-1;
				my_data2[i].total_bits=lowbits;
				my_data2[i].bits_per_pass=my_data->bits_per_pass+1;
				my_data2[i].result0=my_data->result1+s;
				my_data2[i].result1=my_data->result0+s;
				my_data2[i].tuples=my_data->result0+s;
				rc = pthread_create(&threads[i], &attr, _parallel_radix_clustering_naive, (void *)&my_data2[i]); 
			        if (rc){
					printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
			      	}
				s=borders[i];
			}
			my_data2[i].number_tuples=borders[i]-s;
			my_data2[i].passes=my_data->passes-1;
			my_data2[i].total_bits=lowbits;
			my_data2[i].bits_per_pass=my_data->bits_per_pass+1;
			my_data2[i].result0=my_data->result1+s;
			my_data2[i].result1=my_data->result0+s;
			my_data2[i].tuples=my_data->result0+s;
			_parallel_radix_clustering_naive((void*)&my_data2[i]);
			/* Free attribute and wait for the other threads */
			pthread_attr_destroy(&attr);
			//Do the joins...
			void* status;
			for(i=0; i<nthreads; i++){
			      rc = pthread_join(threads[i], &status);
			      if (rc) {
			         printf("ERROR; return code from pthread_join() is %d\n", rc);
			         exit(-1);
         		      }
			}
	

		}
		else {	
			pthread_t threads[nthreads];
			pthread_attr_t attr;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			if(thread_affinity[total_passes-my_data->passes]){//setaffinity
				cpu_set_t currentCPU;
				CPU_ZERO(&currentCPU);
				CPU_SET(sched_getcpu(), &currentCPU);
				pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &currentCPU);
			}
			int rc; 
			uint n=nclusters/(nthreads+1);
			struct block_thread_data my_data2[nthreads+1];
			uint j;
			for(i=0; i<nthreads+1; i++){
				my_data2[i].number_tuples=(uint*)malloc(n*sizeof(uint));
				my_data2[i].passes=(uint*)malloc(n*sizeof(uint));
				my_data2[i].total_bits=(uint*)malloc(n*sizeof(uint));
				my_data2[i].bits_per_pass=(uint**)malloc(n*sizeof(uint*));
				my_data2[i].result0=(uint**)malloc(n*sizeof(uint*));
				my_data2[i].result1=(uint**)malloc(n*sizeof(uint*));
				my_data2[i].tuples=(uint**)malloc(n*sizeof(uint*));
				my_data2[i].n=n;
				for (j=0; j<n; j++){
					my_data2[i].number_tuples[j]=borders[i*n+j]-s;
					my_data2[i].passes[j]=my_data->passes-1;
					my_data2[i].total_bits[j]=lowbits;
					my_data2[i].bits_per_pass[j]=my_data->bits_per_pass+1;
					my_data2[i].result0[j]=my_data->result1+s;
					my_data2[i].result1[j]=my_data->result0+s;
					my_data2[i].tuples[j]=my_data->result0+s;
					s=borders[i*n+j];
				}

				if(i!=nthreads){
					rc = pthread_create(&threads[i], &attr, _per_block_parallel_radix_clustering_naive, (void *)&my_data2[i]); 
				        if (rc){
						printf("ERROR; return code from pthread_create() is %d\n", rc);
						exit(-1);
				      	}
				}
				else{ //Recursive call
					_per_block_parallel_radix_clustering_naive((void *)&my_data2[i]);
				}

			}
			/* Free attribute and wait for the other threads */
			pthread_attr_destroy(&attr);
			//Do the joins...
			void* status;
			for(i=0; i<nthreads; i++){
			      rc = pthread_join(threads[i], &status);
			      if (rc) {
			         printf("ERROR; return code from pthread_join() is %d\n", rc);
			         exit(-1);
         		      }
			}

		}
	}


}

/*Launching function for parallel radix clustering for skewed data (when clusters from a given pass have different number of tuples).
It only performs initiallization calculations and data packing, and then calls: _parallel_radix_clustering_naive*/
void* parallel_radix_clustering_naive(uint nthreads, uint number_tuples, uint passes, uint total_bits, uint bits_per_pass[], uint* result0, uint* result1, uint * tuples) {
	
	//FIRST WE CALCULATE THE BEST THEORETICAL PARALLELIZATION

	/*Initialization*/
	threads_per_pass=(uint*)malloc(passes*sizeof(uint));
	threads_per_pass[0]=1;
	uint currentThreads=1;
	uint unusedThreads=nthreads-1;
	uint clusters_per_pass;
	

	/*The code itself for deciding the proper parallelization. */
	uint i;
	for (i=1; i<passes && unusedThreads>0; i++){
		clusters_per_pass=1<<bits_per_pass[i-1];
		int val= (int) floor(       log2((double)(unusedThreads/currentThreads) +1)         ); //Could this also be written using binary shifts?
		threads_per_pass[i]=(val==0)?0:(((1<<val)-1)>=clusters_per_pass-1)?clusters_per_pass-1:(1<<val)-1;
		if(threads_per_pass[i]==clusters_per_pass-1){ //Per cluster parallelization achieved (i.e. one thread per each cluster)
			unusedThreads-=	threads_per_pass[i]*currentThreads;
			currentThreads+=(threads_per_pass[i]*currentThreads);
		}
		else{ //No parallelization or Per block parallelization (i.e. one thread will handle a number of clusters, this number will be a power of 2). 
		//If there are still unused threads, they won't be usable in the next passes.
			unusedThreads=0;
			currentThreads+=(threads_per_pass[i]*currentThreads);
		}
	}
	while (i<passes){
		threads_per_pass[i]=0;
		i++;
	}

	struct thread_data my_data;
	my_data.number_tuples=number_tuples;
	my_data.passes=passes;
	my_data.total_bits=total_bits;
	my_data.bits_per_pass=bits_per_pass;
	my_data.result0=result0;
	my_data.result1=result1;
	my_data.tuples=tuples;
	_parallel_radix_clustering_naive((void*)&my_data);
	
}




int main(int argc, char* argv[]) {

	/*Variable declarations....*/
	char filename[MAXLINESIZE];

	uint nthreads, number_tuples, total_bits, t, repeats;
	uint* bits_per_pass;
	uint* tuples;

	strcpy(filename, argv[1]); //We read the filename.

	/*To begin with, we load and validate the data.*/	
	int loadFlag=load_Data(filename, nthreads, repeats, total_passes, bits_per_pass, total_bits, number_tuples, tuples);
	if (loadFlag==0){
		printf("Invalid input.\n");
		return -1;
	}
	else if (loadFlag==-1){
		printf("Missing tuples.\n");
		return -1;
	}
	else if (loadFlag==-2){
		printf("Too many tuples.\n");
		return -1;
	}
	
	uint *result = (uint*) malloc(number_tuples*sizeof(uint)); //This will be the result of the clustering if passes%2!=0.
	uint *result1 = (uint*) malloc(number_tuples*sizeof(uint)); //Auxiliary array, result of the clustering if passes%2==0.

/*	printf("Repeats=%d passes=%d total_bits=%d number_tuples=%d\n", repeats, total_passes, total_bits, number_tuples);
	for (uint i = 0; i < total_passes; i++) {
		printf("Pass Number:%d bits:%d\n", i+1, bits_per_pass[i]);		
	}
	for (uint i=0; i<number_tuples; i++){
		printf("Tuple:%d Value=%d\n", i+1, tuples[i]);
	}

	printf("\n ******************** \n");*/


	/*Execution of the clustering for a given number of repetitions.*/
//	t = timer();
	while(repeats--) {
		parallel_radix_clustering_naive(nthreads, number_tuples, total_passes, total_bits, bits_per_pass, result, result1, tuples); 
	}
//	printf("\n Execution finished. Tuples: %d,  Time(miliseconds): %d\n", number_tuples, (timer()-t));

	/*Printing of results...*/
/*	if(total_passes%2==0){
		for (uint i=0; i<number_tuples; i++){
			printf("RESULTING Tuple pos:%d Value=%d\n", i+1, result1[i]);
		}
	}
	else{
		for (uint i=0; i<number_tuples; i++){
			printf("RESULTING Tuple pos:%d Value=%d\n", i+1, result[i]);
		}

	}
	printf("\n");*/
	return 0;
}
