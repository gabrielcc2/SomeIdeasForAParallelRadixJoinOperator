/*Radix_clustering for skewed data:
This C program implements a serial and recursive version of radix clustering for skewed data. 
It is coded after Boncz et al Radix join implementation,* expanding upon it by allowing input 
(a series of unsigned integers) and skewness (variable cluster sizes). The last of which was
achieved by performing an extra pass for counting the number of tuples per cluster, and 
doing a prefix sum over this result. Because of this our implementation would take about 1.5 
the time of the Boncz et al code.


OBSERVATIONS
- We considered implementing to allow tuples of type char, but conversion took too much time. 
An alternative implementation for chars and floats could be easily done.

INPUT:
To be passed to the program as a file, with the following format:
On the first line: Number of MHz: Used for estimating the number of repetitions of the execution, for performance evaluation purposes.
On the second line: Number of passes
On the third line: Number of bits per pass, separated by spaces or commas
On the forth line: Number of tuples
On the fifth line, onwards: Values of the tuples (unsigned integers).

OUTPUT:
This program prints on console the resulting clustering of the tuples.

FUNCTIONS INCLUDED:
- uint timer(void): Returns the current time in miliseconds.
- int load_Data (char [], uint&, uint&, uint*&, uint&, uint&, uint* &): 
	Loads the data from the file in the given format. Returns 1 if successful, 0 if invalid input, -1 if missing tuples, -2 if too many tuples.
- void prefix_sum(uint*&, uint): Utility function for calculating a prefix sum over an array.
- void radix_cluster_for_skewed_data(uint, uint, uint, uint[], uint*, uint*, uint *): 
	Following Boncz, if the number of passes is even, the result will be stored in the second to last parameter, otherwise in the third to last.
- int main(int, char*[])

OvGU Magdeburg June 2014.

*Available in http://homepages.cwi.nl/~boncz/mimuw/. 

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
typedef unsigned int uint;
#define MAXLINESIZE 128  /* or other suitable maximum line size for reading from file, and maximum tuple length as strings*/
uint timer(void) {
        struct timeval tp;
        gettimeofday(&tp, NULL);
        return tp.tv_sec * 1000 + tp.tv_usec / 1000;
}

/*The following function loads the parameters for the clustering from a file, with the following format:
On the first line: Number of MHz
On the second line: Number of passes
On the third line: Number of bits per pass, separated by spaces or commas
On the forth line: Number of tuples
On the fifth line, onwards: Values of the tuples.
*/
int load_Data (char filename [], uint& MHz, uint& passes, uint*& bits_per_pass, uint& total_bits, uint& number_tuples, uint* &tuples){
	FILE *file = fopen ( filename, "r" );
	if ( file != NULL )
		{
			char line [ MAXLINESIZE ];
			fgets ( line, sizeof line, file );
			if (atoi(line)<=0)
				return 0;
			MHz=atoi(line);
			fgets ( line, sizeof line, file );
			if (atoi(line)<=0)
				return 0;
			passes=atoi(line);
			bits_per_pass=(uint*)malloc(sizeof(uint)*passes);
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
			if (atoi(line)<=0)
				return 0;
			number_tuples=atoi(line);
			tuples=(uint*)malloc(sizeof(uint)*number_tuples);
			uint currentTuple=-1;
			while ( fgets ( line, sizeof line, file ) != NULL ) /* read a line */
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
			perror ( filename ); /* why didn't the file open? */
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
/*Implementation of radix clustering for skewed data (when clusters from a given pass have different number of tuples).*/
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

/*The only parameter of this program will be a filename, with the following elements:
On the first line: Number of MHz
On the second line: Number of passes
On the third line: Number of bits per pass, separated by spaces or commas
On the forth line: Number of tuples
On the fifth line, onwards: Values of the tuples.
*/
int main(int argc, char* argv[]) {

/*Variable declarations....*/
	char filename[MAXLINESIZE];

	uint MHz=0; //Used to decide number of repeats for testing.

	uint passes;
	uint* bits_per_pass;
	uint total_bits;
	uint number_tuples;
	uint* tuples;

	uint t;

	strcpy(filename, argv[1]); //We read the filename.



/*To begin with, we load and validate the data.*/	
	int loadFlag=load_Data(filename, MHz, passes, bits_per_pass, total_bits, number_tuples, tuples);
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
	
/*Initialization of additional data and printing of input*/
	uint repeats = (100000LL * MHz) / number_tuples; /* repeats. idea: repeats*number_tuples = Hz of the machine.. Repetitions are for performance evaluation.*/
	uint *result = (uint*) malloc(number_tuples*sizeof(uint)); //This will be the result of the clustering if passes%2!=0.
	uint *result1 = (uint*) malloc(number_tuples*sizeof(uint)); //Auxiliary array, result of the clustering if passes%2==0.

	printf("MHz=%d passes=%d total_bits=%d number_tuples=%d repeats=%d\n", MHz, passes, total_bits, number_tuples, repeats);
	for (uint i = 0; i < passes; i++) {
		printf("Pass Number:%d bits:%d\n", i+1, bits_per_pass[i]);		
	}
	for (uint i=0; i<number_tuples; i++){
		printf("Tuple:%d Value=%d\n", i+1, tuples[i]);
	}

	printf("\n ******************** \n");


/*Execution of the clustering for a given number of repetitions.*/
	t = timer();
	while(repeats--) {
		radix_cluster_for_skewed_data(number_tuples, passes, total_bits, bits_per_pass, result, result1, tuples); 
	}
	printf("\n Execution finished. Tuples: %d,  Time(miliseconds): %d\n", number_tuples, (timer()-t));

/*Printing of results...*/
	for (uint i=0; i<number_tuples; i++){
		printf("RESULTING Tuple pos:%d Value=%d\n", i+1, (passes%2==0)?result1[i]:result[i]);
	}
	printf("\n");
	return 0;
}
