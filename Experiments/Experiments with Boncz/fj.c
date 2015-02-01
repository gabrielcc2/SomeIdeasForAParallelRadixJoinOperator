#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#define PRIME 1442968193

typedef unsigned int uint;
typedef struct { uint key, pos; } tuple;

/* after VLDB'04 Manegold,Boncz,Nes,Kersten */
void radix_decluster(uint b, uint lowbits, uint window, uint *__restrict result, tuple *__restrict input) {
	uint i, n = 1 << b, clustersize = 1 << lowbits, nborders = 2 << (b-lowbits);
	uint inserted = 0, limit = window;

	/* border structure, this time a hard start- and end-pointer */
	tuple **borders = (tuple**) alloca(nborders*sizeof(tuple*));
	borders[1] = (borders[0] = input) + clustersize;
	for(i=2; i<nborders; i+=2) borders[i+1] = (borders[i] = borders[i-1]) + clustersize;

	/* the radix-decluster algorithm */
	while(1)
	for(i=0; i<nborders; i+=2) {
		while (borders[i] < borders[i+1] && borders[i]->pos < limit) {
			result[borders[i]->pos] = borders[i]->key;
			if (++inserted == window) {
				if ((n -= window) == 0) return;
				limit += window;
				inserted = 0;
			} 
			borders[i]++;
		}
	}
}

/* after VLDB'99 Boncz,Manegold,Kersten */
//radix_cluster(number of tuples, p number of passes, b=log2 of the number of tuples or the number of bits of the numbers, nbits: the array with bits per pass, (p&1)?result0:result1, (p&1)?result1:result0, input); 
void radix_cluster(uint n, uint npasses, uint totbits, uint bits[], tuple*__restrict result0, tuple *__restrict result1, tuple *__restrict input);
void radix_cluster(uint n, uint npasses, uint totbits, uint bits[], tuple*__restrict result0, tuple *__restrict result1, tuple *__restrict input) {
	uint i, clustersize = n >> bits[0], nclusters = 1 << bits[0], mask = nclusters-1, lowbits = totbits-bits[0];
	int s;

	/* as the data is a permutation of a power of two, all groups are of equal size (KISS) */
	uint *borders = (uint*) alloca(nclusters*sizeof(uint)); //Dont pay attention to this, this is only algebra and it is the head position of each cluster.
	for(borders[0]=0, i=1; i<nclusters; i++) 
		borders[i] = borders[i-1]+clustersize;

	/* do the partitioning work of this pass */
	for(i=0; i < n; i++) {
		uint clustnr = (input[i].key >> lowbits) & mask;
		result0[borders[clustnr]++] = input[i]; /* most pain is here */
	}
	if (npasses > 1) /* more passes to do? call recursively */
		for(s = n-clustersize; s >= 0; s -= clustersize)
			radix_cluster(clustersize, npasses-1, lowbits, bits+1, result1+s, result0+s, result0+s);
}

/* simple prefetch */
void fetchjoin(uint n, uint*__restrict result, tuple*__restrict input, uint*__restrict fetch) {
	uint i; 
	for(i=0; i<n; i++) result[i] = fetch[input[i].key];
}

/* same, but also propagate the 'pos' in the result tuple (decluster needs pos), and prefetch */
uint fetchjoin_clustered(uint n, uint clustersize, tuple*__restrict result, tuple*__restrict input, uint*__restrict fetch) {
	uint dummy=0, j=0, i=0; 
	while(i < n) {
		/* prefetch input[] sequentially */
		uint limit = j + clustersize;
		if (limit > n) {
			limit = n;
		} else if (input[i].key < limit) { 
			for(; j < limit; j += 64)
				dummy += fetch[j] + fetch[j+16] + fetch[j+32] + fetch[j+48];
		} else {
			j = limit; /* no values in this cluster */
		}
	
		/* do the fetch work */
		for(; i < j; i++) {
			result[i].key = fetch[input[i].key];
			result[i].pos = input[i].pos;
		}
	}
	return dummy;
}

/* some helper functions */
uint hash(uint x) {
	unsigned long long v = ((unsigned long long) x)*PRIME;
	return (v >> 32) ^ v;  
}
tuple *permute(uint n) { /* create a permutation of [0..n> (actually, a random cycle) */
	uint *tmp = (uint*) malloc(n*sizeof(uint));
	tuple *res = (tuple*) malloc(n*sizeof(tuple));
	uint i,chain=n-1; /* create a chain that starts at the last pos */
	for(i=0; i<n; i++) tmp[i] = i;
	n--;
	for(i=0; i<n; i++) {
		uint j = hash(i) % (n-i); /* select a random location from what is left in tmp */
		uint pos = tmp[j]; /* get an unvisited random position */
		res[chain].key = pos; 
		res[chain].pos = chain; 
		chain = pos;
		tmp[j] = tmp[n-(i+1)];    /* move the last of tmp into the selected location */
	}
	res[chain].key = n;
	res[chain].pos = chain;
	free(tmp);
	return res;
}
void print(uint v, uint* fetch, tuple* input, tuple* result0, tuple* result1, uint* fetched) {
	uint i;
	printf("          +--------------------------------+---------------------+\n");//---------------------+----------+
	printf("          +   fetch* |               input | after radix-cluster |\n");//    after fetchjoin2 |          |\n");
	printf("          +----------+----------+----------+----------+----------+\n");//----------+----------+----------+\n");
	printf("          |   fetch* |      key |      pos |      key |      pos |\n");//      key |      pos | fetched* |\n");
	printf("          +----------+----------+----------+----------+----------+\n");//----------+----------+----------+\n");
	for(i=0; i<v; i++)
		printf(" % 8u | % 8u | % 8u | % 8u | % 8u | % 8u |\n", i, 
			fetch[i], input[i].key, input[i].pos, result0[i].key, result0[i].pos);
	printf("          +----------+----------+----------+----------+----------+\n");//----------+----------+----------+\n");
}

uint timer(void) {
        struct timeval tp;
        gettimeofday(&tp, NULL);
        return tp.tv_sec * 100 + tp.tv_usec / 10000;
}

/* usage: prog <datasize_bits> <MHz> <windowsize_bits> <verbose_ntuples> <ncols> (<pass_bits>)* */
int main(uint argc, char* argv[]) { //NOTE: This algorithm does a most significant bit first clustering.
	uint b = atoi(argv[1]); /* total significant bits (ie datasize). This implies 2 to the b tuples. */
	uint n = 1 << b; /* number of tuples, given as log2... This is equivalent to raising 2 to the power of b. */
	uint r = (100000LL * atoi(argv[2])) / n; /* repeats. idea: r*n = Hz of the machine.. The repetitions are only done for performance evaluation.*/
	uint w =0; // 1 << atoi(argv[3]); /* decluster window (0=run cluster-only) */

	uint v = atoi(argv[4]); /* verbose output of up to v tuples */
	uint c = atoi(argv[5]); /* verbose output of up to v tuples */

	uint *nbits = (uint*) alloca(argc*sizeof(uint)); /* an array, because each pass can have a different number of bits */

	/*UNUSED*/
	uint *fetch = (uint*) malloc(c*n*sizeof(uint));
	uint *fetched = (uint*) malloc(c*n*sizeof(uint));
	/*UNUSED*/

	tuple *result0 = (tuple*) malloc(n*sizeof(tuple)); //This will be the result of the clustering.
	tuple *result1 = (tuple*) malloc(n*sizeof(tuple)); //This wont be used.
	
	tuple *input = permute(n); //This will be the attribute to cluster. 
	uint p, i, t, lowbits = b; /* low bits are the lowermist bits that are left unordered by radix-cluster */
	for(p=0; p<argc-6; p++) {
		nbits[p] = atoi(argv[6+p]); lowbits -= nbits[p];
		printf("radix-cluster-pass%d=%d ", p+1, nbits[p]);
	} //Filling the array with the number of bits per pass.

	/*UNUSED*/
	for(i=0; i<n; i++) { 
		for(t=0; t<c; t++) fetch[i+t*n] = i&255;
	}
	/*UNUSED*/

	printf("datasize=%d fetchcols=%d decluster-window=%d repeats=%d\n", n, c, w, r);
	t = timer();
	while(r--) { /* do the experiment */

	/*UNUSED*/
		if (w == 1) {
			for(i=0; i<c; i++) {
				fetchjoin(n, fetched+i*n, input, fetch+i*n);
			}

		/*UNUSED*/
		} else {
			radix_cluster(n, p, b, nbits, (p&1)?result0:result1, (p&1)?result1:result0, input); 

			
			//for(i=0; i<c; i++) {
			//	fetchjoin_clustered(n, 1<<lowbits, result1, result0, fetch+i*n);
			//	radix_decluster(b, lowbits, w, fetched+i*n, result1);
			//}
		}
		//input[r & (n-1)].key &= (1<<30)-1;  /* confuse the compiler that we modify input (but we don't) */
	}
	if (v) print((v > n)?n:n, fetch, input, result0, result1, fetched);
	printf("%d %d cycles\n", n, (timer()-t)/10);
	return fetched[0] + result0[0].pos;  /* dummy value to force intelligent compilers to do the work */
}
