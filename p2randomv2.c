#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <time.h>

typedef struct {
	size_t index;
	uint32_t num[625];
} rand32_t;


rand32_t *rand32_init(uint32_t x)
{
	rand32_t *s = malloc(sizeof(rand32_t));
	uint32_t *n = s->num;
	size_t i = 1;
	n[0] = x;
	do {
		x = 0x6c078965 * (x ^ (x >> 30));
		n[i] = x;
	} while (++i != 624);
	s->index = i;
	return s;
}

uint32_t rand32_next(rand32_t *s)
{
	uint32_t x, *n = s->num;
	size_t i = s->index;
	if (i == 624) {
		i = 0;
		do {
			x = (n[i] & 0x80000000) + (n[i + 1] & 0x7fffffff);
			n[i] = (n[i + 397] ^ (x >> 1)) ^ (0x9908b0df & -(x & 1));
		} while (++i != 227);
		n[624] = n[0];
		do {
			x = (n[i] & 0x80000000) + (n[i + 1] & 0x7fffffff);
			n[i] = (n[i - 227] ^ (x >> 1)) ^ (0x9908b0df & -(x & 1));
		} while (++i != 624);
		i = 0;
	}
	x = n[i];
	x ^= (x >> 11);
	x ^= (x <<  7) & 0x9d2c5680;
	x ^= (x << 15) & 0xefc60000;
	x ^= (x >> 18);
	s->index = i + 1;
	return x;
}

int int32_cmp(const void *x, const void *y)
{
	int32_t a = * (const int*) x;
	int32_t b = * (const int*) y;
	return a < b ? -1 : a > b ? 1 : 0;
}

int32_t *generate(size_t n, rand32_t *gen)
{
	size_t i;
	int32_t *a = malloc(n << 2);
	for (i = 0 ; i != n ; ++i)
		a[i] = rand32_next(gen);
	return a;
}

int32_t *generate_sorted_unique(size_t n, rand32_t *gen)
{
	size_t i = 0;
	size_t m = n / 0.7;
	uint8_t z = 0;
	uint32_t *a = malloc(n << 2);
	uint32_t *b = calloc(m, 4);
	while (i != n) {
		uint32_t k = rand32_next(gen);
		if (k != 0) {
			size_t h = (uint32_t) (k * 0x9e3779b1);
			h = (h * (uint64_t) m) >> 32;
			while (b[h] != k) {
				if (b[h] == 0) {
					b[h] = a[i++] = k;
					break;
				}
				if (++h == m) h = 0;
			}
		} else if (z == 0) {
			a[i++] = 0;
			z = 1;
		}
	}
	free(b);
	qsort(a, n, 4, int32_cmp);
	return (int32_t*) a;
}

void ratio_per_bit(const int32_t *a, size_t n)
{
	size_t i, j, *c = calloc(32, sizeof(size_t));
	for (i = 0 ; i != n ; ++i) {
		int32_t x = a[i];
		for (j = 0 ; j != 32 ; ++j)
			c[j] += (a[i] >> j) & 1;
	}
	for (j = 0 ; j != 32 ; ++j)
		;//fprintf(stderr, "%2ld: %.2f%%\n", j + 1, c[j] * 100.0 / n);
	free(c);
}

void create_tree(int32_t *tree[], int *index, int32_t *keys, int n, char **level, int numLevels) {
	// TODO: error checking that level[i] is a long (strtol function)
	// TODO: free extra space in array
	// TODO: use memset or something more efficient to initialize index array
	// TODO: free memory
	// TODO: check for off by one errors, especially w/ fanout vs. capacity
	// TODO: make sure strtol returns > 2

	int error = 0;
	
	// create array which holds pointers to the start of each level
	// tree[0] = root level; tree[1] = next ... tree[numLevels - 1] = leaves
	//int32_t *tree[size];// = (int32_t *)malloc(sizeof(int32_t *) * numLevels);

	int i = 0, j = 0;
	long fanout = 1;
	long arraySize = 1;
	char *ptr;

	// index used to keep track of how to fill up structure
	//int index[numLevels];
	// sizes used to keep track of # keys per node per level (level 0 = root)
	int sizes[numLevels];
	// which level are you inserting on
	int tlevel = numLevels-1;
	// capacity of each level
	int maxSizes[numLevels];

	// tells the fanout at each level  - can use strtol(level[i]...) instead
	//int fans[numLevels];

	//printf("%d",numLevels);
	// TODO: fix this?
	for (i = 0; i<numLevels; i++) {
		index[i] = 0;
	}

	// allocate space for each level, store pointers in tree[]   
	for (i = 0; i < numLevels; i++) {
	//	printf("");
		arraySize = (strtol(level[i], &ptr, 10) - 1) * fanout;
		fanout = fanout * (strtol(level[i], &ptr, 10));

		sizes[i] = strtol(level[i], &ptr, 10) - 1;
		maxSizes[i] = arraySize;
	//	printf("%d\n", maxSizes[i]);
	//	printf("%ld",strtol(level[i],&ptr,10));
	//	printf("Array size at level %d is %ld\n", i, arraySize);
		int32_t *l = malloc(sizeof(int32_t) * arraySize);
		tree[i] = l;
	//	printf("Address of level %d is %p\n", i, tree[i]);
	//	printf("Size at level %d is %d\n", i, sizes[i]);

	}

	// Now need to fill in array!
	// n = number of keys
	for (i = 0; i < n; i++) {
		// figure out which level 
	//	printf("Index is %d at level %d\n", index[tlevel], tlevel);
	//	printf("%p\n",tree[tlevel]+i);
	//	if (index[tlevel] > maxSizes[tlevel]) {
			// have too many keys
	//		error = 1;
	//		break;
	//	}

		// index[tlevel] is the offset for that array
		*(tree[tlevel]+index[tlevel]) = *(keys + i);
	//	printf("loading %d into the tree\n", *(tree[tlevel]+index[tlevel]));
		index[tlevel]++;

		// TODO: how to know when to change levels?

		// if not at leaf level, drop back down
		if (tlevel != (numLevels-1)) {
			tlevel = numLevels - 1;
		}
		//todo  greater than length
		// deal with case if sizes[tlevel] is 1
		else {
			if (sizes[tlevel] == 1) {
				tlevel--;
			}
			else {
				int k = sizes[tlevel];
				int tl = tlevel;
				while (index[tlevel] % k == 0) {// && index[tlevel] % (sizes[tlevel]*2) != 0 && index[tlevel] != 0) {
					tl--;
					if (tlevel < 0) {
						tl = numLevels - 1;
						break;
					}
					k = k * (strtol(level[tl], &ptr, 10));
					if (tlevel < 0) {
						tlevel = numLevels - 1;
						break;
					}
				}
				tlevel = tl;
				
			}
		}
	}
		
	//check if leaf nodes are partially filled, pad with MAXINT if so
/*	int m = numLevels-1;
	while(m+1>0){
		int mx = sizes[m];
		int idx = index[m] % mx;
		int k = m;
		if(idx != 0) {
			printf("need to fill with maxInt\n");
			while(idx != mx) {
				*(tree[m] + index[k]) = INT_MAX;
				idx++;
				k++;
			}
		}
		m--;
	}	*/
//	printf("%d\n", INT_MAX);
		//*(tree[tlevel]+index[tlevel]) = *(keys + i);
	//(index[tlevel] > maxSizes[tlevel]) 
	//pad each array with MAX INT
	int lvl = numLevels - 1;
/*	int idx;
	lvl = 0;
	while(lvl < numLevels) {
		idx = index[lvl]; 
		while(idx < maxSizes[lvl]+1) {
			*(tree[lvl] + idx) = INT_MAX;
			idx++;
		}
		lvl++;
	}
	*/
	int idx;

	while(lvl >= 0) {
		// TODO check for off by 1
	//	printf("level is %d\n",lvl);
		idx = index[lvl];
		while (idx < maxSizes[lvl]) {
					//printf("adding 2 to tree, index is %d and maxSizes is %d\n",index[lvl], maxSizes[lvl]);
			*(tree[lvl]+idx) = INT_MAX;
	//		printf("%d\n", *(tree[lvl]+index[lvl]));
			idx++;
		}
		//index[lvl]--;
		lvl--;
	}


}
//return 0 for exact key or 0 for move left on node, 1 on move right
int check_node(int key, int val) {

	printf("\n%d vs\n%d\n", key, val);
	
	if(key <= val)
		return 0;
	else 
		return 1;
}

int binary_search(int32_t *tree [], int fan [], int probe, int max){

	//right now func will output the index of tree, not the range, need a way to calc range

	int idx,  out, right, left;
	int lvl = 0; //keeps track of which level of tree you are on
	int node = 0;//keeps track of which node you are on
	int range = 0;

	//iterate through each level of tree, searching for range
	while(lvl < max){
		idx = (int)(((fan[lvl]-1)/2)-1); //initial pt of search, fanout tells # of slots per node
		if(idx<0) {
			idx = 0;
		}
		right = 0;
		left = 0;
		while((idx >= 0) && (idx < (fan[lvl])-1)){ //while on same node
			printf("\n%d   %d   %d\n", lvl, node, idx);
			out = check_node(probe, *(tree[lvl]+idx+node));//compare node
			//printf("out:%d\n", out);
			if(out == 0 && right == 0){ //move left but never moved right- go left
				left = 1;
				idx--;	
			}
			else if (out == 0 && right == 1){ //move left but moved right already- go down tree from cur idx
				if(lvl == max-1) {
					//set range
					//printf("lvl: %d node: %d idx %d\n", lvl, node, idx);
					range = idx+node;
					return range;
				}
				node = idx*fan[lvl+1];
				break;
			}	
			else if (out == 1 && left == 0){ //move right, but never moved left- go right
				right = 1;
				idx++;
			}
			else if(out == 1 && left == 1){ //move right, but moved left already- go down tree from cur idx
				if(lvl == max-1) {
					//set range
					//printf("lvl: %d node: %d idx %d\n", lvl, node, idx);
					range = idx+node;
					return range;
				}
				node = idx*fan[lvl+1];
				break;
			}	
		}
		//printf("idx after break: %d\n", idx);
		if (idx < 0) { //move down left
			//printf("lvl %d\n", lvl);
			if(lvl == max-1) {
				//set range
				//printf("lvl: %d node: %d idx %d\n", lvl, node, idx);
				range = idx+node+1;
				return range;
			}
			//printf("move down left\n");
			node = node*fan[lvl+1];
		}
		else if (idx >= fan[lvl]-1) { //move down right
			//printf("lvl %d\n", lvl);
			if(lvl == max-1) {
				//set range
				//printf("lvl: %d node: %d idx %d\n", lvl, node, idx);
				range = idx+node-1;
				return range;
			}
			//printf("move down right\n");
			node = (idx*(fan[lvl+1]))-1;
		}	
		
		lvl++;
	}

	return -1;
}

int main(int argc, char **argv)
{
	clock_t p1_begin, p1_end, p2_begin, p2_end, p3_begin, p3_end;

// begin phase 1
	p1_begin = clock();

	// arg0: a.out
	// arg1: K
	// arg2: P
	// arg3-x: depth at root-leaves
	rand32_t *gen = rand32_init(time(NULL));
	size_t i, n = argc > 1 ? atoll(argv[1]) : 10;
	size_t n2 = argc > 1 ? atoll(argv[2]) : 10;
	// keys
	int32_t *a = generate_sorted_unique(n, gen);
	// probes
	int32_t *p = generate(n2, gen);
	free(gen);
	for (i = 1 ; i < n ; ++i)
		assert(a[i - 1] < a[i]);
	ratio_per_bit(a, n);

	int level = argc-3;
	int32_t *tree[level];
	int *index = malloc(level*sizeof(int));
	int fan[level];
	int maxSize[level];
	int maxKeys = 0;
	int k, j;
	int f = 1;
	/*
	printf("keys\n");
	for (i = 0 ; i < n ; ++i) {
		printf("%d\n", a[i]);
	}
	printf("probes\n");
	for (i=0; i< n2; ++i) {
		printf("%d\n", p[i]);
	}
	*/

	// TODO: make sure arc2 is > 1
	printf("n: %zu, argc-3: %d\n", n, level);

	//fill array of fanouts
	for(k = 0; k < level; k++) {
		fan[k] = atoi(argv[3+k]);
	}
	
	//fill array of max number of nodes at each level
	for(k = 0; k < level; k++) {
		maxSize[k] = f * (fan[k]-1);
		//printf("%d\n", maxSize[k]);
		f = f * fan[k];
		maxKeys = maxKeys + maxSize[k];
	}

	if(n > maxKeys) {
		fprintf(stderr, "ERROR: Too many keys!\n\n");
		return EXIT_FAILURE;
	}

	//fill tree
	create_tree(tree, index, a, n, &argv[3], level);

// end phase 1, begin phase 2
	p1_end = p2_begin = clock();

	//print index
	for(k = 0; k < level; k++) {
		printf("%d\n",index[k]);
	}
	
	if(index[0] == 0) {
		fprintf(stderr, "ERROR: Not enough keys! Root is empty!\n\n");
		return EXIT_FAILURE;
	}
	
	//print tree
	for (k = 0; k < level; k++) {
		printf("Level%d\n", k);
		for( j = 0; j < index[k]; j++) {
			if(*(tree[k]+j)){
				printf("%d\t", *(tree[k]+j));
			}
		}
		printf("\n\n\n");
	}

// end phase 2, begin phase 3
	p2_end = p3_begin = clock();

	int *out = malloc(n2*sizeof(int));

	printf("probes\n");
	for(k = 0; k < n2; k++) {
		printf("%d\n", p[k]);
		out[k] = binary_search(tree, fan, p[k], level);
	}

	printf("out\n");
	for(k = 0; k < n2; k++) {
		printf("%d\n", out[k]);
	}

// end phase 3
	p3_end = clock();

	printf("Phase 1: %f s\n",(double)(p1_end - p1_begin) / CLOCKS_PER_SEC);

	printf("Phase 2: %f s\n",(double)(p2_end - p2_begin) / CLOCKS_PER_SEC);

	printf("Phase 3: %f s\n",(double)(p3_end - p3_begin) / CLOCKS_PER_SEC);


	free(a);
	free(p);
	free(index);
	free(out);
	return EXIT_SUCCESS;
}
