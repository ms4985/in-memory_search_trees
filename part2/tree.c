#include "tree.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSE4.2
#include <ammintrin.h> //SSE4A
#include <x86intrin.h>
#include <immintrin.h>


extern int posix_memalign(void** memptr, size_t alignment, size_t size);
size_t alignment = 16;

Tree* build_index(size_t num_levels, size_t fanout[], size_t num_keys, int32_t key[]) {
	// return null pointer for invalid tree configuration
	size_t min_num_keys = 1;
	for (size_t i = 0; i < num_levels - 1; ++i) {
		min_num_keys *= fanout[i];
	}
	size_t max_num_keys = min_num_keys * fanout[num_levels - 1] - 1;
	if (num_keys < min_num_keys || num_keys > max_num_keys) {
		fprintf(stderr, "Error: incorrect number of keys, min %zu, max %zu\n", min_num_keys, max_num_keys);
		return NULL;
	}

	// initialize the tree index
	Tree* tree = malloc(sizeof(Tree));
	assert(tree != NULL);
	tree->num_levels = num_levels;
	tree->node_capacity = malloc(sizeof(size_t) * num_levels);
	assert(tree->node_capacity != NULL);
	for (size_t i = 0; i < num_levels; ++i) {
		tree->node_capacity[i] = fanout[i] - 1;
	}
	tree->key_array = malloc(sizeof(int32_t*) * num_levels);
	assert(tree->key_array != NULL);
	size_t* key_count = malloc(sizeof(size_t) * num_levels);
	assert(key_count != NULL);
	size_t* array_capacity = malloc(sizeof(size_t) * num_levels);
	assert(array_capacity != NULL);
	for (size_t i = 0; i < num_levels; ++i) {
		size_t size = sizeof(int32_t) * tree->node_capacity[i];         // allocate one node per level
		int error = posix_memalign((void**) &(tree->key_array[i]), alignment, size);
		assert(error == 0);
		key_count[i] = 0;
		array_capacity[i] = tree->node_capacity[i];     // array_capacity[i] is always a multiple of node_capacity[i]
	}

	// insert sorted keys into index
	for (size_t i = 1; i < num_keys; ++i) {
		assert(key[i - 1] < key[i]);
	}
	for (size_t i = 0; i < num_keys; ++i) {
		size_t level = num_levels - 1;
		while (key_count[level] == array_capacity[level])
			level -= 1;
		tree->key_array[level][key_count[level]] = key[i];
		key_count[level] += 1;
		while (level < num_levels - 1) {
			level += 1;
			size_t new_capacity = array_capacity[level] + tree->node_capacity[level];
			size_t size = sizeof(int32_t) * new_capacity;           // allocate one more node
			int32_t* new_array = NULL;
			int error = posix_memalign((void**) &new_array, alignment, size);
			assert(error == 0);
			memcpy(new_array, tree->key_array[level], sizeof(int32_t) * key_count[level]);
			free(tree->key_array[level]);
			tree->key_array[level] = new_array;
			array_capacity[level] = new_capacity;
		}
	}

	// pad with INT32_MAXs
	for (size_t i = 0; i < num_levels; ++i) {
		for (size_t j = key_count[i]; j < array_capacity[i]; ++j)
			tree->key_array[i][j] = INT32_MAX;
		key_count[i] = array_capacity[i];
	}

	// print the tree
	for (size_t i = 0; i < num_levels; ++i) {
		printf("Level %zu:", i);
		for (size_t j = 0; j < key_count[i]; ++j)
			printf(" %d", tree->key_array[i][j]);
		printf("\n");
	}

	free(array_capacity);
	free(key_count);
	return tree;
}

uint32_t probe_index(Tree* tree, int32_t probe_key) {
	/* previous implementation:*/
	size_t result1 = 0;
	for (size_t level = 0; level < tree->num_levels; ++level) {
		size_t offset = result1 * tree->node_capacity[level];
		size_t low = 0;
		size_t high = tree->node_capacity[level];
		while (low != high) {
			size_t mid = (low + high) / 2;
			if (tree->key_array[level][mid + offset] < probe_key)
				low = mid + 1;
			else
				high = mid;
		}
		size_t k = low;       // should go to child k
		result1 = result1 * (tree->node_capacity[level] + 1) + k;
	}
	printf("Answer should be %zu\n",result1);



	/* ROOT */
/*	printf("probe is %d\n",probe_key);
	printf("root is %d %d %d %d\n",tree->key_array[0][0],tree->key_array[0][1],tree->key_array[0][2],tree->key_array[0][3]);
	size_t low = 0;
	size_t high = tree->node_capacity[0];
	printf("high is %zu\n",high);
	while (low != high) {
		size_t mid = (low + high) / 2;
		printf("key_array[0][mid] is %d when mid is %zu\n",tree->key_array[0][mid], mid);
		if (tree->key_array[0][mid] < probe_key)
			low = mid + 1;
		else
			high = mid;
		printf("high is %zu and low is %zu\n",high, low);
	}
	size_t result = low;

	printf("Root result is %zu\n", low);

	//	int res = 0;
	//	int r_0 = 0; //0 because we start at first item in index array and bitshift over after each load????
	//	int r_1 = 0;
	//	int r_2 = 0;
	int r_3 = 0;

	int rprev = result;
	int r = 0;*/

	int rprev = 0;
	int r = 0;

	for (size_t level = 0; level < tree->num_levels; ++level) {
		printf("rprev is %d at level %d\n", rprev, level);
		if (tree->node_capacity[level] == 4) {
			/* 5-way */
			// access level 1 (non-root) of the index (5-way)
			int32_t *index_L1 = tree->key_array[level];
			__m128i lvl_1 = _mm_load_si128((__m128i*)&index_L1[rprev << 2]);
			
			int *val = (int*) &lvl_1;
			printf("Numerical: %d %d %d %d \n", 
			               val[0], val[1], val[2], val[3]);//, val[4], val[5], 

			__m128i key = _mm_loadl_epi64((__m128i*)&probe_key);
			key = _mm_shuffle_epi32(key, 0);

			
			val = (int*)&key;
			printf("mm load si: %d  \n", val[0]);
			__m128i cmp_1i = _mm_cmpgt_epi32(lvl_1, key);
			val = (int*) &cmp_1i;
			printf("mm cmpgt: %d  \n", val[0]);
			__m128 cmp_1 = _mm_castsi128_ps(cmp_1i);
			r = _mm_movemask_ps(cmp_1); // ps: epi32
			printf("r is %d first\n", r);
			
			
			int t1 = _bit_scan_forward(r); // it seems that we don't need ^ 0x1FF
			r =t1;
			r += (rprev << 2) + rprev;
			rprev = r;
			printf("r is %d\n", r);
		}
		else if (tree->node_capacity[level] == 8) {
			/* 9-way */
			// access level 2 of the index (9-way)
			int32_t *index_L2 = tree->key_array[level] + rprev*tree->node_capacity[level];
			__m128i lvl_2_A = _mm_load_si128((__m128i*)&index_L2[ r << 3]);
			__m128i lvl_2_B = _mm_load_si128((__m128i*)&index_L2[(r << 3) + 4]);
			
			__m128i key = _mm_loadl_epi64((__m128i*)&probe_key);
			key = _mm_shuffle_epi32(key, 0);


			__m128i cmp_2_A = _mm_cmpgt_epi32(lvl_2_A, key);
			int *val = (int*) &cmp_2_A;
			printf("cmp_2_A: %d  \n", val[0]);
			__m128i cmp_2_B = _mm_cmpgt_epi32(lvl_2_B, key);
			val = (int*) &cmp_2_B;
			printf("cmp_2_B: %d  \n", val[0]);
			__m128i cmp_2 = _mm_packs_epi32(cmp_2_A, cmp_2_B);
			val = (int*) &cmp_2;
			printf("cmp_2: %d  \n", val[0]);
			cmp_2 = _mm_packs_epi16(cmp_2, _mm_setzero_si128());
			val = (int*) &cmp_2;
			printf("cmp_2 again: %d  \n", val[0]);
			r = _mm_movemask_epi8(cmp_2);
			int t1 = _bit_scan_forward(r ^ 0x1FFFF);
			int t2 = _bit_scan_forward(r);
			printf("t1 is %d and t2 is %d\n", t1, t2);
			r = _bit_scan_forward(r ^ 0x1FFFF);
			r += (rprev << 3) + rprev;
			printf("r2 is %d\n", r);
			rprev = r;
		}
		else if (tree->node_capacity[level] == 16) {
			/* 17-way */
			// broadcast 1 32-bit key to all SIMD lanes
			//__m128i key = _mm_loadl_epi32(input_keys++); // asm: movd
			//key = _mm_shuffle_epi32(key, 0);
			// from piazza:
			//__m128i key = _mm_cvtsi32_si128(input_keys[i++]);
			__m128i key = _mm_cvtsi32_si128(probe_key);
			key = _mm_shuffle_epi32(key, _MM_SHUFFLE(0,0,0,0));

			//store 16 delimiters in 4 registers
			int32_t *index_level = tree->key_array[level];
			__m128i del_ABCD = _mm_load_si128((__m128i*)&index_level[r_3 << 4]);
			__m128i del_EFGH = _mm_load_si128((__m128i*)&index_level[r_3 << 4]);
			__m128i del_IJKL = _mm_load_si128((__m128i*)&index_level[r_3 << 4]);
			__m128i del_MNOP = _mm_load_si128((__m128i*)&index_level[r_3 << 4]);

			// compare with 16 delimiters stored in 4 registers
			//__m128i tmp = _mm_load_si128( (__m128i*)&probe_key);
			__m128i cmp_ABCD = _mm_cmpgt_epi32(key, del_ABCD);
			__m128i cmp_EFGH = _mm_cmpgt_epi32(key, del_EFGH);
			__m128i cmp_IJKL = _mm_cmpgt_epi32(key, del_IJKL);
			__m128i cmp_MNOP = _mm_cmpgt_epi32(key, del_MNOP);
			// pack results to 16-bytes in a single SIMD register
			__m128i cmp_A_to_H = _mm_packs_epi32(cmp_ABCD, cmp_EFGH);
			__m128i cmp_I_to_P = _mm_packs_epi32(cmp_IJKL, cmp_MNOP);
			__m128i cmp_A_to_P = _mm_packs_epi16(cmp_A_to_H, cmp_I_to_P);
			// extract the mask the least significant bit
			int mask = _mm_movemask_epi8(cmp_A_to_P);
			int res17 = _bit_scan_forward(mask | 0x10000); // asm: bsf
		}
		else {
			printf("Please check node capacity - trying with %zu\n", tree->node_capacity[level]);
			return -1;
		}
		//rprev = rprev * tree->node_capacity[level];
	}
	return (uint32_t) r;//return (uint32_t) result
}

void cleanup_index(Tree* tree) {
	free(tree->node_capacity);
	for (size_t i = 0; i < tree->num_levels; ++i)
		free(tree->key_array[i]);
	free(tree->key_array);
	free(tree);
}


