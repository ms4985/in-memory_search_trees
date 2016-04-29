#define _BSD_SOURCE 

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include "p2random.h"
#include "tree.h"

int main(int argc, char* argv[]) {
        // parsing arguments
        assert(argc > 3);
        size_t num_keys = strtoull(argv[1], NULL, 0);
        size_t num_probes = strtoull(argv[2], NULL, 0);
        size_t num_levels = (size_t) argc - 3;
        size_t* fanout = malloc(sizeof(size_t) * num_levels);
        assert(fanout != NULL);
        for (size_t i = 0; i < num_levels; ++i) {
                fanout[i] = strtoull(argv[i + 3], NULL, 0);
                assert(fanout[i] >= 2 && fanout[i] <= 17);
        }

	struct timeval p1_begin, p1_end, p2_1_begin, p2_1_end, p2_2_begin, p2_2_end, p3_begin, p3_end, 
			p1res, p2_1res, p2_2res, p3res;
        
	gettimeofday(&p1_begin, NULL);

        // building the tree index
        rand32_t* gen = rand32_init((uint32_t) time(NULL));
        assert(gen != NULL);
        int32_t* delimiter = generate_sorted_unique(num_keys, gen);
        assert(delimiter != NULL);
        Tree* tree = build_index(num_levels, fanout, num_keys, delimiter);
        free(delimiter);
        free(fanout);
        if (tree == NULL) {
                free(gen);
                exit(EXIT_FAILURE);
        }

        // generate probes
        int32_t* probe = generate(num_probes, gen);
        assert(probe != NULL);
        free(gen);
        uint32_t* result = malloc(sizeof(uint32_t) * num_probes);
        assert(result != NULL);

        gettimeofday(&p1_end, NULL);
        gettimeofday(&p2-1_begin, NULL);

        // perform index probing (Phase 2 Part 1) 
        for (size_t i = 0; i < num_probes; ++i) {
                result[i] = probe_index_1(tree, probe[i]);
        }

        gettimeofday(&p2-1_end, NULL);
        gettimeofday(&p2-2_begin, NULL);

        // perform index probing (Phase 2 Part 2) 
        for (size_t i = 0; i < num_probes; ++i) {
                result[i] = probe_index_2(tree, probe[i]);
        }

        gettimeofday(&p2-2_end, NULL);
        gettimeofday(&p3_begin, NULL);

        // output results
        for (size_t i = 0; i < num_probes; ++i) {
                fprintf(stdout, "%d %u\n", probe[i], result[i]);
        }

        gettimeofday(&p3_end, NULL);

        timersub(&p1_end, &p1_begin, &p1res);
        timersub(&p2-1_end, &p2-1_begin, &p2-1res);
        timersub(&p2-2_end, &p2-2_begin, &p2-2res);
        timersub(&p3_end, &p3_begin, &p3res);


        printf("Time elapsed for Phase1: %ld.%06ld s\n", (long int)p1res.tv_sec, (long int)p1res.tv_usec);
        printf("Time elapsed for Phase2 Part1: %ld.%06ld s\n", (long int)p2-1res.tv_sec, (long int)p2-1res.tv_usec);
        printf("Time elapsed for Phase2 Part2: %ld.%06ld s\n", (long int)p2-2res.tv_sec, (long int)p2-2res.tv_usec);
        printf("Time elapsed for Phase3: %ld.%06ld s\n", (long int)p3res.tv_sec, (long int)p3res.tv_usec);

        // cleanup and exit
        free(result);
        free(probe);
        cleanup_index(tree);
        return EXIT_SUCCESS;
}
