#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

	clock_t p1_begin, p1_end, p2_begin, p2_end, p3_begin, p3_end;

        p1_end = p2_begin = clock();

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

        p1_end = p2_begin = clock();

        // perform index probing (Phase 2)
        for (size_t i = 0; i < num_probes; ++i) {
                result[i] = probe_index(tree, probe[i]);
        }

        p2_end = p3_begin = clock();

        // output results
        for (size_t i = 0; i < num_probes; ++i) {
                fprintf(stdout, "%d %u\n", probe[i], result[i]);
        }

        p3_end = clock();

        printf("Phase 1: %f s\n",(double)(p1_end - p1_begin) / CLOCKS_PER_SEC);

        printf("Phase 2: %f s\n",(double)(p2_end - p2_begin) / CLOCKS_PER_SEC);

        printf("Phase 3: %f s\n",(double)(p3_end - p3_begin) / CLOCKS_PER_SEC);

        // cleanup and exit
        free(result);
        free(probe);
        cleanup_index(tree);
        return EXIT_SUCCESS;
}
