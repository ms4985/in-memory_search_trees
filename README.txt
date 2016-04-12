Carolyn Fine crf2133
Megan Skrypek ms4985

Project 2 - Traversal of In-Memory Search Trees

Part 1:

- Instructions to run: 	make
						./build K P N0 ... Ni
	- K = number of keys to generate
	- P = number of probes to generate
	- N0 ... Ni fanout at each level (N0 = root; Ni = leaves)

- We used the p2randomv2.c code as a skeleton for our search tree.  We generated the keys using generate_sorted_unique() and the probes using generate().  We then check that there are not more keys than we have capacity for.  Next we call create_tree(), which allocates memory and fills in the array for each level.  Though there are not pointers between levels (i.e. from parent to child), we do store a pointer to the head of each array (level).  The completed tree is a valid key-only B+ tree.  When that function returns, we check that the root is not empty.  Then, for each probe we run binary_tree() and store the output in the out[] variable.  As a last step, we print out the out[] array.

- We used the clock() function and clock_t variables to measure the time in each phase.  Six variables were initialized, corresponding to the start and end of each phase.  Before main() exits we printed out each phase's execution time.
	NOTE: times are always 0, unless run with valgrind because the methods do not take that much time, valgrind slows them down enough to see the differences.
