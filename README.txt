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



- Example runs:
./build 50 3 9 5 9
Level0
1449016084	


Level1
-1051168240	-472517129	237491495	541601439	


Level2
-1810036580	-1599780957	-1573766484	-1514568422	-1505812806	-1331822453	-1229313089	-1155030660	-104382466-981599069	-958206951	-812417474	-651452045	-642828372	-593907325	-552767594	-471374245	-316377796-260089207	-176981921	-97677375	-26528976	-9175992	134907146	258088631	354936737	391954834	393075426	402973709	435894283	493642829	494866564	622759867	633585167	698693805	820590273	956615958	1154493465	1265501589	1291610599	1461479270	1751855557	1845312575	1980224785	2081700233


-355213889	17
630285875	35
-263088311	17
Phase 1: 0.000066 s
Phase 2: 0.000113 s
Phase 3: 0.000014 s




./build 11 3 3 2 3
Level0
232602094	


Level1
-1251147284	1100846421	


Level2
-1707383052	-1445936979	-958639932	-707725292	746596123	1045385710	1613769224	1945037450	


-643477281	3
-1108076324	2
-1579730557	1
Phase 1: 0.000041 s
Phase 2: 0.000086 s
Phase 3: 0.000012 s




./build 40 3 17 2 17
Level0
1290957356	


Level1
-546706059	


Level2
-2057350520	-1974399504	-1934591011	-1842337521	-1790143040	-1548634749	-1402943542	-1377607196	-118885380-1100023865	-1064717560	-1045190764	-1002942309	-942006334	-938742520	-656609891	-519880049	-457877766-305869827	-226168246	-7740387	17493515	21549018	214290401	243565388	536747416	673532960	1054254703	1109365972	1135111768	1182606707	1236053646	1422830359	1596826650	1622912848	19359817852022349359	2117333771	


-44866442	19
-1509346172	5
1325968993	32
Phase 1: 0.000055 s
Phase 2: 0.000106 s
Phase 3: 0.000012 s
