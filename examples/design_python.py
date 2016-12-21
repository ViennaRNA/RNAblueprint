#!/usr/bin/python
# This script is an example implementation on how to use the Python
# interface. It generates 1000 neighbors of an initially sampled
# random sequence.

import RNAblueprint as rd

# define structures
structures = ['(((((....)))))', '(((....)))....']
# construct dependency graph with these structures
dg = rd.DependencyGraphMT(structures)

# print this sequence
print dg.get_sequence()

# mutate globally for 1000 times and print
for i in range(0, 1000):
    dg.sample_clocal()
    print dg.get_sequence()
    # revert to the previous sequence
    dg.revert_sequence();

# print the amount of solutions
print 'Maximal number of solutions: ' + str(dg.number_of_sequences())
# print the amount of connected components
print 'Number of Connected Components: ' + str(dg.number_of_connected_components())

# make a deep copy of the dependency graph
dg1 = rd.DependencyGraphMT(dg)
