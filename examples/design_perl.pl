#!/usr/bin/perl
# This script is an example implementation on how to use the Perl
# interface. It generates 1000 neighbors of an initially sampled
# random sequence.

use RNAblueprint;

# define structures
@structures = ['(((((....)))))', '(((....)))....'];
# construct dependency graph with these structures
$dg = new RNAblueprint::DependencyGraphMT(@structures);

# print this sequence
print $dg->get_sequence()."\n";

# mutate globally for 1000 times and print
for($i=0; $i<1000; $i++) {
    $dg->sample_clocal();
    print $dg->get_sequence()."\n";
    # revert to the previous sequence
    $dg->revert_sequence();
}
# print the amount of solutions
print 'Maximal number of solutions: '.$dg->number_of_sequences()."\n";
# print the amount of connected components
print 'Number of Connected Components: '.$dg->number_of_connected_components()."\n";

# make a deep copy of the dependency graph
$dg1 = new RNAblueprint::DependencyGraphMT($dg)
