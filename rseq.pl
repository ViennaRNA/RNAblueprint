#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <2000-03-10 16:33:03 xtof>

use Getopt::Long;
use strict;
use vars qw/$opt_a $opt_l $opt_n @A/;

usage() unless GetOptions("a=s" => \$opt_a,
                           "l=i" => \$opt_l,
			   "n=i" => \$opt_n);
srand;
@A = (defined $opt_a) ? (split //, $opt_a) : qw/A U G C/;
$opt_l = ($opt_l > 1) ? $opt_l : 1;
$opt_n = ($opt_n > 1) ? $opt_n : 1; 

for(1..$opt_n) {
  print join '', @A[map {rand @A} 1..$opt_l], "\n";
}

sub usage {
  printf STDERR "\nusage: rseq [-l length] [-n number] [-a alphabet]\n";
  printf STDERR "default: length = 1 number = 1 alphabet = AUGC\n\n";  
  exit;
}

# End of file
