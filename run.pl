#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <2012-12-20 12:27:53 jango>

use Getopt::Long;
use Pod::Usage;
use warnings;
use FileHandle;
use IPC::Open2;

# global variables
my $bin_subopt = "RNAsubopt";
my $bin_barriers = "../Barriers-install/bin/barriers";
my $bin_treekin = "../treekin-install/bin/treekin";

# option parser
my $man = 0;
my $help = 0;
my $infile;
my $outfile = "out.txt";
my $verbose = 0;

pod2usage(-verbose => 0)
    unless GetOptions (	"in=s" => \$infile,
    			"out=s" => \$outfile,
			"verbose" => \$verbose,
			"man"   => sub{pod2usage(-verbose => 2)},
			"help"  => sub{pod2usage(-verbose => 1)});
	
print "Input file is: $infile\n" if($verbose);


my $cmd_subopt = "$bin_subopt -e 10 -s";
my $cmd_barriers = "$bin_barriers -G RNA -noLP --bsize --max 1000 --rates";
my $cmd_treekin = "$bin_treekin -m I --p0 4=1 --ratesfile=rates.out";

my @input;
my @out_subopt;
my @out_barriers;


print STDERR "$command\n" if($verbose);


open(INPUT, "$infile") || die "can't open input file: $!";
 while (<INPUT>) { 
    	push(@input, $_);
    }
close(INPUT) || die "can't close input file: $!";
print @input;
print ">>>>>>>>>>>>>>>>>>  Input file written!  <<<<<<<<<<<<<<<<<<<<<<<\n";

$pid = open2(*Reader, *Writer, "$cmd_subopt");
print Writer @input;
@got = <Reader>;

print @got;

#open(SUBOPT, "$cmd_subopt") || die "can't run subopt: $!";
#    print {SUBOPT} @input;
#    while (<SUBOPT>) {
#    	push(@out_subopt, $_);
#    }
#close(SUBOPT) || die "can't close subopt: $!";
#print @out_subopt;
#print ">>>>>>>>>>>>>>>>>>  output subopt written  <<<<<<<<<<<<<<<<<<<<<<<\n";

#open(BARRS, "$cmd_barriers") || die "can't run barriers: $!";
#	print {BARRS} @out_subopt;
#    	while (<BARRS>) { 
#    		push(@out_barriers, $_);
#    	}
#close(BARRS) || die "can't close barriers: $!";
#print @out_barriers;
#print ">>>>>>>>>>>>>>>>>>  output barriers written  <<<<<<<<<<<<<<<<<<<<<<<\n";

#my @output = `$command`;

#print STDERR "Written to $outfile\n\n" if($verbose);

#open (OUTFILE, '>>.$outfile');
#print OUTFILE @output;
#close (OUTFILE); 

#`xmgrace -log x -nxy $outfile`;



__END__

=head1 NAME

design sequence which folds into specific intermediate states

=head1 SYNOPSIS

run.pl -in I<string> [-out I<string>]

=head1 DESCRIPTION

The program designs a sequence.

=head1 OPTIONS

=over 4

=item B<-in> I<STRING>

Provides the name of the input file

=item B<-out> I<STRING>

Provides the name of the output file (default: out.txt)

=back

=head1 AUTHORS

Stefan Hammer

=head1 BUGS

Please send comments and bug reports to jango@tbi.univie.ac.at

=cut


