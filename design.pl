#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <2013-01-30 12:53:12 jango>
# Please adopt the path variables to your system.

# path variables
my $bin_subopt = "RNAsubopt";
my $bin_sort = "sort";
my $bin_barriers = "/home/mescalin/jango/Software/Barriers-install/bin/barriers";
my $bin_treekin = "/home/mescalin/jango/Software/treekin-install/bin/treekin";
my $bin_gracebat = "gracebat";



# includes
use RNA;
use Thread;
use Math::Complex;
use Pod::Usage;
use Getopt::Long;

my $man = 0;
my $help = 0;
my $paralell = '4';
my $opt_subopt = "-e 20";
my $opt_barriers = "--bsize --max 1000";
my $opt_treekin = "-m I";
my $treshold = '0.1';
my $noLP = 1;

# option parser
&usage() unless GetOptions(	'treshold=f' => \$treshold,
				'noLP' => \$noLP,
				'subopt=s' => \$opt_subpot,
				'barriers=s' => \$opt_barriers,
				'treekin=s' => \$opt_treekin,
				'paralell=s' => \$paralell,
				'help|?' => \$help,
				man => \$man) or pod2usage(2);
# Parse options and print usage if there is a syntax error, or if usage was explicitly requested.
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# parse options
if ($noLP) {
	$opt_subopt += " --noLP";
	$opt_barriers += " -G RNA-noLP"; 
} else { 
	$opt_barriers += " -G RNA"; 
}


# Shared data containers
my @subopt: shared;
my @barriers: shared;
my @treekin: shared;
my @sort: shared;
my @evaluate: shared;
my %sequences: shared;
my $id: shared = 1;
our %desired_structures;


my $interactiv = init_ia(0);
my $ia = 1 if -t STDIN && -t STDOUT && $#ARGV < 0;
$interactiv->() if $ia;
process_input();

# interactive interface for structure input
sub init_ia {
   my $times = 0;
   return sub {
      if ($times == 0) {
         print "Input start-structure, intermediate-structure and end-structure\nwith their desired inflection times (space separated).\n";
         print "e.g.: ((((((((....)))))))) 0 0.333\n";
         print "    @ to quit, enter for random start string\n";
      }
      for (1..8) { print '....,....', $_;}
      print "\n";
      $times=1;
   }
}

# read input file with structures (sequence, rel time first inflection point, second point)
sub process_input {
	$_ = <>;
	chomp;
	$desired_structures{"start"} = 	[split(/\s+/,$_)];
	chomp($_ = <>);
	$desired_structures{"inter"} = 	[split(/\s+/,$_)];
	chomp($_ = <>);
	$desired_structures{"end"} = 	[split(/\s+/,$_)];
	
	die "structures have unequal length"
		if ((length($desired_structures{"start"}[0]) != length($desired_structures{"inter"}[0])) || 
		(length($desired_structures{"inter"}[0]) != length($desired_structures{"end"}[0])));
	print "Your input is:\n";
	for ("start", "inter", "end") {
		print $_."\t".$desired_structures{$_}[0]."\t".$desired_structures{$_}[1]."\t".$desired_structures{$_}[2]."\n";
		die "$_ structure: relative times of inflection points missing or invalid value!"
			if (($desired_structures{$_}[1] < 0) || ($desired_structures{$_}[2] < 0) || 
			($desired_structures{$_}[1] > 1) || ($desired_structures{$_}[2] > 1))
	}
}

# add a switch.pl sequence for the fist iteration
add_sequence("0", get_switch_sequence());

# Start the queuing system
my $queue = Thread->new(\&queue_manager, $paralell);


END {
	$queue->join();
}

# The queue manager checks every queue and starts the jobs
sub queue_manager {
	my $end = 0;
	my %thisistheend;
	
	($p) = @_;
	for (1 .. $p) {
		${"job".$_} = Thread->new(\&idle, $_);
	}
	# queue manager should always check if processes in threads are still running.
	while($end < $p) {
		for (1 .. $p) {
			if (${"job".$_}->done) {
				# if one job is finished, evaluate the shortest queue and start a new job!
				#print "Job number $_ is already done.\n";
				${"job".$_}->join;
				
				lock(@subopt,@barriers,@treekin,@sort,@evaluate);
				my %num = ("subopt", scalar(@subopt),
					   "barriers", scalar(@barriers),
					   "treekin", scalar(@treekin),
					   "sort", scalar(@sort),
					   "evaluate", scalar(@evaluate));
				# show how many items are in the queues at the moment
				print "$_: Queue:\t";
				foreach (keys %num) { print $_." ".$num{$_}."\t";};
				print "\n";
				my ($key) = sort { $num{$b} cmp $num{$a} } keys %num;
				# start a job from the longest queue
				if ($num{$key} != 0) {
					${"job".$_} = Thread->new(\&{"start_".$key});
					$thisistheend{$_} = 0;
				# or start an "idle" job doing nothing
				} else {
					${"job".$_} = Thread->new(\&idle, $_);
					$thisistheend{$_} = 1;
				}
			} else {
				#print "$_: job is still running...\n";
			}
			sleep 1;
		}
		$end = 0;
		foreach (values %thisistheend) { $end += $_; }
	}
	# in the end join all processes
	for (1 .. $paralell) {
		${"job".$_}->join;
	}
}

# job is idle. important if nothing has to be run
sub idle {
	($cid) = @_;
	print "$cid: is idle.\n";
}

# sub function to run RNAsubopt
sub start_subopt {
	my $cid;
	{
		lock(@subopt);
		$cid = shift(@subopt);
	}
	warn("cid is undef in @{[(caller($cid))[3]]}"), return undef unless defined($cid);
	print "start_subopt with id: $cid.\n";
	system("mkdir $cid");
	system("echo '$sequences{$cid}' > '$cid/$names{$cid}.fa'");
	system("echo '$bin_subopt $opt_subopt < '$cid/$names{$cid}.fa' > $cid/subopt.out' >> $cid/commands.txt");
	system("$bin_subopt $opt_subopt < '$cid/$names{$cid}.fa' > $cid/subopt.out");
	if ($? == -1) {
		print "failed to execute: $!\n";
	} elsif ($? & 127) {
		printf "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
	} else {
		printf "child exited with value %d\n", $? >> 8;
	}
	
	lock(@sort);
	push(@sort, $cid);
}

# sub function to run sort
sub start_sort {
	my $cid;
	{
		lock(@sort);
		$cid = shift(@sort);
	}
	warn("cid is undef in @{[(caller($cid))[3]]}"), return undef unless defined($cid);
	print "start_sort with id: $cid.\n";
	system("mkdir $cid/tmp");
	system("echo '$bin_sort -k2 -n $cid/subopt.out -T $cid/tmp > $cid/sort.out' >> $cid/commands.txt");
	system("$bin_sort -k2 -n $cid/subopt.out -T $cid/tmp > $cid/sort.out");
	if ($? == -1) {
		print "failed to execute: $!\n";
	} elsif ($? & 127) {
		printf "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
	} else {
		printf "child exited with value %d\n", $? >> 8;
	}
	
	lock(@barriers);
	push(@barriers, $cid);
}

# sub function to run barriers
sub start_barriers {
	my $cid;
	{
		lock(@barriers);
		$cid = shift(@barriers);
	}
	warn("cid is undef in @{[(caller($cid))[3]]}"), return undef unless defined($cid);
	print "start_barriers with id: $cid.\n";
	system("echo 'cd $cid; $bin_barriers $opt_barriers --rates < sort.out > barriers.out; cd ..;' >> $cid/commands.txt");
	system("cd $cid; $bin_barriers $opt_barriers --rates < sort.out > barriers.out; cd ..;");
	if ($? == -1) {
		print "failed to execute: $!\n";
	} elsif ($? & 127) {
		printf "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
	} else {
		printf "child exited with value %d\n", $? >> 8;
	}
		
	lock(@treekin);
	push(@treekin, $cid);
}

# sub function to run treekin
sub start_treekin {
	my $cid;
	{
		lock(@treekin);
		$cid = shift(@treekin);
	}
	warn("cid is undef in @{[(caller($cid))[3]]}"), return undef unless defined($cid);
	print "start_treekin with id: $cid.\n";
	
	# get id of start structure and start kinetics there
	my $start = [get_structure_ids($cid)]->[0];
	if ($start == "") { $start = 1; warn "Could not find ID of start structure for cid: $cid! Results are WRONG!"; }
	$opt_treekin += " --p0 $start=1";
	
	# run treekin
	system("echo '$bin_treekin $opt_treekin --ratesfile=$cid/rates.out < $cid/barriers.out > $cid/treekin.out' >> $cid/commands.txt");
	system("$bin_treekin $opt_treekin --ratesfile=$cid/rates.out < $cid/barriers.out > $cid/treekin.out");
	if ($? == -1) {
		print "failed to execute: $!\n";
	} elsif ($? & 127) {
		printf "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
	} else {
		printf "child exited with value %d\n", $? >> 8;
	}

	# make beautiful graphs with gracebat
	system("echo '$bin_gracebat -log x -nxy $cid/treekin.out -hdevice PostScript -hardcopy -printfile $cid/treekin.ps' >> $cid/commands.txt");
	system("$bin_gracebat -log x -nxy $cid/treekin.out -hdevice PostScript -hardcopy -printfile $cid/treekin.ps");
	lock(@evaluate);
	push(@evaluate, $cid);
}

# evaluation of your sequence
sub start_evaluate {
	my $cid;
	{
		lock(@evaluate);
		$cid = shift(@evaluate);
	}
	warn("cid is undef in @{[(caller($cid))[3]]}"), return undef unless defined($cid);
	
	# curve scetching and calculate fitness function
	$score = analyse_curve($cid);
	lock(%sequences);
	my @sequence = (split(';',$sequences{$cid}));
	my $father = $sequence[1];
	$sequences{$cid} = join(';', ($score,$father,$sequence[2]));
	
	open (MYFILE, '>solution_tree.out');
	print MYFILE $_."\t".$sequences{$_}."\n", 
		for (sort {$sequences{$b} <=> $sequences{$a}} keys %sequences);
	close (MYFILE);
	print "ID $cid is evaluated and completed! Score: $score\n";
	
	# decide and mutate or output!
	if ($score >= [split(';',$sequences{$father})]->[0]) {
		add_sequence($cid, mutate_sequence());
		add_sequence($cid, mutate_sequence());
		add_sequence($cid, mutate_sequence());
		add_sequence($cid, mutate_sequence());
	} else {
		#this branch is dead!
	}
}

# read kinetics from file and get special points
# whole sub to analyse the curve with id $cid
sub analyse_curve {
	my ($cid) = @_;
	my $score = 0;
	
	my @treekin;
	open TREEKIN, "<$cid/treekin.out";
	while (<TREEKIN>)
	{
		chomp;
		if (substr($_,0,1) ne '#') {
			push(@treekin, [split(' ', $_)]);
			
		}
	}
	close TREEKIN;
	#print $treekin[1]->[3]."\n";
	#foreach (@{$treekin[2]}) {
	#print $_."\n";
	#}
	
	# get ids of our desired structures
	my ($start,$intermediate,$end) = get_structure_ids($cid);
	print "s,i,e ids: ".$start."\t".$intermediate."\t".$end."\n";
	
	#Get max points and equilibrium distribution
	my %maximum;
	my %equilibrium_dist;
	my $counttimes = scalar(@treekin);
	#iterate over all values and remember the maximal points
	for (my $j = 0; $j < $counttimes; $j++) {
		for (my $i = 1; $i < scalar(@{$treekin[0]}); $i++) {
			#print $i." i\t".$j." j\t".$treekin[$j]->[$i]." value\n";
			if ($maximum{($i)} < $treekin[$j]->[$i]) {
				# $i because first structure is always labeled 1 (=MFE)
				$maximum{($i)} = $treekin[$j]->[$i];
			}
			
			#remember values in last column = equilibrium distribution
			if ($j == $counttimes-1) {
				$equilibrium_dist{$i} = $treekin[$j]->[$i];
			}
		}
	}
	
	# get inflection point x-values
	# iterate overall maximum points and define treshold
	my %inflection_points;
	
	for (my $i = 1; $i <= keys(%maximum); $i++) {
		# TODO: just need inflection points of desired structures?!
		if (($maximum{$i} >= $treshold) || ($i==$start) || ($i==$intermediate) || ($i==$end)) {
			print "\n".$i."\t".$maximum{($i)}."\t".$equilibrium_dist{($i)}."\n";
			# generate array with single curve
			my @curve;
			for (my $j = 0; $j < $counttimes; $j++) {
				push(@curve, $treekin[$j]->[$i]);
			}
			$inflection_points{$i} = get_inflection_points(@curve);
			
			foreach (@{$inflection_points{$i}}) {
				print $_->[0]."\t".int($treekin[$_->[0]]->[0])."\t".$_->[1]."\n";
			}
		}
	}
	
	# calculate fittness function score
	# start structure(height, inflection point down, equil. dist height)
	# intermediate str(height, inflection point up and down, equili. dist height)
	# end structure (height/equili. dist height, inflection point up
	# max points:
	# $maximum{$i}
	# desired structures
	# $start,$intermediate,$end
	# inflection points (0-~x-value) (1-up/down)
	# $inflection_points{$i}->[0]	$inflection_points{$i}->[1]
	# calculate the x values in percent: 
	# $value/$counttimes
	# get equilibrium distribution
	# $equilibrium_dist{$i}
	# desired first inflection point of structures: start,inter,end
	# $desired_structures{"start"}->[1]	$desired_structures{"inter"}->[1]	$desired_structures{"end"}->[1]
	# desired second inflection point of structure 0,1,2 = start,intermed,end
	# $desired_structures{"start"}->[2]
	lock(%desired_structures);
	
	# iterate over curves and add scores
	for (my $i = 1; $i < scalar(@{$treekin[0]}); $i++) {
		# if curve is one of the desired ones, do something special
		if ($i == $start) {
			$score += $maximum{$i};
			$score -= $equilibrium_dist{$i};
			# take first inflection point from a down curve
			foreach (@{$inflection_points{$i}}) {
				if ($_->[1] eq "down") {
					$score -= abs($desired_structures{"start"}->[2] - $_->[0]/$counttimes);
					last;
				}
			}
		} elsif ($i == $intermediate) {
			$score += $maximum{$i};
			$score -= $equilibrium_dist{$i};
			# take first inflection point from a up curve
			foreach (@{$inflection_points{$i}}) {
				if ($_->[1] eq "up") {
					$score -= abs($desired_structures{"start"}->[2] - $_->[0]/$counttimes);
					last;
				}
			}
			# take last inflection point from a down curve
			for ($n = scalar(@{$inflection_points{$i}}); $n >= 0; $n--) {
				if (${$inflection_points{$i}}[$n]->[1] eq "down") {
					$score -= abs($desired_structures{"end"}->[2] - ${$inflection_points{$i}}[$n]->[0]/$counttimes);
					last;
				}
			}
			
		} elsif ($i == $end) {
			$score += $maximum{$i};
			$score += $equilibrium_dist{$i};
			# take last inflection point from a upwards curve
			for ($n = scalar(@{$inflection_points{$i}}); $n >= 0; $n--) {
				if (${$inflection_points{$i}}[$n]->[1] eq "up") {
					$score -= abs($desired_structures{"end"}->[2] - ${$inflection_points{$i}}[$n]->[0]/$counttimes);
					last;
				}
			}
		} else {
			if ($maximum{$i} >= $treshold) {
				$score -= $equilibrium_dist{$i};
				$score -= $maximum{$i};
			}
		}
	}
	
	return $score;
}

# gets the ids of our desired structures from the barriers.out file
sub get_structure_ids {
	my ($cid) = @_;
	my (%structures,$start,$intermediate,$end);
	open BARRIERS, "<$cid/barriers.out";
	while (<BARRIERS>)
	{
		chomp;
		$_ =~ s/^\s+//;
		my ($key, $val) = split / /;
		$structures{$key} = $val;
	}
	close BARRIERS;
	
	# find desired structures in our hash.
	lock(%desired_structures);
	for (keys %structures) {
		$start = $_ if ($structures{$_} eq $desired_structures{"start"}->[0]);
		$intermediate = $_ if ($structures{$_} eq $desired_structures{"inter"}->[0]);
		$end = $_ if ($structures{$_} eq $desired_structures{"end"}->[0]);
	}
	return ($start,$intermediate,$end);
}

# get inflection point of a curve
sub get_inflection_points {
	# function takes a array of numbers (one curve only)
	my @curve = @_;
	my @inflection_points;
	my $counttimes = scalar(@curve);
	my @delta1;
	my @delta2;
	
	# calculate first derivative
	for (my $j = 0; $j < scalar(@curve)-1; $j++) {
		#print "$j curve is:\t".$curve[$j]."\n";
		push(@delta1, ($curve[$j+1] - $curve[$j]));
	}
	# calculate second derivative
	for (my $j = 0; $j < scalar(@delta1)-1; $j++) {
		#print "$j delta1 is:\t".$delta1[$j]."\n";
		push(@delta2, ($delta1[$j+1] - $delta1[$j]));
	}
	
	# calculate mean value of the first derivative
	my $sum;
	foreach (@delta1) {
		$sum += $_;
	}
	my $mean_d1 = ($sum / scalar(@delta1));
	print "mean of delta1 is: ".$mean_d1."\n";
	
	# calculate standard deviation from first derivative
	my $intermediate;
	foreach (@delta1) {
		$intermediate += ($_-$mean_d1)*($_-$mean_d1);
	}
	my $std_deviation = sqrt($intermediate / scalar(@delta1));
	print "std deviation of delta1 is: ".$std_deviation."\n";
	
	# treshold to exclude too flat inflection points (e.g.: in a horizontal line)
	# treshold is mean plus standard deviation for delta1
	$zerotreshold = $mean_d1 + $std_deviation;
	print "zerotreshold is: ".$zerotreshold."\n";
	for (my $j = 0; $j < scalar(@delta2)-1; $j++) {
		#print "$j delta2 is:\t".$delta2[$j]."\n";
		if (abs($delta1[$j]) >= $zerotreshold) {
			if (($delta2[$j] > 0) && ($delta2[$j+1] < 0)) {
				push(@inflection_points, [$j,"up"]);
				print $delta2[$j]."\t> 0 >\t".$delta2[$j+1]."\n";
			} elsif (($delta2[$j] < 0) && ($delta2[$j+1] > 0)) {
				push(@inflection_points, [$j,"down"]);
				print $delta2[$j]."\t< 0 <\t".$delta2[$j+1]."\n";
			} elsif (($delta2[$j-1] < 0) && ($delta2[$j] == 0) && ($delta2[$j+1] > 0)) {
				push(@inflection_points, [$j,"down"]);
				print $delta2[$j-1]."\t< 0=j <\t".$delta2[$j+1]."\n";
			} elsif (($delta2[$j-1] > 0) && ($delta2[$j] == 0) && ($delta2[$j+1] < 0)) {
				push(@inflection_points, [$j,"up"]);
				print $delta2[$j-1]."\t> 0=j >\t".$delta2[$j+1]."\n";
			}
		}
	}

	
	#foreach (@inflection_points) {
		#print $_->[0]."\t".$_->[1]."\n";
	#}

	
	return \@inflection_points;
}


# add a sequence to the queue. remember its sequence, father and id in a shared hash "sequences" ($cid)[score;father;sequence]
sub add_sequence {
	my ($father, $s) = @_;
	lock(%sequences, $id, @subopt);
	$sequences{$id} = join(';', (0,$father,$s));

	push(@subopt, $id);
	print "Added new sequence: id \t$id\tfather $father\n$s\n";
	$id++;
}

# get switch.pl sequence
sub get_switch_sequence {
	my $length = length($desired_structures{"start"}[0]);
	my $sequence = `perl rseq.pl -l $length`;
	return $sequence;
}

# mutate the sequence
sub mutate_sequence {
	
	
	
	my $length = length($desired_structures{"start"}[0]);
	my $sequence = `perl rseq.pl -l $length`;
	return $sequence;
}



__END__

=head1 NAME

Pipline to calculate RNAsubopt -> Sort -> Barriers -> Treekin and draw a beautiful chart.

=head1 SYNOPSIS

treekin_pipline.pl [options] Sequence-File

Try treekin_pipline.pl --man for more information!

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-paralell=(string)>

Defines how many calculations should be done in paralel (default: 4)

=item B<-subopt=(string)>

Command-line options for subopt (default: "-e 20 --noLP")

=item B<-barriers=(string)>

Command-line options for barriers (default: "-G RNA-noLP --bsize --max 1000")

=item B<-treekin=(string)>

Command-line options for treekin (default: "-m I --p0 1=1")

=back

=head1 DESCRIPTION

This program reads Sequences and calcuclates the Folding Kinetics using RNAsubopt, Barriers and Treekin.

=cut

