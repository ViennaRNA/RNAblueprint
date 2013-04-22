#!/usr/bin/perl
use Math::Complex;


my $treshold = '0.1';
# holds the secondary structure, the x-value in percent of the first inflection point and the second inflection point
my %desired_structures: shared;
$desired_structures{"start"} = 	["......((((((......))))))(((....)))((((......))))........",0,0.333];
$desired_structures{"inter"} = 	["......((((((......))))))........((((((......))))))......",0.333,0.666];
$desired_structures{"end"} = 	["((((((((((((((((((((................))))))))))))))))))))",0.666,0];

print "score is: ".analyse_curve('4')."\n";



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
	lock(%desired_structures);
	# $desired_structures{"start"}->[1]	$desired_structures{"inter"}->[1]	$desired_structures{"end"}->[1]
	# desired second inflection point of structure 0,1,2 = start,intermed,end
	# $desired_structures{"start"}->[2]
	
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

