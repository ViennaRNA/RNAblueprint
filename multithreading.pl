#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <2012-12-20 12:28:16  jango>
#
# Program counting up and test for primes using multithreaded programming
# As the only shared variable is like a pinhole, this program is not really calculating primes in paralell.
# However it shows very nice how multithreading can be done!

use Thread;

# only one variable which is shared between the threads
my $number: shared = 1;

# start n threads with a unique name and pass the id along as an option
for (1..100) {
	${"thread".$_} = Thread->new(\&count, $_);
}
# join the threads again in the end so that the program is able to end
for (1..100) {
	${"thread".$_}->join();
}

#main program ends here!


#subroutine which counts up in every thread
#this routine also starts the primefinder
sub count {
	for (1..20000) {
		#print thread-id times a space
		for (0..$_[0]) { print " "; }
		#print thread id
		print "($_[0]) ";
		# lock shared variable, otherwise funny stuff happens!
		lock($number);
		print $number;
		# test if number is a prime
		if (testprime($number)) { print ": >>>>>>>> is a PRIME! <<<<<<<<";}
		print "\n";
		#count up!
		$number++;
	}
}

# subrouitine which tests if the number is a prime - primefinder
sub testprime {
	my $test = $_[0];
	my $i = 2;
	while ($i < $test) {
		# a prime should never return 0 when number modolo something smaller
		return 0 unless ($test % $i++);
	}
	# return 1 if it is a prime
	return 1;
}

__END__
