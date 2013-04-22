#!/usr/bin/perl
# -*-Perl-*-


system("ls -al >> test.out");
system("ls >> test.out");
     if ($? == -1) {
    print "failed to execute: $!\n";
    }
    elsif ($? & 127) {
    printf "child died with signal %d, %s coredump\n",
    ($? & 127), ($? & 128) ? 'with' : 'without';
    }
    else {
    printf "child exited with value %d\n", $? >> 8;
    }



    push() - adds an element to the end of an array.
    unshift() - adds an element to the beginning of an array.
    pop() - removes the last element of an array.
    shift() - removes the first element of an array.
push(@coins, "Penny");

print scalar(@nums)."<br />";
print scalar(@alpha)."<br />";

# REDEFINE TO SCALAR
$nums = @nums;
$alpha = @alpha;

# SORT 'EM
@foods = sort(@foods);
@Foods = sort(@Foods);

# DEFINE A HASH
%coins = ( "Quarter" , 25,
           "Dime" ,    10,
           "Nickel",    5 );
# FOREACH LOOP
foreach $key (sort keys %coins) {
     print "$key: $coins{$key}<br />";
}


# DEFINE A HASH
%coins = ( "Quarter" , .25,
           "Dime" , .10,
           "Nickel", .05 );
# FOREACH LOOP
foreach $value (sort {$coins{$a} cmp $coins{$b} }
           keys %coins)
{
     print "$value $coins{$value}<br />";
}

# ADD NEW ELEMENT PAIRS
$coins{Penny} = .01;
$coins{HalfDollar} = .50;

# DELETE THE ELEMENT PAIRS
delete($coins{Penny});
delete($coins{HalfDollar});


# GET STRUCTURE OF DATA
use Data::Dumper

Dumper(\@array);
