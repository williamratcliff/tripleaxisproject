#!/usr/bin/perl

$file = $ARGV[0];
open(P,$file) or die "Can't open $file ";

# Ignore things in the header
while(<P>) {
    last if /\#Columns/;
}

@head = split;
shift @head;

# Nothing but columns now
while (<P>) {
    chomp;
    @parts = split;
    for ($i=0;$i<=$#parts;$i++) {
	$sum{$head[$i]} += 1.0 * $parts[$i]; # Coerce to be number
    }
}
close(P);
# Spit out the results
foreach $i (0..47) {
    $key = "PSDC".$i;
    #print "$key  $sum{$key}\n";
    printf "%2d\t%6d\n",$i,$sum{$key};
}
