#!/usr/bin/perl -wT

#A short program to integrate Histograms with a given period
#TODO: Various Histogram lengths with a Command line option.
#TODO: Various File Formats
#Written by Thomas Prosser 02/12/2004


use strict;

use warnings;
my $pi = 4*atan2 1., 1.;
my $file = $ARGV[0];
my @data;
my @vector;
my @keyvaluepairs;
my %parameters;
my @calculation;
my @means;
my @variances;
my $periods;

open(FILE, $file) or die "File $file not found";
while(<FILE>){
	chomp;
	push @data, $_;
}

@vector=@data;
my @temp = splice @vector, 0, 3;

foreach (@temp){
	chomp;
	push @keyvaluepairs, [split /\t/];
}
foreach (@{$keyvaluepairs[0]}){
	tr/a-z/A-Z/;
}
my $keynumber = @{$keyvaluepairs[0]};

for(my $i=0; $i<$keynumber; ++$i){
	$parameters{$keyvaluepairs[0][$i]}=$keyvaluepairs[1][$i];
}

#up to here, we saved the parameters, the vector is ready to be taken apart.

while(@vector){
	shift @vector;
	my @workspace = splice @vector, 0, $parameters{"POINTS"};	
	push @means, mean(@workspace);
	push @variances, sqrt(variance(@workspace));
}

$periods = @means;
printf("Points: %d\n", $parameters{"POINTS"});
printf("Kick strength: %d\n", $parameters{"KICK STRENGTH"});
printf("Coupling: %d\n", $parameters{"COUPLING"});
printf("Periods: %d\n", $periods);
print "\n";
for(my $k=0; $k<$periods; ++$k){
	printf("%.9f\t%.9f\n", $means[$k], $variances[$k]);
}

###############Subroutines############################################

sub mean{
	local @points=@_;
	local $norm=0;
	#first the normalisation of the means...
	for(my $j=0; $j<$parameters{"POINTS"}; ++$j){
		$norm+=$points[$j];
	}
	#then the weighted sum...			
	local $sum=0;
	for(my $i = 0; $i< $parameters{"POINTS"}; ++$i){
		local $momentum = -1*$pi*$parameters{"POINTS"}/$parameters{"SYSTEM SIZE"}+$i*2.*$pi/$parameters{"SYSTEM SIZE"};
		$sum+=$points[$i]*$momentum;
	}
	$sum/=$norm;
	return $sum;
}

sub meansq{
	local @points=@_;
	local $norm=0;
	#first the normalisation of the means...
	for(my $j=0; $j<$parameters{"POINTS"}; ++$j){
		$norm+=$points[$j];
	}
	#then the weighted sum...			
	local $sum=0;
	for(my $i = 0; $i< $parameters{"POINTS"}; ++$i){
		local $momentum = -1*$pi*$parameters{"POINTS"}/$parameters{"SYSTEM SIZE"}+$i*2.*$pi/$parameters{"SYSTEM SIZE"};
		$sum+=$points[$i]*$momentum*$momentum;
	}
	$sum/=$norm; 
	return $sum;
}

sub variance{
	local @points=@_;
	local $me=mean(@points);
	return meansq(@points)-$me*$me;
}
