#!/usr/bin/perl -wT

#A short program to integrate Histograms with a given period
#TODO: Various Histogram lengths with a Command line option.
#TODO: Various File Formats
#Written by Thomas Prosser 02/12/2004


my $pi = 4*atan2 1., 1.;
my $file = $ARGV[0];
my @data;
my @vector;
my @keyvaluepairs;
my $sum;
my $sumsum;

open(FILE, $file) or die "File $file not found";
while(<FILE>){
	chomp;
	push @data, $_;
}

my @temp = splice @data, 0, 5;

foreach (@temp){
	chomp;
	push @keyvaluepairs, [split /\s/];
}

my $keynumber = @{$keyvaluepairs[0]};

foreach(@data){
	chomp;
	push @vector, [split /\t/];
}

my $sets = @vector;

for(my $i = 0; $i < $sets; ++$i){
	$sum += $vector[$i][1];
	$sumsum += $vector[$i][1]*$vector[$i][1];
}
$sum/=$sets; $sumsum /= $sets;

printf("%.5f \t %.5f\n", $sum, ($sumsum-$sum*$sum) );

	
