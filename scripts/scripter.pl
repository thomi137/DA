#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my $shebang="#!/bin/bash\n";
my $program_path = "$ENV{'HOME'}/DA/codes/";
my $output_path = "$ENV{'HOME'}/DA/data/";

our ($opt_s, $opt_w, $opt_p, $opt_n);
getopts('sw:p:n:');

print "Walltime: ".$opt_w ."\n" if $opt_w or die "No Walltime set";
print "ppn: ".$opt_p."\n" if $opt_p or die "No ppn set";
print "n: ".$opt_n."\n" if $opt_n or die "No ppn set";

my ($program, $output, $scriptfile) = @ARGV;
my $jobdetails = "#PBS -l nodes=".$opt_n.":ppn=".$opt_p.",walltime=".$opt_w."\n";
my $job_to_run = "$program_path$program > $output_path$output\n";

open(FILE,"> $scriptfile");
print FILE $shebang;
print FILE $jobdetails;
print FILE "#job_m\n";
print FILE "#\n";
print FILE $job_to_run;
close(FILE);

chmod(0755, $scriptfile);

`qsub $scriptfile` if $opt_s;
