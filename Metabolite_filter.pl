#!/usr/bin/env perl

use strict;
use warnings;

my $csvFile = shift;

my $goodPeakCount= 24;

my $maxQuality = 0.80;


open IN, "$csvFile" or die $!;
open OUT, ">$csvFile.filtered.csv" or die $!;
my $firstLine = <IN>;
#my @tmp = split(",",$firstLine);
#my @r = ($tmp[2],$tmp[9]);
#my @r = ($tmp[9]);
#push @r,@tmp[14 .. @tmp-1];
#my $t = join("\t",@r);
print OUT $firstLine;
#my $n = 0;
while(<IN>)
{
	my @tmp = split(",",$_);
	if($tmp[1] >= $goodPeakCount and $tmp[4] >= $maxQuality)
	{
		print OUT $_;
	}
}
close IN;
close OUT;
#END
