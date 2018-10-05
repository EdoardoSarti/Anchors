#!/usr/bin/perl -w

use strict;
use Math::FFT;


use constant window => 10;
use constant ds     => 3;

my $series1;
my $series2;
my $series;
open (READ, "<$ARGV[0]");

my $k = 0;
while (<READ>) {
	chomp;
	my @mas = split/\s+/;
	next unless defined $mas[1];
	next unless defined $mas[2];
	$series1->[$k] = $mas[1];
	$series2->[$k] = $mas[2];
	$series->[$k]  = ($mas[1] + 2) * ($mas[2] + 2) - 4;
	++$k;
	}

close READ;

my $anchors = searchPeaks(convolve(defineCore(window), $series));

print "",scalar @{$anchors},"\n";
foreach my $arg (@{$anchors}) {
	print "",getPos($series1, $arg),"\t",getPos($series2, $arg),"\t-1500.0\n";
	}

#for (my $k = 0; $k < scalar @{$out}; $k++) {
#	print "$k\t$series1->[$k]\t$series2->[$k]\t$out->[$k]\n";
#	}

sub getPos {
	my $series = shift;
	my $pos    = shift;
	my $k = 1;
	for (my $i = 0; $i < scalar @{$series}; $i++) {
		last if $i eq $pos;
		++$k unless $series->[$i] eq '-2';
		}
	return $k;
	}

sub searchPeaks {
	my $series = shift;
	my @peaks;
	for (my $k = ds; $k < scalar @{$series} - ds; $k++) {
		next if $series->[$k] < -2;
		next if $series->[$k] < $series->[$k + ds];
		next if $series->[$k] < $series->[$k - ds];
		$k += ds;
		push (@peaks, $k); 
		}
	return \@peaks;
	}

sub convolve {
	my $core   = shift;
	my $series = shift;
	my $out;
	for (my $k = 0; $k < scalar @{$series}; $k++) {
		$out->[$k] = 0;
		for (my $l = (-1)*int((scalar @{$core})/2); $l < scalar @{$core} - int((scalar @{$core})/2); $l++) {
			next unless defined $series->[$k + $l];
			$out->[$k] = $out->[$k] + $series->[$k + $l]*$core->[$l];
			}
		}
	return $out;
	}

sub defineCore {
	my $window = shift;
	my $core;
	for (my $k = 0; $k < $window; $k++) {
		$core->[$k] = 0.1;
		}
	my $sum = 0;
#	for (my $k = int($window/2); $k < $window; $k++) {
#		$core->[$k] = exp(((-1)*$k) + int($window/2));
#		$sum += $core->[$k];
#		}
	for (my $k = 0; $k < $window; $k++) {
		$sum += $core->[$k];
		}
	for (my $k = 0; $k < $window; $k++) {
		$core->[$k] = $core->[$k] / $sum;
		}
	
	return $core;
	}

