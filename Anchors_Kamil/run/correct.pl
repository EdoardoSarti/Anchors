use strict;
use warnings;

open (READ,"<$ARGV[0]");

my $READ_line = <READ>;
$READ_line = <READ>;chomp $READ_line;my @s_ref = split //, $READ_line;
$READ_line = <READ>;
$READ_line = <READ>;chomp $READ_line;my @q_ref = split //, $READ_line;

my $i = 0;
my $s = 0;
my $q = 0;
my %ref;
while (defined $s_ref[$i]) {
	++$s unless ($s_ref[$i] eq "-");
	++$q unless ($q_ref[$i] eq "-");
	$ref{$s} = $q if (($s_ref[$i] ne "-")and($q_ref[$i] ne "-"));
	$ref{$s} = "-1" if (($s_ref[$i] ne "-")and($q_ref[$i] eq "-"));
#	print "$s\t$ref{$s}\t$q\n";
	++$i;
	}
#print "\n";
close READ;

#print "S: $s\n";

open (READ,"<$ARGV[1]");

$READ_line = <READ>;
$READ_line = <READ>;chomp $READ_line;my @s_test = split //, $READ_line;
$READ_line = <READ>;
$READ_line = <READ>;chomp $READ_line;my @q_test = split //, $READ_line;

$i = 0;
$s = 0;
$q = 0;
my $shift = 0;
my $correct = 0;
my $length;
while (defined $s_test[$i]) {
	++$s unless ($s_test[$i] eq "-");
	++$q unless ($q_test[$i] eq "-");
	if ($s_test[$i] ne "-") {
		if ($q_test[$i] ne "-") {
			++$length;
			$shift += abs($ref{$s} - $q);
			if ($ref{$s} eq $q) {
				++$correct;
				}
			}
		}
	++$i;
	}
#print "\n";

close READ;

print "",$shift/$length,"\n",$correct/$length,"\n";




































