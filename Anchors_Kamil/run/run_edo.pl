#!/usr/bin/perl -w

use strict;

my $length = 10;
my $start = 0;
my $max = 1000.5;
my $min = 1000;
my $number = 20;
my $cycles = 1;
my $score = "correct.pl";
my $skaPath = "ska.fasta";
my $anchPath = "anch.positions.v2";

my %ac_number =("C" => "1", "A" => "5", "E" => "9",  "K" => "13", "V" => "17",
                "S" => "2", "G" => "6", "Q" => "10", "M" => "14", "F" => "18",
                "T" => "3", "N" => "7", "H" => "11", "I" => "15", "Y" => "19",
                "P" => "4", "D" => "8", "R" => "12", "L" => "16", "W" => "20"
                );
my %blosum62 = (
                "C" => ['9' ,'-1','-1','-3','0' ,'-3','-3','-3','-4','-3','-3','-3','-3','-1','-1','-1','-1','-2','-2','-2'],
                "S" => ['-1','4' ,'1' ,'-1','1' ,'0' ,'1' ,'0' ,'0' ,'0' ,'-1','-1','0' ,'-1','-2','-2','-2','-2','-2','-3'],
                "T" => ['-1','1' ,'5' ,'-1','0' ,'-2','0' ,'-1','-1','-1','-2','-1','-1','-1','-1','-1','0' ,'-2','-2','-2'],
                "P" => ['-3','-1','-1','7' ,'-1','-2','-2','-1','-1','-1','-2','-2','-1','-2','-3','-3','-2','-4','-3','-4'],
                "A" => ['0' ,'1' ,'0' ,'-1','4' ,'0' ,'-2','-2','-1','-1','-2','-1','-1','-1','-1','-1','0' ,'-2','-2','-3'],
                "G" => ['-3','0' ,'-2','-2','0' ,'6' ,'0' ,'-1','-2','-2','-2','-2','-2','-3','-4','-4','-3','-3','-3','-2'],
                "N" => ['-3','1' ,'0' ,'-2','-2','0' ,'6' ,'1' ,'0' ,'0' ,'1' ,'0' ,'0' ,'-2','-3','-3','-3','-3','-2','-4'],
                "D" => ['-3','0' ,'-1','-1','-2','-1','1' ,'6' ,'2' ,'0' ,'-1','-2','-1','-3','-3','-4','-3','-3','-3','-4'],
                "E" => ['-4','0' ,'-1','-1','-1','-2','0' ,'2' ,'5' ,'2' ,'0' ,'0' ,'1' ,'-2','-3','-3','-2','-3','-2','-3'],
                "Q" => ['-3','0' ,'-1','-1','-1','-2','0' ,'0' ,'2' ,'5' ,'0' ,'1' ,'1' ,'0' ,'-3','-2','-2','-3','-1','-2'],
                "H" => ['-3','-1','-2','-2','-2','-2','1' ,'-1','0' ,'0' ,'8' ,'0' ,'-1','-2','-3','-3','-3','-1','2' ,'-2'],
                "R" => ['-3','-1','-1','-2','-1','-2','0' ,'-2','0' ,'1' ,'0' ,'5' ,'2' ,'-1','-3','-2','-3','-3','-2','-3'],
                "K" => ['-3','0' ,'-1','-1','-1','-2','0' ,'-1','1' ,'1' ,'-1','2' ,'5' ,'-1','-3','-2','-2','-3','-2','-3'],
                "M" => ['-1','-1','-1','-2','-1','-3','-2','-3','-2','0' ,'-2','-1','-1','5' ,'1' ,'2' ,'1' ,'0' ,'-1','-1'],
                "I" => ['-1','-2','-1','-3','-1','-4','-3','-3','-3','-3','-3','-3','-3','1' ,'4' ,'2' ,'3' ,'0' ,'-1','-3'],
                "L" => ['-1','-2','-1','-3','-1','-4','-3','-4','-3','-2','-3','-2','-2','2' ,'2' ,'4' ,'1' ,'0' ,'-1','-2'],
                "V" => ['-1','-2','0' ,'-2','0' ,'-3','-3','-3','-2','-2','-3','-3','-2','1' ,'3' ,'1' ,'4' ,'-1','-1','-3'],
                "F" => ['-2','-2','-2','-4','-2','-3','-3','-3','-3','-3','-1','-3','-3','0' ,'0' ,'0' ,'-1','6' ,'3' ,'1'],
                "Y" => ['-2','-2','-2','-3','-2','-3','-2','-3','-2','-1','2' ,'-2','-2','-1','-1','-1','-1','3' ,'7' ,'2'],
                "W" => ['-2','-3','-2','-4','-3','-2','-4','-4','-3','-2','-2','-3','-3','-1','-3','-2','-3','1' ,'2' ,'11'],
                );
sub sp {
	my @s = split//, shift;
	my @q = split//, shift;
	my $i = 0;
	my $sp = 0;
	my $gap = 0;
	my $al = 0;
	while (defined $s[$i]) {
		if (($s[$i] eq "-") or ($q[$i] eq "-")) {
			++$gap;
			} else {
			++$al;
			$sp += ($blosum62{"$s[$i]"})->[$ac_number{"$q[$i]"} - 1];
			}
		++$i;
		}
	return $sp;
	}

sub getseq {
#	my @seq = split//, shift;
#	my $position = shift;
#	my $result = '';
#	my $n = 0;
#	my $i = 0;
#	my $flag = 0;
#	OUTER: while (defined($seq[$i])) {
#		++$n unless ($seq[$i] eq '-');
#		if ($n eq $position) {
#			my $j = $i - $length;
#			while (defined($seq[$j])) {
#				last OUTER if $j > $i +  $length;
#				$result = $result . $seq[$j];
#				++$j;
#				}
#			}
#		++$i;
#		}
#	return $result;
	my @seq1 = split//, shift;
	my @seq2 = split//, shift;
	my $pos1 = getpos(join("",@seq1), shift);
	my $pos2 = getpos(join("",@seq2), shift);
	my @s = @seq1[(int(($pos1 + $pos2) / 2) - 10)..(int(($pos1 + $pos2) / 2) + 10)];
	my @q = @seq2[(int(($pos1 + $pos2) / 2) - 10)..(int(($pos1 + $pos2) / 2) + 10)];
	return [join("", @s), join("", @q)];
	}

sub getpos {
	my @seq = split//, shift;
	my $position = shift;
	my $n = 0;
	my $i = 0;
	while (defined ($seq[$i])) {
		++$n unless ($seq[$i] eq '-');
		if ($n eq $position) {
			return $i;
			}
		++$i;
		}
	}

sub print_anch {
	my @positions = @{$_[0]};
	open (WRITE,">anchors");
	print WRITE "$number\n";
	foreach my $arg (@positions) {
		print WRITE join("\t",@{$arg}),"\n";
		}
	close WRITE;
	}

sub getMaxMin {
	my $seq1 = shift;
	my $seq2 = shift;
	my @s1 = split//, $seq1;
	my @s2 = split//, $seq2;
	die unless length($seq1) eq length($seq2);
	my $max = -1000;
	my $min = 1000;
	for (my $i = $length + 1; $i < length($seq1) - $length - 1; $i++) {
		my $s = join("",@s1[($i - 10)..($i + 10)]);
		my $q = join("",@s2[($i - 10)..($i + 10)]);
		my $sp = sp($s, $q);
		if ($sp > $max) {$max = $sp}
		if ($sp < $min) {$min = $sp}
		}
	return [$max, $min];
	}

my %history;

sub set_anchors {
	my @positions = @{$_[0]};
	open (READ,"<out.fasta");
	my $seq1 = <READ>;
	$seq1 = <READ>; chomp $seq1;
	my $seq2 = <READ>;
	$seq2 = <READ>; chomp $seq2;
	close READ;

	my $i = 0;
	my $limit = getMaxMin($seq1, $seq2);
	while (defined($positions[$i])) {
		my $s = getseq($seq1, $seq2, ($positions[$i])->[0], ($positions[$i])->[1]);
		my $sp = sp(($s)->[0], ($s)->[1]);
		($positions[$i])->[2] = ($max - $min) * (($sp - ($limit)->[1])/(($limit)->[0] - ($limit)->[1])) + $min;
		$history{($positions[$i])->[0]} = [@{$history{($positions[$i])->[0]}},$sp, ($positions[$i])->[2]];
		++$i;
		}
	}

sub launch {
	`../../anchors_lucy/alignme1.2/alignme1.2.exe -fasta_file1 betp -fasta_file2 leut -similarity_score_file scale.txt -output_aligned_sequences out -anchors anchors -gap_opening_penalty 21 -gap_extension_penalty 2.1 -termini_gap_opening_penalty 0 -termini_gap_extension_penalty 1`;
	`echo ">betp" > out1.fasta`;
#	`grep -v "*" out | grep -v "Align" | sed '/^\\s*\$/d' | awk '{print \$2}' | sed 'n; d' >> out1.fasta`;
	`grep -v "*" out | grep -v "Align" |  awk '{FS=""; if(\$1!="") {if(a) printf "%s", a}; FS=" "} {a=\$2} END {print}' >> out1.fasta`;
	`echo ">leut" > out2.fasta`;
#	`grep -v "*" out | grep -v "Align" | sed '/^\\s*\$/d' | awk '{print \$2}' | sed '1d; n; d' >> out2.fasta`;
        `grep -v "*" out | grep -v "Align" |  awk '{FS=""; if(\$1=="") {if(a) printf "%s", a}; FS=" "} {a=\$2} END {print}' >> out2.fasta`;
	`cat out1.fasta out2.fasta > out.fasta`;
#	`sed -e 's/^>/|>/g' out.fasta | tr "\n" ">" | tr "|" "\n" | sed -e 's/^>/|/g' | sed -e \$'s/>/\\\\\\n/' | sed -e 's/>//g' | sed -e 's/|/>/g' | sed -e /./!d > out.fasta.tmp`;
#	`mv out.fasta.tmp out.fasta`;
	`rm out`;
	}

open (READ,"<$anchPath");

my @positions;
$number = <READ>; chomp $number;
while (<READ>) {
	chomp;
	my @mas = split/\t/;
	push (@positions, [$mas[0], $mas[1],$start]);
	$history{$mas[0]} = [$mas[0],$mas[1],$start];
	}

close READ;

#open (READ,"<out.fasta");
#
#my $seq1 = <READ>;
#$seq1 = <READ>; chomp $seq1;
#my $seq2 = <READ>;
#$seq2 = <READ>; chomp $seq2;
#
#close READ;

#foreach my $arg (@positions) {
#	my $s = getseq($seq1, $seq2, ($arg)->[0], ($arg)->[1]);
#	print "",($arg)->[0],"\t",($arg)->[1],"\t",sp(($s)->[0], ($s)->[1]),"\n";
#	print "",($s)->[0],"\n";
#	print "",($s)->[1],"\n";
#	}


for (my $i = 1; $i <= $cycles; $i++) {
	print_anch([@positions]);
	launch;
#	my $shift = `perl $score $skaPath out.fasta`;
        my @shiftall = `perl $score $skaPath out.fasta | awk '{printf "%s ", \$1}'`;
#	if ($shift =~ /^(\d+.?\d*)/) {
#		$shift = $1;
#		}
	set_anchors([@positions]);
#	print "$shift\n";
	print "@shiftall\n";
	}

foreach my $key (keys %history) {
#	print STDERR join("\t", @{$history{$key}}),"\n";
	}


