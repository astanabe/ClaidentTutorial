use strict;

my $nsample = $ARGV[0];
my $nblank = $ARGV[1];
my $nspecies = $ARGV[2];
my $ndiff = $ARGV[3];
my $inputfile = $ARGV[-1];
my %seq;
my @seqname;
my @community;
my $seed = time^$$;
my $gen;
my $inputfilehandle;
my $outputfilehandle;

# read input file
{
	my $seqname;
	open($inputfilehandle, "< $inputfile") || die(__LINE__ . "error.\n");
	while (<$inputfilehandle>) {
		if (/^>(.+)\r?\n?/) {
			$seqname = $1;
		}
		elsif ($seqname && s/[\s\r\n]//g) {
			$seq{$seqname} .= $_;
		}
	}
	close($inputfilehandle);
	@seqname = keys(%seq);
}

# initialize random number generator
eval "use Math::Random::MT::Auto";
if ($@) {
	eval "use Math::Random::MT::Perl";
	if ($@) {
		&errorMessage(__LINE__, "Perl module \"Math::Random::MT::Auto\" and \"Math::Random::MT:Perl\" is not available.");
	}
	else {
		$gen = Math::Random::MT::Perl->new($seed);
	}
}
else {
	$gen = Math::Random::MT::Auto->new();
	$gen->srand($seed);
}

# make sample community
{
	my @temp = @seqname;
	my @previous;
	for (my $i = 0; $i < $nsample; $i ++) {
		my @tempcommunity;
		if ($i == 0) {
			for (my $j = 0; $j < $nspecies; $j ++) {
				push(@tempcommunity, splice(@temp, int($gen->rand(scalar(@temp))), 1));
			}
		}
		else {
			for (my $j = 0; $j < $ndiff; $j ++) {
				push(@tempcommunity, splice(@temp, int($gen->rand(scalar(@temp))), 1));
			}
			for (my $j = 0; $j < $nspecies - $ndiff; $j ++) {
				push(@tempcommunity, splice(@previous, int($gen->rand(scalar(@previous))), 1));
			}
		}
		push(@{$community[$i]}, @tempcommunity);
		push(@temp, @previous);
		@previous = @tempcommunity;
	}
}

# make blank sample community
{
	my %allspecies;
	for (my $i = 0; $i < scalar(@community); $i ++) {
		for (my $j = 0; $j < scalar(@{$community[$i]}); $j ++) {
			$allspecies{$community[$i][$j]} = 1;
		}
	}
	my @allspecies = keys(%allspecies);
	for (my $i = 0; $i < $nblank; $i ++) {
		my @temp = @allspecies;
		my @tempcommunity;
		for (my $j = 0; $j < $nspecies; $j ++) {
			push(@tempcommunity, splice(@temp, int($gen->rand(scalar(@temp))), 1));
		}
		push(@{$community[($nsample + $i)]}, @tempcommunity);
	}
}

# output reference files
{
	for (my $i = 0; $i < scalar(@community); $i ++) {
		if ($i < $nsample) {
			open($outputfilehandle, "> Sample" . sprintf("%02d", $i + 1) . ".fasta");
			for (my $j = 0; $j < scalar(@{$community[$i]}); $j ++) {
				print($outputfilehandle ">" . $community[$i][$j] . "\n" . $seq{$community[$i][$j]} . "\n");
			}
			close($outputfilehandle);
		}
		else {
			open($outputfilehandle, "> Blank" . sprintf("%02d", $i - $nsample + 1) . ".fasta");
			for (my $j = 0; $j < scalar(@{$community[$i]}); $j ++) {
				print($outputfilehandle ">" . $community[$i][$j] . "\n" . $seq{$community[$i][$j]} . "\n");
			}
			close($outputfilehandle);
		}
	}
}
