use strict;
my $indexfile = $ARGV[0];
my %index;
my $pindexhopping = 0.01;
my @inputfiles;
my %inputfiles;
my $inputfilehandle;
my $outputfilehandle;
my $seed = time^$$;
my $gen;

for (my $i = 1; $i < scalar(@ARGV); $i ++) {
	my @temp = glob($ARGV[$i]);
	if (scalar(@temp) > 0) {
		foreach (@temp) {
			if (!exists($inputfiles{$_})) {
				$inputfiles{$_} = 1;
				push(@inputfiles, $_);
			}
			else {
				die(__LINE__ . " error.\n");
			}
		}
	}
	else {
		die(__LINE__ . " error.\n");
	}
}

{
	my $samplename;
	open($inputfilehandle, "< $indexfile") or die(__LINE__ . " error.\n");
	while (<$inputfilehandle>) {
		if (/^>(\S+)/) {
			$samplename = $1;
		}
		elsif (/^(\S+)/ && $samplename) {
			$index{$samplename} = $1;
		}
		else {
			die(__LINE__ . " error.\n");
		}
	}
	close($inputfilehandle);
}

# initialize random number generator
eval "use Math::Random::MT::Auto";
if ($@) {
	eval "use Math::Random::MT";
	if ($@) {
		eval "use Math::Random::MT::Perl";
		if ($@) {
			die(__LINE__ . ": Perl module \"Math::Random::MT::Auto\", \"Math::Random::MT\" and \"Math::Random::MT:Perl\" are not available.\n");
		}
		else {
			$gen = Math::Random::MT::Perl->new($seed);
		}
	}
	else {
		$gen = Math::Random::MT->new($seed);
	}
}
else {
	$gen = Math::Random::MT::Auto->new();
	$gen->srand($seed);
}

foreach my $inputfile (@inputfiles) {
	my $samplename = $inputfile;
	$samplename =~ s/\_\d\..+$//;
	unless ($index{$samplename}) {
		die(__LINE__ . " error.\n");
	}
	my @allsamples = keys(%index);
	open($inputfilehandle, "< $inputfile") or die(__LINE__ . " error.\n");
	my $lineno = 1;
	my $seqno = 1;
	while (<$inputfilehandle>) {
		s/\r?\n?$//;
		if ($lineno % 4 == 1) {
			if (/^\@([A-Za-z0-9_\.\-]+)__([A-Za-z0-9_\.\-]+)__([A-Za-z0-9_\.\-]+)__([A-Za-z0-9_]+)\-\d+/) {
				print(STDOUT '@' . "$1:$2:$3:$4" . '_' . $samplename . '_' . sprintf("%05d", $seqno) . "\n");
			}
			else {
				die(__LINE__ . " error.\n$_\n");
			}
		}
		elsif ($lineno % 4 == 2) {
			if ($gen->rand(1) < $pindexhopping) {
				my $tempindex = $index{$allsamples[int($gen->rand(scalar(@allsamples)))]};
				print(STDOUT $tempindex . "\n");
			}
			else {
				print(STDOUT $index{$samplename} . "\n");
			}
		}
		elsif ($lineno % 4 == 3) {
			print(STDOUT '+' . "\n");
		}
		elsif ($lineno % 4 == 0) {
			foreach (1..length($index{$samplename})) {
				print(STDOUT 'I');
			}
			print(STDOUT "\n");
			$seqno ++;
		}
		$lineno ++;
	}
	close($inputfilehandle);
}
