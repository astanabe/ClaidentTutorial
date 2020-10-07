use strict;
my $indexfile = $ARGV[0];
my %index;
my @inputfiles;
my %inputfiles;
my $inputfilehandle;
my $outputfilehandle;

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

foreach my $inputfile (@inputfiles) {
	my $samplename = $inputfile;
	$samplename =~ s/\_\d\..+$//;
	unless ($index{$samplename}) {
		die(__LINE__ . " error.\n");
	}
	open($inputfilehandle, "< $inputfile") or die(__LINE__ . " error.\n");
	my $lineno = 1;
	my $seqno = 1;
	while (<$inputfilehandle>) {
		s/\r?\n?$//;
		if ($lineno % 4 == 1) {
			print(STDOUT '@' . $samplename . '_' . sprintf("%05d", $seqno) . "\n");
		}
		elsif ($lineno % 4 == 2) {
			print(STDOUT $index{$samplename} . "\n");
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
