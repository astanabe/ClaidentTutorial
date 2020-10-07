use strict;
use LWP::UserAgent;

while (<>) {
	if (/^>.+__([A-Za-z0-9_]+)/) {
		my $acc = $1;
		my $taxid = 7898;
		s/^>(.+__[A-Za-z0-9_]+)/>$acc taxid=$taxid; $1/;
	}
	print(STDOUT $_);
}
