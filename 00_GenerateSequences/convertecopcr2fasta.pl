foreach (@ARGV) {
	if (/^-+(?:fp|forwardprimer)=([ACGT]+)/i) {
		$forwardprimer = uc($1);
	}
	elsif (/^-+(?:rp|reverseprimer)=([ACGT]+)/i) {
		$reverseprimer = uc($1);
	}
}

while (<STDIN>) {
	unless (/^#/) {
		my @entry = split(/ +\| +/);
		chomp($entry[21]);
		print('>' . $entry[0] . ' taxid=' . $entry[2] . '; ' . $entry[21] . "\n");
		if ($forwardprimer && $reverseprimer) {
			print($forwardprimer . $entry[20] . reversecomplement($reverseprimer) . "\n");
		}
		else {
			print($entry[13] . $entry[20] . reversecomplement($entry[16]) . "\n");
		}
	}
}

sub reversecomplement {
	my @temp = split('', $_[0]);
	my @seq;
	foreach my $seq (reverse(@temp)) {
		$seq =~ tr/ACGTMRYKVHDBacgtmrykvhdb/TGCAKYRMBDHVtgcakyrmbdhv/;
		push(@seq, $seq);
	}
	return(join('', @seq));
}
