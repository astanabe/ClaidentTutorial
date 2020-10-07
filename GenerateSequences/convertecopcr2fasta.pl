while (<>) {
	unless (/^#/) {
		my @entry = split(/ +\| +/);
		chomp($entry[21]);
		print('>' . $entry[0] . ' taxid=' . $entry[2] . '; ' . $entry[21] . "\n");
		print($entry[13] . $entry[20] . reversecomplement($entry[16]) . "\n");
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
