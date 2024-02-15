my $lineno = 1;
while (<>) {
	if ($lineno % 4 == 1 && /^\@/) {
		print('@M99999:999:999999999-Z9ZZZ:1:9999:' . ($lineno + 10000) . ':' . ($lineno + 2000) . ' 1:N:0:0' . "\n");
	}
	else {
		print($_);
	}
	$lineno ++;
}
