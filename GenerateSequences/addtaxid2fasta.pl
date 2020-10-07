use strict;
use LWP::UserAgent;

while (<>) {
	if (/^>.+__([A-Za-z0-9_]+)/) {
		my $acc = $1;
		print(STDERR "Now downloading taxonomy information of $acc...\n");
		my $taxid = &gettaxid($acc);
		s/^>(.+__[A-Za-z0-9_]+)/>$acc taxid=$taxid; $1/;
	}
	print(STDOUT $_);
}

sub gettaxid {
	my $acc = shift(@_);
	my $taxid;
	my $ua = LWP::UserAgent->new;
	$ua->timeout(30);
	$ua->agent('addlabel2fasta');
	$ua->env_proxy;
	my $baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id=' . &encodeURL($acc);
	my $req = HTTP::Request->new(POST => $baseurl . '&retmode=xml');
	my $res = $ua->request($req);
	if ($res->is_success) {
		foreach (split(/\n/,$res->content)) {
			if (/<Item Name="TaxId" Type="Integer">(\d+)<\/Item>/i) {
				$taxid = $1;
				last;
			}
		}
		sleep(5);
	}
	else {
		sleep(30);
		my $req = HTTP::Request->new(POST => $baseurl . '&retmode=xml');
		my $res = $ua->request($req);
		if ($res->is_success) {
			foreach (split(/\n/,$res->content)) {
				if (/<Item Name="TaxId" Type="Integer">(\d+)<\/Item>/i) {
					$taxid = $1;
					last;
				}
			}
		}
		else {
			die(__LINE__ . ": $acc error.\n");
		}
	}
	if ($taxid) {
		return($taxid);
	}
	else {
		sleep(30);
		my $req = HTTP::Request->new(POST => $baseurl . '&retmode=xml');
		my $res = $ua->request($req);
		if ($res->is_success) {
			foreach (split(/\n/,$res->content)) {
				if (/<Item Name="TaxId" Type="Integer">(\d+)<\/Item>/i) {
					$taxid = $1;
					last;
				}
			}
		}
		else {
			die(__LINE__ . ": $acc error.\n");
		}
		if ($taxid) {
			return($taxid);
		}
		else {
			die(__LINE__ . ": $acc error.\n");
		}
	}
}

sub encodeURL {
	my $str = shift;
	$str =~ s/([^\w ])/'%'.unpack('H2', $1)/eg;
	$str =~ tr/ /+/;
	return $str;
}
