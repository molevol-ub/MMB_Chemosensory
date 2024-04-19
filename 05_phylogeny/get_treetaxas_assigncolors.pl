#!usr/bin/perl
use strict;
use warnings;


my ($line, $name);

open (ResultsBranch, ">", "$ARGV[0]_colors_styles_branches.txt");
open (ResultsLabel, ">", "$ARGV[0]_colors_styles_labels.txt");
#open (ResultsSom, ">", "$ARGV[0]_colors_styles_species_shadow.txt");

print ResultsBranch "TREE_COLORS
SEPARATOR TAB
DATA\n";

print ResultsLabel "TREE_COLORS
SEPARATOR TAB
DATA\n"; 

#print ResultsSom "TREE_COLORS
#SEPARATOR TAB
#DATA\n"; 

open (File, "<", $ARGV[0]);
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\(|,|\)/, $line);
	foreach my $field (@subline){
		if ($field =~ /(Dmel\S+):/){ # Green
			print ResultsBranch "$1\tbranch\t#339966\tnormal\n";
			print ResultsLabel "$1\tlabel\t#339966\n";
		}
		elsif ($field =~ /(Amel\S+):/){ # Dark red
			print ResultsBranch "$1\tbranch\t#b01414\tnormal\n"; 
			print ResultsLabel "$1\tlabel\t#b01414\n";
		}
		elsif ($field =~ /(Sinv\S+):/){ # Dark blue
			print ResultsBranch "$1\tbranch\t#0043ff\tnormal\n";
			print ResultsLabel "$1\tlabel\t#0043ff\n";
		}

		# example for other species names
		elsif ($field =~ /(Dpul\w+):/){ # Light blue
			print ResultsBranch "$1\tbranch\t#0dfff0\tnormal\n"; 
			print ResultsLabel "$1\tlabel\t#0dfff0\n";
		}
		elsif ($field =~ /(Smar\w+):/){ # Orange
			print ResultsBranch "$1\tbranch\t#ff6600\tnormal\n";
			print ResultsLabel "$1\tlabel\t#ff6600\n";
		}
		elsif ($field =~ /(Isca\w+):/){ # Purple
			print ResultsBranch "$1\tbranch\t#990000\tnormal\n";
			print ResultsLabel "$1\tlabel\t#990000\n";
		}
		elsif ($field =~ /(Erow\w+):/){ # Dark blue 0043ff
			print ResultsBranch "$1\tbranch\t#0043ff\tnormal\n";
			print ResultsLabel "$1\tlabel\t#0043ff\n";
			#print ResultsSom "$1\trange\t#cccccc\n";			
		}

		elsif ($field =~ /(Hduj\w+):/){ # Dark blue 0043ff
			print ResultsBranch "$1\tbranch\t#FF66FF\tnormal\n";
			print ResultsLabel "$1\tlabel\t#FF66FF\n";
		}
		elsif ($field =~ /(Rvar\w+):/){ # Dark blue 0043ff
			print ResultsBranch "$1\tbranch\t#9900CC\tnormal\n";
			print ResultsLabel "$1\tlabel\t#9900CC\n";
		}		
		
	}
}
close File;
#close Results;

