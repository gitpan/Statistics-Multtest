use strict;
use warnings;
use Test::More tests => 3;

BEGIN { use_ok('Statistics::Multtest') }

my $p1 = [0.1041, 0.0012, 0.2014, 0.5132, 0.0024, 0.1451, 0.0001, 0.1541];
my $p2= {"a1" => 0.1041,
         "a2" => 0.0012,
		 "a3" => 0.2014,
		 "a4" => 0.5132,
		 "a5" => 0.0024,
		 "a6" => 0.1451,
		 "a7" => 0.0001,
		 "a8" => 0.1541,};

my $adjp1 = Statistics::Multtest::bonferroni($p1);
my $adjp2 = Statistics::Multtest::bonferroni($p2);

is( ($adjp1->[1] - 0.0130 > 0.001) + 0, 0);
is( ($adjp2->{"a4"} - 1.0000 > 0.001) + 0, 0);
