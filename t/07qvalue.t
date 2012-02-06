use strict;
use warnings;
use Test::More tests => 2;

BEGIN { use_ok('Statistics::Multtest') }

my $p1 = [0.80234001, 0.08046698, 0.22305289, 0.44636930, 0.01076263, 0.97254271];

my $adjp1 = Statistics::Multtest::qvalue($p1);

is( ($adjp1->[1] - 0.2414 > 0.001) + 0, 0);
