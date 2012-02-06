package Statistics::Multtest;

use List::Vectorize;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(bonferroni holm hommel hochberg BH BY qvalue);
our %EXPORT_TAGS = (all => [qw(bonferroni holm hommel hochberg BH BY qvalue)]);

our $VERSION = '0.10';

1;

sub initial {
	my $p = shift;
	
	unless(ref($p) eq "ARRAY" or ref($p) eq "HASH") {
		die "ERROR: P-values should be array ref or hash ref.\n";
	}
	
	my $name = [];
	if(ref($p) eq "HASH") {
		$name = [ keys %{$p} ];
		$p = [ values %{$p} ];
	}
	
	if(max($p) > 1 or min($p) < 0) {
		die "ERROR: P-values should between 0 and 1.\n";
	}
	
	return ($name, $p);
}


sub get_result {
	my ($adjp, $name) = @_;
	
	if(len($name) == 0) {
		return $adjp;
	}
	else {
		my $result;
		for (0..$#$name) {
			$result->{$name->[$_]} = $adjp->[$_];
		}
		return $result;
	}
}

sub bonferroni {
    my $p = shift;
    my $name;
    ($name, $p) = initial($p);
    my $adjp = _bonferroni($p);
	return get_result($adjp, $name);
}

sub holm {
    my $p = shift;
    my $name;
    ($name, $p) = initial($p);
    my $adjp = _holm($p);
	return get_result($adjp, $name);
}

sub hommel {
    my $p = shift;
    my $name;
    ($name, $p) = initial($p);
    my $adjp = _hommel($p);
	return get_result($adjp, $name);
}

sub hochberg {
    my $p = shift;
    my $name;
    ($name, $p) = initial($p);
    my $adjp = _hochberg($p);
	return get_result($adjp, $name);
}

sub BH {
    my $p = shift;
    my $name;
    ($name, $p) = initial($p);
    my $adjp = _BH($p);
	return get_result($adjp, $name);
}

sub BY {
    my $p = shift;
    my $name;
    ($name, $p) = initial($p);
    my $adjp = _BY($p);
	return get_result($adjp, $name);
}

sub qvalue {
    my $p = shift;
    my $name;
    ($name, $p) = initial($p);
    my $adjp = _qvalue($p);
	return get_result($adjp, $name);
}

sub _qvalue {
	my $p = shift;
	my $o = order($p);
	$p = subset($p, $o);
	my $ro = order($o);
	
	my $m = scalar(@$p);
	my @lambda;
	for(my $i = 0; $i < 19; $i++) {
		$lambda[$i] = $i*0.05;
	}
	
	#
	my $pi0;
	if(length(@lambda) == 1) {
		if($lambda[0] >0 and $lambda[0] < 1) {
			$pi0 = scalar(@{&subset($p, sub{$_[0] >= $lambda[0]})})/$m/(1-$lambda[0]);
			$pi0 = min([$pi0, 1]);
		}
	}
	else {
		# pi0[lamda]
		my @pi0;
		for(my $i = 0; $i < scalar(@lambda); $i ++) {
			$pi0[$i] = scalar(@{&subset($p, sub{$_[0] >= $lambda[0]})})/$m/(1-$lambda[$i]);
		}
		
		# bootstrap
		my $minpi0 = min(\@pi0);
		my @mse;
		my @pi0_boot;
		
		for(my $i = 0; $i < 100; $i ++) {
			my @p_boot = sample($p, $m, "replace" => 1);
			my @pi0_boot;
			for(my $j = 0; $j < scalar(@lambda); $j ++) {
				$pi0_boot[$j] = scalar(@{subset(\@p_boot, sub{$_[0] > $lambda[$j]})})/$m/(1-$lambda[$j]);
				$mse[$j] += ($pi0_boot[$j] - $minpi0)**2;
			}
		}
		
		my $minmse = min(\@mse);
		my @tmp_pi0;
		for(my $j = 0; $j < scalar(@lambda); $j ++) {
			if(abs($mse[$j]-$minmse) < 0.000001) {
				push(@tmp_pi0, $pi0[$j]);
			}
		}
		$pi0 = min(\@tmp_pi0);
		$pi0 = min([$pi0, 1]);
	}
	print "pi0 = $pi0\n";
	my $q = [];
	if($pi0 <= 0) {
		return $q;
	}
	
	$q->[$#$p] = min([$pi0*$p->[$#$p], 1]);     # q(p[m]) = pi0 * p[m]
	for(my $i = $#$p - 1; $i >= 0; $i --) {
		$q->[$i] = min([$pi0*$m*$p->[$i]/($i+1), $q->[$i+1]]);
		$q->[$i] = min([$q->[$i], 1]);
	}
	
	return subset($q, $ro);
}

# R code: pmin(1, n * p)
sub _bonferroni {
    my $p = shift;
    my $n = len($p);
    
    my $adjp = multiply($n, $p);
    
    return pmin(1, $adjp);
}

# R code: i = 1:n
#         o = order(p)
#         ro = order(o)
#         pmin(1, cummax((n - i + 1) * p[o]))[ro]
sub _holm {
    my $p = shift;
    my $n = len($p);
	
	my $i = seq(1, $n);
	my $o = order($p);
	my $ro = order($o);
    
    my $adjp = multiply(minus($n + 1, $i), subset($p, $o));
    $adjp = cumf($adjp, \&max);
	$adjp = pmin(1, $adjp);
    
    return subset($adjp, $ro);
}

# R code: i = 1:n
#         o = order(p)
#         p = p[o]
#         ro = order[o]
#         q = pa = rep.int(min(n * p/i), n)
#         for (j in (n - 1):2) {
#             ij = 1:(n - j + 1)
#             i2 = (n - j + 2):n
#             q1 <- min(j * p[i2]/(2:j))
#             q[ij] <- pmin(j * p[ij], q1)
#             q[i2] <- q[n - j + 1]
#             pa <- pmax(pa, q)
#         }
#         pmax(pa, p)[ro]
sub _hommel {
    my $p = shift;
    my $n = len($p);
	
	my $i = seq(1, $n);
	my $o = order($p);
	$p = subset($p, $o);
	my $ro = order($o);
	
	my $pa = rep(min(divide(multiply($n, $p), $i)), $n);
	my $q = copy($pa);
	
	# set the first index as 1
    unshift(@$p, 0);
    unshift(@$q, 0);
    unshift(@$pa, 0);
	
	my $ij;
	my $i2;
	my $q1;
    for my $j (@{seq($n - 1, 2)}) {
        
		$ij = seq(1, $n - $j + 1);
		$i2 = seq($n - $j + 2, $n);
		$q1 = min(divide(multiply($j, subset($p, $i2)), seq(2, $j)));
		subset_value($q, $ij, pmin(multiply($j, subset($p, $ij)), $q1));
		subset_value($q, $i2, $q->[$n - $j + 1]);
        $pa = pmax($pa, $q);
    }
    
    shift(@$p);
    shift(@$q);
    shift(@$pa);
	
	my $adjp = pmax($pa, $p);
	return subset($adjp, $ro);    
}

# R code: i = n:1
#         o <- order(p, decreasing = TRUE)
#         ro <- order(o)
#         pmin(1, cummin((n - i + 1) * p[o]))[ro]
sub _hochberg {
    
    my $p = shift;
    my $n = len($p);
    my $i = seq($n, 1);
    
    my $o = order($p, sub {$_[1] <=> $_[0]});
	my $ro = order($o);
	
    my $adjp = multiply(minus($n+1, $i), subset($p, $o));
    $adjp = cumf($adjp, \&min);
	$adjp = pmin(1, $adjp);
    return subset($adjp, $ro);
}

# R code: i <- n:1
#         o <- order(p, decreasing = TRUE)
#         ro <- order(o)
#         pmin(1, cummin(n/i * p[o]))[ro]
sub _BH {
    my $p = shift;
    my $n = len($p);
    my $i = seq($n, 1);
    
    my $o = order($p, sub {$_[1] <=> $_[0]});
	my $ro = order($o);
	
    my $adjp = multiply(divide($n, $i), subset($p, $o));
    $adjp = cumf($adjp, \&min);
	$adjp = pmin(1, $adjp);
    return subset($adjp, $ro);

}

# R code: i <- n:1
#         o <- order(p, decreasing = TRUE)
#         ro <- order(o)
#         q <- sum(1/(1L:n))
#         pmin(1, cummin(q * n/i * p[o]))[ro]
sub _BY {
    
    my $p = shift;
    my $n = len($p);
    my $i = seq($n, 1);
    
    my $o = order($p, sub {$_[1] <=> $_[0]});
	my $ro = order($o);
	
    my $q = sum(divide(1, seq(1, $n)));
    my $adjp = multiply(divide($q*$n, $i), subset($p, $o));
    $adjp = cumf($adjp, \&min);
    $adjp = pmin(1, $adjp);
    return subset($adjp, $ro);
}

sub pmin {
	my $array1 = shift;
	my $array2 = shift;
	
	return mapply($array1, $array2, sub {min(\@_)});
}

sub pmax {
	my $array1 = shift;
	my $array2 = shift;
	
	return mapply($array1, $array2, sub {max(\@_)});
}

__END__

=pod

=head1 NAME

Statistics::Multtest - Control false discovery rate in multiple test problem

=head1 SYNOPSIS

  use Statistics::Multtest qw(bonferroni holm hommel hochberg BH BY qvalue);
  use strict;
  
  my $p;
  # p-values can be stored in an array by reference
  $p = [0.01, 0.02, 0.05,0.41,0.16,0.51];
  # @$res has the same order as @$p
  my $res = BH($p);
  print join "\n", @$res;
  
  # p-values can also be stored in a hash by reference
  $p = {"a" => 0.01,
        "b" => 0.02,
        "c" => 0.05,
        "d" => 0.41,
        "e" => 0.16,
        "f" => 0.51 };
  # $res is also a hash reference which is the same as $p
  $res = qvalue($p);
  foreach (sort {$res->{a} <=> $res->{$b}} keys %$res) {
      print "$_ => $res->{$_}\n";
  }

=head1 DESCRIPTION

For statistical test, p-value is the probability of false positives. While there
are many hypothesis testing simultaneously, the probability of getting at least one
false positive would be large. Therefore the origin p-values should be adjusted to decrease
the false discovery rate.

Seven procedures to controlling false positive rates is provided. 
The name of the methods are derived from C<p.adjust> in 
C<stat> package and C<qvalue> in C<qvalue> package in R.
Code is translated directly from R to Perl using L<List::Vectorize> module.

All seven subroutine receive one argument which can either be an array reference
or a hash reference.

=head2 Subroutines

=over 4

=item C<bonferroni($pvalue)>

Bonferroni single-step process.

=item C<hommel($pvalue)>

Hommel singlewise process.

Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika, 75, 383每386. 

=item C<holm($pvalue)>

Holm step-down process.

Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6, 65每70. 

=item C<hochberg($pvalue)>

Hochberg step-up process.

Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika, 75, 800每803. 

=item C<BH($pvalue)>

Benjamini and Hochberg, controlling the FDR.

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289每300. 

=item C<BY($pvalue)>

Use Benjamini and Yekutieli.

Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165每1188. 

=item C<qvalue($pvalue)>

Storey and Tibshirani.

Storey JD and Tibshirani R. (2003) Statistical significance for genome-wide experiments. Proceedings of the National Academy of Sciences, 100: 9440-9445. 

=back

=head1 AUTHOR

Zuguang Gu E<lt>jokergoo@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2012 by Zuguang Gu

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.1 or,
at your option, any later version of Perl 5 you may have available.

=cut
