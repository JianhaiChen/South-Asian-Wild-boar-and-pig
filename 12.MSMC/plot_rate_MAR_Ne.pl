#!/usr/bin/perl -w
# plot MSMC result 
# effective population size: `cat CNN.final.txt CNS.final.txt >merge`
# ssh node106 perl plot_rate_MAR.pl -M CNN,CNS -u 0.75e-08 -g 5 -w 5 -R Ne merge
# Author: lh3

use strict;
use warnings;
use Getopt::Std;

my $version = "0.2.0";

my %opts = (u=>2.5e-8, 's'=>100, Y=>0, m=>5, X=>0, M=>'',
            x=>100, y =>1000,n=>20, g=>25, f=>"Helvetica,22",w=>4, P=>"right top", T=>'');
getopts('x:y:u:s:X:Y:RGpm:n:M:N:g:f:w:P:T:', \%opts);
die("
        Usage:   msmc_plot.pl [options] <out.prefix> <in.msmc(merge)>\n
        Options: 
        -u FLOAT   absolute mutation rate per nucleotide [$opts{u}]
        -s INT     skip used in data preparation [$opts{s}]
        -X FLOAT   maximum generations, 0 for auto [$opts{X}]
        -x FLOAT   minimum generations, 0 for auto [$opts{x}]
        -y FLOAT   minimum popsize, 0 for auto [$opts{y}]
        -Y FLOAT   maximum popsize, 0 for auto [$opts{Y}]
        -m INT     minimum number of iteration [$opts{m}]
        -n INT     take n-th iteration (suppress GOF) [$opts{n}]
        -M titles  multiline mode [null]
        -f STR     font for title, labels and tics [$opts{f}]
        -g INT     number of years per generation [$opts{g}]
        -w INT     line width [$opts{w}]
        -P STR     position of the keys [$opts{P}]
        -T STR     figure title [null]
        -p         convert to PDF (with epstopdf)
        -R         do not remove temporary files
        -G         plot grid
        \n") if (@ARGV < 2);

my $prefix = shift(@ARGV);
my $result = shift(@ARGV);
my (@data, $d, $N0, $skip, $Mseg, $Msize, $id, $min_ri, $do_store, $gof, $round, @FN, @nscale, @tscale, @alpha, $dt);
# load data
my ($max_seg, $max_size, $fh);
$max_seg = $max_size = 0;
open IN,$result or die $!;
my $i='N';
while(<IN>){
    chomp;
    my @a=split;
    my $line=$_;
    my ($index,$left_time,$right_time,$lambda)=(@a)[0,1,2,3];
    if($line=~/time_index/ && $i eq 'N') {
        $i=0;
        if( -e "$prefix.$i.txt") {`rm $prefix.$i.txt`;}
        open TEM,'>>',"$prefix.$i.txt" || die;
        next;
    }
    elsif($line=~/time_index/){
        $i++;
        close TEM;
        if( -e "$prefix.$i.txt") {`rm $prefix.$i.txt`;}
        open TEM,'>>',"$prefix.$i.txt" || die;
        next;
    }
    else{;}

    my $real_time = $left_time/$opts{u}*$opts{g};    
    my $effective_pop_size = (1/$lambda)/($opts{u}*2);
    print TEM "$real_time\t$effective_pop_size\n";
}
close TEM;
close IN;
my $i_sum=$i;
# plot

my $yran = ($opts{Y} > 0)? $opts{Y} : '*';
my $xran = ($opts{X} > 0)? $opts{X} : '*';
my $keyconf = $opts{M}? "set key $opts{P}" : "set key off";
my $grid = $opts{G}? "set grid" : 'unset grid';
my $afont = qq/font "$opts{f}"/;
my $lw = qq/lw $opts{w}/;
my $ylab_aux = sprintf("%.2fx10^{-8}", $opts{u}/1e-8);

open($fh, "| tee $prefix.gp | gnuplot") || die;
print $fh qq(
 set size 1, 0.8;
 set xran [$opts{x}:$xran];
 set log x;
 set format x "10^{\%L}";
 set mxtics 10;
 set ytics 10;
 set mytics 10;
 $grid;
 $keyconf;
 set xtics axis $afont;
 set ytics nomirror $afont;

 set xtics nomirror
 unset x2tics

 set xlab "Years (g=$opts{g}, {/Symbol m}=$ylab_aux)" $afont;
 set t po eps enhance so co "Helvetica,14";
 set key reverse Left samplen 1

);

print $fh qq/set title "$opts{T}";/ if ($opts{T});

print $fh qq(
  set yran [$opts{y}:$yran];
  set log y;
  set format y "10^{\%L}";
  set ylab "Effective population size" $afont;

  set style line 1 lt 1 lc rgb "#B8860B" $lw;
  set style line 2 lt 2 lc rgb "#458B74" $lw;
  set style line 3 lt 3 lc rgb "#0000FF" $lw;
  set style line 4 lt 4 lc rgb "#00008B" $lw;
  set style line 5 lt 5 lc rgb "#8A2BE2" $lw;
  set style line 6 lt 6 lc rgb "#A52A2A" $lw;
  set style line 7 lt 7 lc rgb "#FF4040" $lw;
  set style line 8 lt 8 lc rgb "#8B7355" $lw;
  set style line 9 lt 9 lc rgb "#7FFF00" $lw;
  set style line 10 lt 10 lc rgb "#458B00" $lw;
  set style line 11 lt 11 lc rgb "#D2691E" $lw;
  set style line 12 lt 12 lc rgb "#FF7256" $lw;
  set style line 13 lt 13 lc rgb "#B0E0E6" $lw;
  set style line 14 lt 14 lc rgb "yellow" $lw;
  set style line 15 lt 15 lc rgb "yellow4" $lw;
  set style line 16 lt 16 lc rgb "#8B0000" $lw;
  set style line 100 lt 100 lc rgb "grey" $lw;
  set output "$prefix.epss" ;
  plot );

my $zwj_i=2;  
  if ( $opts{M}) {
          my @titles = split(/\,/, $opts{M});
          print $fh qq("$prefix.0.txt" u 1:2 t "$titles[0]" w st ls 1);
          foreach my $i (1 .. $i_sum) {
                  print $fh qq(,"$prefix.$i.txt" u 1:2 t "$titles[$i]" w st ls $i + 1);
          }
  } else {
      print $fh qq("$prefix.0.txt" u 1:2 t "popsize" w st ls 1);
          foreach my $i (1 .. $i_sum-1) { 
              print $fh qq(, "$prefix.$i.txt" u 1:2 t "MAR" w  st ls $zwj_i); 
              $zwj_i++;

          }
 }
 
  
  close($fh);

  if (defined $opts{p}) {
          system("epstopdf $prefix.eps");
  }

    system ("convert $prefix.epss $prefix.epss.pdf  ") ;
# remove files
  unless (defined($opts{R})) {
          unlink <$prefix.*.txt>; unlink "$prefix.gp"; 
          unlink "$prefix.epss" ;
  }
if ($opts{M}) 
	  {
            unlink <*.par>;
          }
