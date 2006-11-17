# A simple execution time benchmark script. Thanks to Vasilis Christaras.
use Benchmark;
chomp($ARGV[0]);
$start= new Benchmark;
system($ARGV[0]);
$end=new Benchmark;
$df=timediff($end,$start);
print "The time: ",timestr($df,all)."\n";

