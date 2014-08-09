#!/usr/bin/Perl

use strict;
use threads;
use Benchmark qw(:hireswallclock);

my $starttime = Benchmark->new;
my $finishtime;
my $timespent;
my $num_of_threads = 2;

my @threads = initThreads();
foreach(@threads){
		$_ = threads->create(\&doOperation);
	}
foreach(@threads){
	$_->join();
}
$finishtime = Benchmark->new;
$timespent = timediff($finishtime,$starttime);
print "\nDone!\nSpent ". timestr($timespent);

print "\n\nNow trying without threading:\n\n";

my $starttime = Benchmark->new;

doWithoutThread();

$finishtime = Benchmark->new;
$timespent = timediff($finishtime,$starttime);
select(STDOUT);
print "\nDone!\nSpent ". timestr($timespent);

print "\nProgram Done!\nPress Enter to exit";
$a = <>;

sub initThreads{
	my @initThreads;
	for(my $i = 1;$i<=$num_of_threads;$i++){
		push(@initThreads,$i);
	}
	return @initThreads;
}
sub doOperation{
	# Get the thread id. Allows each thread to be identified.
	my $id = threads->tid();
	my $i = 0;
	while($i < 100000000){
			$i++
	}
	print "Thread $id done!\n";
	# Exit the thread
	threads->exit();
}
sub doWithoutThread{
	my $c = 0;
	for(my $i=0;$i<$num_of_threads;$i++){
		while($c < 100000000){
			$c++;
		}
		$c=0;
		print "Count $i done!\n";
	}
}

