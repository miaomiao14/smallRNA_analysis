#!/usr/bin/perl
if(@ARGV<2) {print "usage: $0 inputmatchout outputnormbed\n"; exit;}

$inputmatchout = shift @ARGV;
$outputnormbed = shift @ARGV;

$tempbed="$outputnormbed.temp";
for(0..10) {$tempbed.=(0..9)[rand(10)];}

`match2all2bed $inputmatchout > $tempbed`;
`count_ntm_from_match2allbed.pl $tempbed > $outputnormbed`;
`rm -f $tempbed`;