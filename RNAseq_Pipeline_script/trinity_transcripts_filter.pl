#!/usr/bin/perl

use File::Basename;
$psl=$ARGV[0];
$rsemfasta=$ARGV[1];
$rrna=$ARGV[2];
$f=$ARGV[3];
if(scalar(@ARGV)==0)
{
	usage();
}

open RRNA, "$rrna" or die "can not find file $rrna: $!";

while(my $line=<RRNA>)
{
	chomp $line;
	$rrna{$line}=1;
}
close(RRNA);

`grep ">" $rsemfasta | gawk -F" " '{gsub(/>/,"",\$1);print \$1}' >fpkm1above.transcripts.list`;
open RSEM, "fpkm1above.transcripts.list" or die "can not find file fpkm1above.transcripts.list: $!";
while(my $line=<RSEM>)
{
	chomp $line;	
	$rsem{$line}=1;
}
close(RSEM);

open PSL, "$psl" or die "can not find file $psl: $!";
open PSLOUT, ">$psl.xrRNA2.xlowq.$f.psl" or die "can not write to file $psl.xrRNA2.xlowq.$f.psl: $!";
while (my $line=<PSL>)
{
	chomp $line;
	@fileds=split(/\t/,$line);
	my $name=$fileds[9];
	$psl{$name}=1;
	if( (!exists $rrna{$name}) && $rsem{$name})
#	if ($rsem{$name})
	{
		print PSLOUT "$line\n";
	}
}
close(PSLOUT);
close(PSL);

`fasta_formatter -i $rsemfasta -o $rsemfasta.tab -t`;
open FASTA, "$rsemfasta.tab" or die "can not find $rsemfasta.tab: $!";
open FFASTA, ">$rsemfasta.xrRNA2.refmap" or die "can not write to file $rsemfasta.xrRNA2.refmap: $!";
while (my $line=<FASTA>)
{
	chomp $line;
	@f=split(/\s/,$line);
	my $name=$f[0];
	#print "$name\t$f[1]\n";
	if( $psl{$name} && (!exists $rrna{$name}) )
	{
		print FFASTA "$line\n"; 
	}
}

close(FFASTA);
close(FASTA);
`awk 'BEGIN{FS="\\t";OFS="\\n";}{print ">"\$1,\$2}' $rsemfasta.xrRNA2.refmap >$rsemfasta.xrRNA2.refmap.fasta`;
#`fasta_formatter -i $rsemfasta.xrRNA2.refmap.temp -o $rsemfasta.xrRNA2.refmap.fasta -w 60`;

#foreach $rrnac (keys %rsem)
#{
#	print "$rrnac\n";
#}
`rm fpkm1above.transcripts.list`;
`rm $rsemfasta.tab`;
`rm $rsemfasta.xrRNA2.refmap`;
#`rm $rsemfasta.xrRNA2.refmap.temp`;
sub usage
{
        print "\nUsage:$0\n\n\t";
        print "REQUIRED\n\t";
	print "PSL RSEM_fasta rRNA_list rpkm\n";
	
	exit(1);
}

