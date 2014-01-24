#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
#BEGIN { unshift @INC,"/home/wangw1/bin/";}
#require "sort_hash_key.pm";
use File::Basename;
use Compress::Zlib;

# the distribution of 3'-5' end distance of piRNAs from the same strand.
#09/14/2012 fixed a bug...for those have multiple minimal distance, report all not just first two.
#separate length. For each length of piRNAs, check the existence of the peak at -22.
#WEI WANG
#09/17/2012 fixed another leaky part...


#10/17/2012
#find the nearest 5' end for each 3' end, no overlap between 3'end of last piRNA and 5'end of the next piRNA
#length specific


##12/11/2012
## extend the analysis to the upstream of each 3?end
   $file1=fileparse($ARGV[0]);  

   
   
   
   %piRNA=();
   #%piRNA_read=();
    #my $gz = gzopen($ARGV[$i], "rb") or die "Cannot open $ARGV[0]: $gzerrno\n" ;
    #while($gz->gzreadline($_) > 0)
   open IN, $ARGV[0];
    while(<IN>)
    {
        chomp; split(/\t/);  
        next if (length($_[4])>29 || length($_[4])<23);
        next if (/data/);
        $len=length($_[4]);
        $r=$_[5]/$_[6];  
        if ($_[3] eq '+')
        {
            $piRNA{$_[0]}{plus}{$_[1]}{$_[2]}+=$r;   #start, end, reads#
            #$piRNA_read{$_[0]}{plus}{$len}{$_[1]}{$_[2]}+=$r;
         
        }
        else
        {
            $piRNA{$_[0]}{minus}{$_[1]}{$_[2]}+=$r;  #start, end, reads# for minus strand, the start is the 3'end
            #$piRNA_read{$_[0]}{minus}{$len}{$_[1]}{$_[2]}+=$r;
            
        }
    }

    open OUT, ">$file1.distance.product.distribution";
    open O, ">$file1.distance.min.distribution";
    %h_no_dis_product=();
    %h_no_dis=();
    foreach $chr (keys %piRNA)
    {
        foreach $strand (keys %{$piRNA{$chr}})
        {
            %plus_chr = %{$piRNA{$chr}{$strand}};
            @plus_sort = &sort_double_key_hash( %plus_chr ); ##sort by the numerical value of the key
                
            %dis_value_no=();
            %dis_value_no_product=();
                
                
            foreach  ($k=0;$k<$#plus_sort;$k++)
            {
                foreach ($j=-10;$j<30;$j++)
                {
                    $next_start=$plus_sort[$k]->[1]+$j;  ##the start of next piRNA species
                    if(%{$piRNA{$chr}{$strand}{$next_start}})
                    {
                            foreach $pair_end (keys %{$piRNA{$chr}{$strand}{$next_start}})
                            {
                                $dis_value_no{$j}+=&min($plus_sort[$k]->[2],$piRNA{$chr}{$strand}{$next_start}{$pair_end});
                                $dis_value_no_product{$j}+=$plus_sort[$k]->[2]*$piRNA{$chr}{$strand}{$next_start}{$pair_end};
                                $h_no_dis{$j}+=&min($plus_sort[$k]->[2],$piRNA{$chr}{$strand}{$next_start}{$pair_end});
                                $h_no_dis_product{$j}+=$plus_sort[$k]->[2]*$piRNA{$chr}{$strand}{$next_start}{$pair_end};
                                
                            }
                        last;   
                    }
                    
                }    
                            
            }
            foreach $d (sort {$a <=> $b } keys %dis_value_no)
            {
                print O "$chr\t$strand\t$d\t$dis_value_no{$d}\n";
            }
            foreach $d (sort {$a <=> $b } keys %dis_value_no_product)
            {
                print OUT "$chr\t$strand\t$d\t$dis_value_no_product{$d}\n";
            }
            
            
            
            
        }#strand
    }#chr
    close(O);
    close(OUT);
 
    
    
    
    open OUT, ">$file1.distance.min.distribution.summary";

        foreach $dis (sort {$a <=> $b } keys  %h_no_dis)
        {
        print OUT "$dis\t$h_no_dis{$dis}\n";
        }
 
    close(OUT);
    open OUT, ">$file1.distance.product.distribution.summary";

    foreach $dis (sort {$a <=> $b } keys  %h_no_dis_product)
    {
        print OUT "$dis\t$h_no_dis_product{$dis}\n";
    }
 
    close(OUT);

sub sort_double_key_hash
{
my %hash= @_;
my @sortarray=();
foreach $start ( sort {$a <=> $b } keys %hash)
{
    foreach $end ( sort {$a <=> $b} keys %{$hash{$start}})
    {
        #print "$start\t$end\t$hash{$start}{$end}\n";
        push @sortarray, [$start,$end,$hash{$start}{$end}];
    }
}
return @sortarray;
}


sub sort_arrayref_value
{
my @array= @_;
my @sortarray=();
foreach $k ( sort {$a->[0] <=> $b->[0] } @array)
{
#       print "$k->[0],$k->[1],$k->[2]\n";
        push @sortarray, [$k->[0],$k->[1],$k->[2]];
}
return @sortarray;
}