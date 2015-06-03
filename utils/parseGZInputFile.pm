

##parse bed or bedlike format
##use readMapRef as a reference for a hash
##CLONE is SRA or DEG(determines size selection)


sub parseInputFile
{
  		my ($readMapRef, $fileHandle,$format,$CLONE) = @_;
		


        	if($format eq "normbed")
        	{
				while($fileHandle->gzreadline(my $line) > 0)
				{ 	 
					next if ($line=~/data/);
					chomp $line; 
					$line=~s/\s+/\t/g; #change all spaces to tab
					
					my $chr;
					my $bedstart;
					my $bedend;
					my $strand;
					my $seq;
					my $reads;
					my $ntm;
					my $len;
        		
	        		($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm)= split(/\t/,$line);
					$len=$bedend-$bedstart+1;
					
					if($CLONE eq "SRA") #length range to filter for piRNAs; for DEG, no filtering should be processed
	        		{
						next if ($len>29 || $len<23);
	        		}
					if ($strand eq '+')
					{
		            	$readMapRef->{$strand}{$chr}{$bedstart}+=$reads/$ntm; 
		            	#$plus_end{$chr}{$bedend}=$l[2];
		        	}
		        	else
		        	{
		            	$readMapRef->{$strand}{$chr}{$bedend}+=$reads/$ntm;
		        	}
				}#while
        	}
			if($format eq "bed")
			{
				while($fileHandle->gzreadline(my $line) > 0)
				{ 
					next if ($line=~/data/);
					chomp $line; 
					$line=~s/\s+/\t/g; #change all spaces to tab
					
					my $chr;
					my $bedstart;
					my $bedend;
					my $strand;
					my $seq;
					my $reads;
					my $ntm;
					my $len;
					my $ntmreads;
					
					($chr,$bedstart,$bedend,$seq,$ntmreads,$strand)= split(/\t/,$line);
					$reads=$ntmreads;
					$ntm=1;
					$len=$bedend-$bedstart;
					
					if($CLONE eq "SRA") #length range to filter for piRNAs; for DEG, no filtering should be processed
	        		{
						next if ($len>29 || $len<23);
	        		}
	        		#($reads,$ntm,$dep)=split(/,/,$l[3]);
	        		if ($strand eq '+') #strand information is not included in the data, but in the file name
	        		{
		            	$readMapRef->{$strand}{$chr}{$bedstart}+=$reads/$ntm;
		            	#$plus_end{$l[0]}{$l[1]}=$_[2];
		        	}
		        	else
		        	{
		            	$readMapRef->{$strand}{$chr}{$bedend}+=$reads/$ntm;
		        	}
				}#while
        	}
        	if($format eq "bed2")#bo's definition
			{
				
				while($fileHandle->gzreadline(my $line) > 0)
				{  
					next if ($line=~/data/);
					chomp $line; 
					$line=~s/\s+/\t/g; #change all spaces to tab
					
					my $chr;
					my $bedstart;
					my $bedend;
					my $strand;
					my $seq;
					my $reads;
					my $ntm;
					my $len;
					($chr,$bedstart,$bedend,$reads,$ntm,$strand,$seq)= split(/\t/,$line);
					$len=$bedend-$bedstart;
	
					if($CLONE eq "SRA")
					{
						next if ($len>29 || $len<23);
					}
	        		if($strand eq "+") #strand information is not included in the data, but in the file name
	        		{
		            	$readMapRef->{$strand}{$chr}{$bedstart}+=$reads/$ntm; 
		            	#$plus_end{$l[0]}{$l[1]}=$l[2];
		        	}
		        	else
		        	{
		            	$readMapRef->{$strand}{$chr}{$bedend}+=$reads/$ntm;
		        	}
        	}#while		
        			
		}#if
} #sub
		


1;