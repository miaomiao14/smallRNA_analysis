BEGIN {
    if (!p) { print "p is not set. please pass in p by the -v option"; exit 1;}
    p_value_df1[0.2500] = 1.3200; p_value_df2[0.2500] = 2.7700 
    p_value_df1[0.2000] = 1.6400; p_value_df2[0.2000] = 3.2200 
    p_value_df1[0.1500] = 2.0700; p_value_df2[0.1500] = 3.7900 
    p_value_df1[0.1000] = 2.7100; p_value_df2[0.1000] = 4.6100 
    p_value_df1[0.0500] = 3.8400; p_value_df2[0.0500] = 5.9900 
    p_value_df1[0.0250] = 5.0200; p_value_df2[0.0250] = 7.3800 
    p_value_df1[0.0200] = 5.4100; p_value_df2[0.0200] = 7.8200 
    p_value_df1[0.0100] = 6.6300; p_value_df2[0.0100] = 9.2100 
    p_value_df1[0.0050] = 7.8800; p_value_df2[0.0050] = 10.6000 
    p_value_df1[0.0025] = 9.1400; p_value_df2[0.0025] = 11.9800 
    p_value_df1[0.0010] = 10.8300; p_value_df2[0.0010] = 13.8200 
    p_value_df1[0.0005] = 12.1200; p_value_df2[0.0005] = 15.2000
	p_value_df1[0.000001] = 23.9281; p_value_df2[0.000001] = 27.6310
	p_value_df1[0.0000000001] = 41.8215; p_value_df2[0.0000000001] = 46.0517
	p_value_df1[1e-20] = 87.1617; p_value_df2[1e-20] = 92.1034
	p_value_df1[1e-40] = 178.5593; p_value_df2[1e-40] = 184.2068
}
{
    if      (NF==6) p_value_cutoff = p_value_df1[p]; # df == 1; two sample
    else if (NF==8) p_value_cutoff = p_value_df2[p]; # df == 2; three sample
    else    {print "unrecognized format"; exit 1;} 

	if ($2>0 && $4==0 && (NF==6 || $6==0) ) # if the sequence is only associated with Sample1 protein
	{
		printf "%s\t%d\n", $1, $2 > ARGV[1]".Sample1.Enriched";
	} 
	else if ($2==0 && $4>0 && (NF==6 || $6==0)) # if the sequence is only associated with Sample2 protein
	{
		printf "%s\t%d\n", $1, $4 > ARGV[1]".Sample2.Enriched";
	} 
	else if ($2==0 && $4==0 && NF==8 && $6>0 ) # if the sequence is only associated with Sample3 protein
	{
		printf "%s\t%d\n", $1, $6 > ARGV[1]".Sample3.Enriched";
	} 
	else if ($NF >= p_value_cutoff)  # if the sequence is not uniquely associated, test the chi score and p value
	{
        number_of_sample_more_enriched = 0;
        most_enriched_col = 0;
		for (j=2; j<NF; j+=2) 
		{
			if ($j > $(j+1))  # if > E
			{
				++number_of_sample_more_enriched;
				most_enriched_col = j;				
			}
		}
		if (number_of_sample_more_enriched == 1) # the reason p value is small is because one sample has this sequence enriched and the other two have it depeleted 
		{
			if (most_enriched_col==2)  # Sample1
			{
				printf "%s\t%d\n", $1, $2 > ARGV[1]".Sample1.Enriched";
			}
			else if (most_enriched_col==4) # Sample2
			{
				printf "%s\t%d\n", $1, $4 > ARGV[1]".Sample2.Enriched";
			}
			else if (most_enriched_col==6) # Sample3
			{
				printf "%s\t%d\n", $1, $6 > ARGV[1]".Sample3.Enriched";
			}
		}
	}
}
