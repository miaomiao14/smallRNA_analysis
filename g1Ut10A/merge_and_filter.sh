#!/usr/bin/env bash

export RUN_DIRECTORY=$(dirname `readlink -f $0`)
export PATH=${RUN_DIRECTORY}:$PATH

usage () {    
cat << EOF
This script takes IP data from PIWI, AUB and AGO3 then apply chi-squared test to assign reads to each PIWI protein.
You have to provide at least two IPs samples from different PIWI protein.
It currently does not store the table for converting chi score to p value, thus only 0.05 and 0.005 is supported.

usage:

    merge_and_filter.sh \ 
        -w Aravin.PiwiIP.SRA.MD.trimmed.fq.gz \ 
        -b Aravin.AubIP.SRA.MD.trimmed.fq.gz \ 
        -g Aravin.Ago3IP.SRA.MD.trimmed.fq.gz \ 
        -P 0.005 \ 
        -p 3 \ 
        -o Aravin.MD

OPTIONS:
    -w  Piwi bound piRNAs in fq or fq.gz format
    -b  Aub  bound piRNAs in fq or fq.gz format
    -g  Ago3 bound piRNAs in fq or fq.gz format
    -P  p-value cutoff to be used. Use only 0.05 or 0.005. default: 0.05
    -o  Output file prefix, default: Piwi-Aub-Ago3
    -p  CPU to use; default: 3
    
EOF
}

while getopts "w:b:g:P:p:o:h" OPTION; do
	case $OPTION in
		h)	usage && exit 0 ;;
        w)  PIWI_FQ=$OPTARG;;
        b)  AUB_FQ=$OPTARG;;
        g)  AGO3_FQ=$OPTARG;;
        P)  P_VALUE=$OPTARG;;
        p)  CPU=$OPTARG;;
        o)  PREFIX=$OPTARG;;
		*)	usage && exit 1 ;;
	esac
done

# set P value
[[ -z "$P_VALUE" ]] && P_VALUE=0.05
[[ "$P_VALUE" != "0.05" && "$P_VALUE" != "0.005" ]] && echo "currently -P only takes 0.05 or 0.005. Use 0.05." && P_VALUE=0.05

# test for number of IP
COUNT=0
[[ ! -z "$PIWI_FQ" ]] && COUNT=$((COUNT+1))
[[ ! -z "$AUB_FQ"  ]] && COUNT=$((COUNT+1))
[[ ! -z "$AGO3_FQ" ]] && COUNT=$((COUNT+1))
[[ $COUNT -lt 2 ]] && echo "have to set at least two IPs" && exit 1

# test for file existence
[[ ! -z "$PIWI_FQ" && ! -f "$PIWI_FQ" ]] && echo "cannot find $PIWI_FQ" && exit 1;
[[ ! -z "$AUB_FQ"  && ! -f "$AUB_FQ" ]]  && echo "cannot find $AUB_FQ"  && exit 1;
[[ ! -z "$AGO3_FQ" && ! -f "$AGO3_FQ" ]] && echo "cannot find $AGO3_FQ" && exit 1;

# set CPU
[[ -z "${CPU##*[!0-9]*}" ]] && CPU=3

# set prefix
[[ ! -z "$PIWI_FQ" ]] && PREFIX=${PREFIX}".Piwi"

if [[ ! -z "$AUB_FQ"  ]]; then
    if [[ ! -z "$PIWI_FQ" ]]; then # Piwi is set
        PREFIX=${PREFIX}"-Aub"
    else # Piwi is not set
        PREFIX=${PREFIX}".Aub"
    fi
fi

[[ ! -z "$AGO3_FQ" ]] && PREFIX=${PREFIX}"-Ago3"

PIWI_IST=$PIWI_FQ
AUB_IST=$AUB_FQ
AGO3_IST=$AGO3_FQ

# echo "converting fq to insert file"
#
# PARAFILE=${RANDOM}${RANDOM}${RANDOM}.para
# [[ ! -z "$PIWI_FQ" ]] && PIWI_IST=${PIWI_FQ%.f[aq]*}.insert && echo "piPipes_fastq_to_insert $PIWI_FQ $PIWI_IST" >> $PARAFILE
# [[ ! -z "$AUB_FQ"  ]] && AUB_IST=${AUB_FQ%.f[aq]*}.insert   && echo "piPipes_fastq_to_insert $AUB_FQ  $AUB_IST"  >> $PARAFILE
# [[ ! -z "$AGO3_FQ" ]] && AGO3_IST=${AGO3_FQ%.f[aq]*}.insert && echo "piPipes_fastq_to_insert $AGO3_FQ $AGO3_IST" >> $PARAFILE
# ParaFly -c $PARAFILE -CPU $CPU

echo "Merging insert files"
insertMerge -i $PIWI_IST $AUB_IST $AGO3_IST -c 1 -t 2 > ${PREFIX}.merged_insert # this step doesn't care about number of sample

echo "calculating expected value and chi scores"
chiSquared ${PREFIX}.merged_insert > ${PREFIX}.merged_insert.chi && rm -rf ${PREFIX}.merged_insert # this step doesn't care about number of sample

echo "filter values according to chi-squared test"
awk -v p=$P_VALUE -f $RUN_DIRECTORY/filter_chi_new.awk ${PREFIX}.merged_insert.chi

# echo "convert insert back to gzipped fq"
# PARAFILE=${RANDOM}${RANDOM}${RANDOM}.para
# if [[ ! -z "$PIWI_FQ" ]] ; then # if there is Piwi, it must be the first one
#     echo "insert2fq ${PREFIX}.merged_insert.chi.Sample1.Enriched | gzip > ${PREFIX}.Piwi.fq.gz"  >> $PARAFILE
#     if [[ ! -z "$AUB_FQ" ]] ; then # if there is Piwi and Aub, Aub must be the second one
#         echo "insert2fq ${PREFIX}.merged_insert.chi.Sample2.Enriched | gzip > ${PREFIX}.Aub.fq.gz"  >> $PARAFILE
#         [[ ! -z "$AGO3_FQ" ]] && echo "insert2fq ${PREFIX}.merged_insert.chi.Sample3.Enriched | gzip > ${PREFIX}.Ago3.fq.gz" >> $PARAFILE
#     else
#         # no Aub, must be Piwi and Ago3
#         echo "insert2fq ${PREFIX}.merged_insert.chi.Sample2.Enriched  | gzip > ${PREFIX}.Ago3.fq.gz"  >> $PARAFILE
#     fi
# else
#     # no Piwi, must be Aub and Ago3
#     echo "insert2fq ${PREFIX}.merged_insert.chi.Sample1.Enriched  | gzip > ${PREFIX}.Aub.fq.gz"  >> $PARAFILE
#     echo "insert2fq ${PREFIX}.merged_insert.chi.Sample2.Enriched  | gzip > ${PREFIX}.Ago3.fq.gz" >> $PARAFILE
# fi
# ParaFly -c $PARAFILE -CPU $CPU
