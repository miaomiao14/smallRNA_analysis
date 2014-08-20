
pipeline_Dir=${HOME}/git/smallRNA_analysis/primary_piRNA_biogenesis
function useZuc {
	# ./run_detail.sh  Phil.SRA.AubIP.42A18-zucHM-zucSG63.rep1.ox.ovary.trimmed.x_rRNA.x_hairpin.x_42AB18_masked_by_w1_L1read.dm3v0.all.x_rpmk_MASK.bed2 Phil.SRA.Ago3IP.42A18-zucHM-zucSG63.rep1.ox.ovary.trimmed.x_rRNA.x_hairpin.x_42AB18_masked_by_w1_L1read.dm3v0.all.x_rpmk_MASK.bed2 $1 ${1/Watson/Crick}
	# ./run_detail.sh  Phil.SRA.Ago3IP.42A18-zucHM-zucSG63.rep1.ox.ovary.trimmed.x_rRNA.x_hairpin.x_42AB18_masked_by_w1_L1read.dm3v0.all.x_rpmk_MASK.bed2 Phil.SRA.AubIP.42A18-zucHM-zucSG63.rep1.ox.ovary.trimmed.x_rRNA.x_hairpin.x_42AB18_masked_by_w1_L1read.dm3v0.all.x_rpmk_MASK.bed2 $1 ${1/Watson/Crick}
	
	#bash -x ./run_detail1.sh Ago3IPuniq_w1_unox_.AubIPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.cis+trans.bed _ $1 ${1/Watson/Crick} $2
	# bash -x ./run_detail1.sh Ago3IPuniq_w1_unox_.AubIPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.trans.bed _ $1 ${1/Watson/Crick} $2
	bash -x ${pipeline_Dir}/run_detail1.sh AubIPuniq_w1_unox_.Ago3IPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.cis.bed _ $1 ${1/Watson/Crick} $2
	# ./run_detail.sh Phil.SRA.AubIP.+_42A18.rep1.ox.ovary.trimmed.x_rRNA.x_hairpin.x_42AB18_masked_by_w1_L1read.dm3v0.all.x_rpmk_MASK.bed2 Phil.SRA.Ago3IP.+_42A18.rep1.ox.ovary.trimmed.x_rRNA.x_hairpin.x_42AB18_masked_by_w1_L1read.dm3v0.all.x_rpmk_MASK.bed2  $1 ${1/Watson/Crick}
}
useZuc Phil.SRA.PiwiIPuniq3.w1.unox.ovary.x_rRNA.x_hairpin.dm3v0.all.piRNA.sorted.Watson.bigWig
useZuc Phil.SRA.Ago3IPuniq3.w1.unox.ovary.x_rRNA.x_hairpin.dm3v0.all.piRNA.sorted.Watson.bigWig
useZuc Phil.SRA.AubIPuniq3.w1.unox.ovary.x_rRNA.x_hairpin.dm3v0.all.piRNA.sorted.Watson.bigWig

