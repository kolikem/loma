#!/bin/bash

time1=`date +%s`

while getopts I:O:b:s:h:d:l:c:r:m:H:K: OPT
do
	case $OPT in
		"I") flg_I="TRUE"; input_fastq_dir=$OPTARG;;
		"O") flg_O="TRUE"; DIR=$OPTARG;;
		"b") flg_b="TRUE"; block=$OPTARG;;
		"s") flg_s="TRUE"; step=$OPTARG;;
		"h") flg_h="TRUE"; hashicut=$OPTARG;;
		"d") flg_d="TRUE"; n_sigma=$OPTARG;;
		"l") flg_l="TRUE"; lr=$OPTARG;;
		"c") flg_c="TRUE"; ess_min_cov_block=$OPTARG;;
		"r") flg_r="TRUE"; ess_lower_x_percent_discard=$OPTARG;;
		"m") flg_m="TRUE"; ess_paf_clean_match_base_number_lower=$OPTARG;;
		"H") flg_H="TRUE"; minimap2=$OPTARG;;
		"K") flg_K="TRUE"; mafft=$OPTARG;;
	esac
done


if [ "$flg_H" = "TRUE" ]; then echo "-H defined: " $minimap2
else minimap2="minimap2"; echo "-H not defined"; fi
if [ "$flg_K" = "TRUE" ]; then echo "-K defined: " $mafft
else mafft="mafft"; echo "-K not defined"; fi
if [ "$flg_I" = "TRUE" ]; then echo "-I defined: " $input_fastq_dir; fi
if [ "$flg_O" = "TRUE" ]; then echo "-O defined: " $DIR; fi
if [ "$flg_b" = "TRUE" ]; then echo "-b defined: " $block
else block=3000; echo "-b not defined. default value is used: " $block; fi
if [ "$flg_s" = "TRUE" ]; then echo "-s defined: " $step
else step=2000; echo "-b not defined. default value is used: " $step; fi
if [ "$flg_h" = "TRUE" ]; then echo "-h defined: " $hashicut
else hashicut=10; echo "-h not defined. default value is used: " $hashicut; fi
if [ "$flg_d" = "TRUE" ]; then echo "-d defined: " $n_sigma
else n_sigma=3; echo "-d not defined. default value is used: " $n_sigma; fi
if [ "$flg_l" = "TRUE" ]; then echo "-l defined: " $lr
else lr="ont"; echo "-l not defined. default value is used: " $lr; fi
if [ "$flg_c" = "TRUE" ]; then echo "-c defined: " $ess_min_cov_block
else ess_min_cov_block=0.7; echo "-c not defined. default value is used: " $ess_min_cov_block; fi
if [ "$flg_r" = "TRUE" ]; then echo "-r defined: " $ess_lower_x_percent_discard
else ess_lower_x_percent_discard=0.5; echo "-r not defined. default value is used: " $ess_lower_x_percent_discard; fi
if [ "$flg_m" = "TRUE" ]; then echo "-m defined: " $ess_paf_clean_match_base_number_lower
else ess_paf_clean_match_base_number_lower=1000; echo "-m not defined. default value is used: " $ess_paf_clean_match_base_number_lower; fi


curdir=`pwd` #; code_dir=${curdir}/src									# script directory
code_dir=`dirname $0`/src/
##abso=$4
cd ${DIR}
dir1=${DIR}/dir1
mkdir ${dir1}
dir2=${DIR}/dir2
mkdir ${dir2}
CONSENSUS=${DIR}/CONSENSUS
mkdir ${CONSENSUS}
input_fastq2nd_dir=${dir2}/fastq2nd
mkdir ${input_fastq2nd_dir}
abso=0


# Making 1st-consensus sequence.
cd ${dir1}
for file in ${input_fastq_dir}/*fastq; do
	file=`basename $file`;
	name=`echo ${file}|rev|cut -c 11-|rev`;
	$minimap2 -x ava-${lr} ${input_fastq_dir}/$file ${input_fastq_dir}/$file > ${dir1}/${name}_out1.paf;
	echo fastq file  : $file;
	echo region name : $name;
	python3 ${code_dir}/EsS.py ${input_fastq_dir}/$file ${dir1}/${name}_out1.paf 0.04 ${block} ${step} ${ess_min_cov_block} ${ess_lower_x_percent_discard} ${dir1} ${ess_paf_clean_match_base_number_lower};
	for file2 in ${dir1}/${name}_*_out2.fa; do
		file2=`basename $file2`;
		name2=`echo ${file2}|rev|cut -c 9-|rev`;
		$mafft --op 0 --ep 1 --thread 8 --threadit 0 ${dir1}/$file2 > ${dir1}/${name2}_out3.txt;
	done;
	for file3 in ${dir1}/${name}*out3*; do
		file3=`basename $file3`;
		name3=`echo ${file3}|rev|cut -c 10-|rev`;
		python3 ${code_dir}/LineGraphScoring.py ${dir1}/${file3} ${dir1} ${dir1} ${block};
	done;
	for file7 in ${dir1}/${name}*out7*; do
		file7=`basename $file7`;
		name7=`echo ${file7}|rev|cut -c 9-|rev`;
		$mafft --op 0 --ep 1 --thread 8 --threadit 0 ${dir1}/$file7 > ${dir1}/${name7}_out8.txt;
	done;
	for file8 in ${dir1}/${name}*out8*; do
		file8=`basename $file8`;
		name8=`echo ${file8}|rev|cut -c 10-|rev`;
		python3 ${code_dir}/ConsMaker.py ${input_fastq_dir}/$file ${dir1}/$file8 0.6 0.5 0.4 ${dir1}/${name8}_out4.txt;
	done;

	het=`find ${dir1} -name "*${name}*h1*out4*"|wc -l`;
	#rm ${dir1}/${name}*out2*
	#rm ${dir1}/${name}*out3*
	#rm ${dir1}/${name}*out7*
	
	if [ ${het} = 0 ]; then
		echo No heteros found;
		python3 ${code_dir}/sort_out4.py ${dir1} $name ${step};
		cnt_F=`python3 ${code_dir}/CutOut4.py $name ${dir1} F ${hashicut} 0`;
		cnt_B=`python3 ${code_dir}/CutOut4.py $name ${dir1} B ${hashicut} 0`;
		cp ${dir1}/${name}_*_out4_${cnt_F}.txt ${dir1}/consensus_${name}_${cnt_F}-${cnt_B}_${cnt_F};
		python3 ${code_dir}/fasta_name_changer.py ${dir1}/consensus_${name}_${cnt_F}-${cnt_B}_${cnt_F} $cnt_F $cnt_B $cnt_F 1 1 0;	### 1022
		for i in `seq ${cnt_F} $((${cnt_B}-1))`; do
			cat ${dir1}/consensus_${name}_${cnt_F}-${cnt_B}_${i} >> ${dir1}/${name}_$((${i}+1)).swap;
			cat ${dir1}/${name}_*_out4_$((${i}+1)).txt >> ${dir1}/${name}_$((${i}+1)).swap;
			python3 ${code_dir}/RenewConsensus.py ${dir1}/${name}_$((${i}+1)).swap ${dir1} ${block} ${step} $mafft ${dir1}/consensus_${name}_${cnt_F}-${cnt_B}_${i} $((${i}+1));
			python3 ${code_dir}/fasta_name_changer.py ${dir1}/consensus_${name}_${cnt_F}-${cnt_B}_$((${i}+1)) $cnt_F $cnt_B $((${i}+1)) 2 1 0;		### 1022
		done;
		cp ${dir1}/consensus_${name}_* ${dir1}/consensus_${name};
		con=`basename ${dir1}/consensus_${name}_*`;
		mv ${dir1}/${con} ${dir1}/pre_${con};
		cp ${dir1}/consensus_${name} ${CONSENSUS}/consensus_${name};
		#rm ${dir1}/${name}*out8*

	elif [ ${het} != 0 ]; then
		echo Hetero found;
		python3 ${code_dir}/Haplotyping2.py ${dir1} ${name};
		#rm ${dir1}/${name}*out8*
		determined=`find ${dir1} -name "*${name}*h1_typed_out4*"|wc -l`;
		if [ ${determined} = 0 ]; then
			echo Haplotypes are NOT separable.;
		elif [ ${determined} != 0 ]; then
			echo Haplotypes are separable.;
			
			python3 ${code_dir}/sort_out4_haplo.py ${dir1} $name ${step};
			cnt_F=`python3 ${code_dir}/CutOut4.py $name ${dir1} F ${hashicut} 1`;
			cnt_B=`python3 ${code_dir}/CutOut4.py $name ${dir1} B ${hashicut} 1`;
			mkdir ${name}_h1;
			cp ${name}*out4* ${name}_h1;
			#rm ${name}_h1/${name}*_h2_*;
			python3 ${code_dir}/dir_rename.py ${name}_h1 1;
			mkdir ${name}_h2; cp ${name}*out4* ${name}_h2/;
			#rm ${name}_h2/${name}*_h1_*;
			python3 ${code_dir}/dir_rename.py ${name}_h2 2;
			for k in `seq 1 2`; do
				cp ${dir1}/${name}_h${k}/${name}_*_out4_${cnt_F}.txt ${dir1}/${name}_h${k}/consensus_${name}_${cnt_F}-${cnt_B}_${cnt_F};
				python3 ${code_dir}/fasta_name_changer.py ${dir1}/${name}_h${k}/consensus_${name}_${cnt_F}-${cnt_B}_${cnt_F} $cnt_F $cnt_B $cnt_F 1 1 0;		### 1022
				for i in `seq ${cnt_F} $((${cnt_B}-1))`; do
                        		cat ${dir1}/${name}_h${k}/consensus_${name}_${cnt_F}-${cnt_B}_${i} >> ${dir1}/${name}_h${k}/${name}_$((${i}+1)).swap;
                        		cat ${dir1}/${name}_h${k}/${name}_*_out4_$((${i}+1)).txt >> ${dir1}/${name}_h${k}/${name}_$((${i}+1)).swap;
                       			python3 ${code_dir}/RenewConsensus.py ${dir1}/${name}_h${k}/${name}_$((${i}+1)).swap ${dir1}/${name}_h${k} ${block} ${step} $mafft ${dir1}/${name}_h${k}/consensus_${name}_${cnt_F}-${cnt_B}_${i} $((${i}+1));
                			python3 ${code_dir}/fasta_name_changer.py ${dir1}/${name}_h${k}/consensus_${name}_${cnt_F}-${cnt_B}_$((${i}+1)) $cnt_F $cnt_B $((${i}+1)) 2 1 0;           ### 1022
				done;
				cp ${dir1}/${name}_h${k}/consensus_${name}_* ${dir1}/consensus_${name};
				cp ${dir1}/${name}_h${k}/consensus_${name}_* ${dir1}/${name}_h${k}/consensus_${name}_hap${k};
				cp ${dir1}/consensus_${name} ${CONSENSUS}/consensus_${name};
			done;
		fi;
	fi;
	#if [ -e ${dir1}/${name}_h1 ]; then
		#rm ${dir1}/${name}_h1/*out4*;
		#rm ${dir1}/${name}_h1/*swap*;
		#rm ${dir1}/${name}_h2/*out4*;
		#rm ${dir1}/${name}_h2/*swap*;
	#fi;
	#rm ${dir1}/*out4*; rm ${dir1}/*swap*;
done

time2=`date +%s`

#<< COMMENTOUT
# Classifying heterozygous fragments.
input_sam2nd_dir=${dir2}/sam2nd
mkdir ${input_sam2nd_dir}
for cons in ${dir1}/consensus*; do
	name=`basename ${cons}`; name=`echo ${name}|cut -c 11-`; echo $name;
        $minimap2 -a ${cons} ${input_fastq_dir}/${name}.sam.fastq > ${input_sam2nd_dir}/${name}_map_on_1st_consensus.sam;	
	python3 ${code_dir}/ReadClassify6.py ${input_sam2nd_dir}/${name}_map_on_1st_consensus.sam 0.1 500 100 100 8 ${input_fastq_dir}/${name}.sam.fastq ${dir2} ${input_fastq2nd_dir} ${n_sigma} ${abso};
done

time3=`date +%s`
 
# Making 2nd-consensus sequence.
cd ${dir2}
nakami=`ls ${input_fastq2nd_dir}`
if [ -z ${nakami} ]; then echo 'dir2のfastq2ndが空のため終了.'; else
for file in ${input_fastq2nd_dir}/*fastq.*; do
	file=`basename $file`;
	work=${dir2}/work_${file};
	mkdir ${work};
	cd ${work};
	name=`echo ${file}|rev|cut -c 19-|rev`;
	HAP=`echo ${file}|rev|cut -c 1|rev`;
	$minimap2 -x ava-${lr} ${input_fastq2nd_dir}/$file ${input_fastq2nd_dir}/$file > ${work}/${name}_out1.paf;
	echo fastq file  : $file;
	echo region name : $name;
	python3 ${code_dir}/EsS.py ${input_fastq2nd_dir}/$file ${work}/${name}_out1.paf 0.04 ${block} ${step} ${ess_min_cov_block} ${ess_lower_x_percent_discard} ${work} ${ess_paf_clean_match_base_number_lower};
	for file2 in ${work}/${name}_*_out2.fa; do
		file2=`basename $file2`;
		name2=`echo ${file2}|rev|cut -c 9-|rev`;
		$mafft --op 0 --ep 1 --thread 8 --threadit 0 ${work}/$file2 > ${work}/${name2}_out3.txt;
	done;
	for file3 in ${work}/${name}*out3*; do
		file3=`basename $file3`;
		name3=`echo ${file3}|rev|cut -c 10-|rev`;
		python3 ${code_dir}/LineGraphScoring.py ${work}/${file3} ${work} ${work} ${block};
	done;
	for file7 in ${work}/${name}*out7*; do
		file7=`basename $file7`;
		name7=`echo ${file7}|rev|cut -c 9-|rev`;
		$mafft --op 0 --ep 1 --thread 8 --threadit 0 ${work}/$file7 > ${work}/${name7}_out8.txt;
	done;
	for file8 in ${work}/${name}*out8*; do
		file8=`basename $file8`;
		name8=`echo ${file8}|rev|cut -c 10-|rev`;
		python3 ${code_dir}/ConsMaker.py ${input_fastq2nd_dir}/$file ${work}/$file8 0.6 0.5 0.4 ${work}/${name8}_out4.txt;
	done;

	het=`find ${work} -name "*${name}*h1*out4*"|wc -l`;
	rm ${work}/${name}*out2*
	rm ${work}/${name}*out3*
	rm ${work}/${name}*out7*
	
	if [ ${het} = 0 ]; then
		echo No heteros found;
		python3 ${code_dir}/sort_out4.py ${work} $name ${step};
		cnt_F=`python3 ${code_dir}/CutOut4.py $name ${work} F ${hashicut} 0`;
		cnt_B=`python3 ${code_dir}/CutOut4.py $name ${work} B ${hashicut} 0`;
		cp ${work}/${name}_*_out4_${cnt_F}.txt ${work}/consensus_${name}_${cnt_F}-${cnt_B}_${cnt_F};
		python3 ${code_dir}/fasta_name_changer.py ${work}/consensus_${name}_${cnt_F}-${cnt_B}_${cnt_F} $cnt_F $cnt_B $cnt_F 1 2 $HAP;
		for i in `seq ${cnt_F} $((${cnt_B}-1))`; do
                        cat ${work}/consensus_${name}_${cnt_F}-${cnt_B}_${i} >> ${work}/${name}_$((${i}+1)).swap;
                        cat ${work}/${name}_*_out4_$((${i}+1)).txt >> ${work}/${name}_$((${i}+1)).swap;
                        python3 ${code_dir}/RenewConsensus.py ${work}/${name}_$((${i}+1)).swap ${work} ${block} ${step} $mafft ${work}/consensus_${name}_${cnt_F}-${cnt_B}_${i} $((${i}+1));
                	python3 ${code_dir}/fasta_name_changer.py ${work}/consensus_${name}_${cnt_F}-${cnt_B}_$((${i}+1)) $cnt_F $cnt_B $((${i}+1)) 2 2 $HAP;
		done;
		cp ${work}/consensus_${name}_* ${CONSENSUS}/consensus_${file};
		rm ${work}/${name}*out8*

	elif [ ${het} != 0 ]; then
		echo Hetero found;
		python3 ${code_dir}/Haplotyping2.py ${work} ${name};
		rm ${work}/${name}*out8*
		determined=`find ${work} -name "*${name}*h1_typed_out4*"|wc -l`;
		if [ ${determined} = 0 ]; then
			echo Haplotypes are NOT separable.;
		elif [ ${determined} != 0 ]; then
			echo Haplotypes are separable.;
			
			python3 ${code_dir}/sort_out4_haplo.py ${work} $name ${step};
			cnt_F=`python3 ${code_dir}/CutOut4.py $name ${work} F ${hashicut} 1`;
			cnt_B=`python3 ${code_dir}/CutOut4.py $name ${work} B ${hashicut} 1`;
			mkdir ${name}_h1;
			cp ${name}*out4* ${name}_h1/;
			rm ${name}_h1/${name}*_h2_*;
			python3 ${code_dir}/dir_rename.py ${name}_h1 1;
			mkdir ${name}_h2; cp ${name}*out4* ${name}_h2/;
			rm ${name}_h2/${name}*_h1_*;
			python3 ${code_dir}/dir_rename.py ${name}_h2 2;
			for k in `seq 1 2`; do
				cp ${work}/${name}_h${k}/${name}_*_out4_${cnt_F}.txt ${work}/${name}_h${k}/consensus_${name}_${cnt_F}-${cnt_B}_${cnt_F};
				python3 ${code_dir}/fasta_name_changer.py ${work}/${name}_h${k}/consensus_${name}_${cnt_F}-${cnt_B}_${cnt_F} $cnt_F $cnt_B $cnt_F 1 2 $HAP;
                                for i in `seq ${cnt_F} $((${cnt_B}-1))`; do
                                        cat ${work}/${name}_h${k}/consensus_${name}_${cnt_F}-${cnt_B}_${i} >> ${work}/${name}_h${k}/${name}_$((${i}+1)).swap;
                                        cat ${work}/${name}_h${k}/${name}_*_out4_$((${i}+1)).txt >> ${work}/${name}_h${k}/${name}_$((${i}+1)).swap;
                                        python3 ${code_dir}/RenewConsensus.py ${work}/${name}_h${k}/${name}_$((${i}+1)).swap ${work}/${name}_h${k} ${block} ${step} $mafft ${work}/${name}_h${k}/consensus_${name}_${cnt_F}-${cnt_B}_${i} $((${i}+1));
					python3 ${code_dir}/fasta_name_changer.py ${work}/${name}_h${k}/consensus_${name}_${cnt_F}-${cnt_B}_$((${i}+1)) $cnt_F $cnt_B $((${i}+1)) 2 2 $HAP;
                                done;
				cp ${work}/${name}_h${k}/consensus_${name}_* ${CONSENSUS}/consensus_${file};
				cp ${work}/${name}_h${k}/consensus_${name}_* ${work}/${name}_h${k}/consensus_${name}_hap${k};
			done;
		fi;
	fi;
	if [ -e ${work}/${name}_h1 ]; then
		rm ${work}/${name}_h1/*out4*;
		rm ${work}/${name}_h1/*swap*;
		rm ${work}/${name}_h2/*out4*;
		rm ${work}/${name}_h2/*swap*;
	fi;
	rm ${work}/*swap*; rm ${dir2}/*out4*;
done
fi
#COMMENTOUT

time4=`date +%s`


time_firstCS=$((time2-time1))
time_readSep=$((time3-time2))
time_secondCS=$((time4-time3))
time_total=$((time4-time1))
echo time to finish 1st CS: $time_firstCS
echo time to finish read classification: $time_readSep
echo time to finish 2nd CS: $time_secondCS
echo total time: $time_total

