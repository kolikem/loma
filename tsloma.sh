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
else step=2000; echo "-s not defined. default value is used: " $step; fi
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


curdir=`pwd`
code_dir=`dirname $0`/tsloma_src
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


# Step1. 1st-CS.
cd ${dir1}
for file in ${input_fastq_dir}/*fastq; do
	file=`basename $file`;
	name=`echo ${file}`;
	$minimap2 -x ava-${lr} ${input_fastq_dir}/$file ${input_fastq_dir}/$file > ${dir1}/${name}.out1;
	echo fastq file  : $file;
	echo region name : $name;
	echo "Runnig command: python3 ${code_dir}/EsS.py ${input_fastq_dir}/$file ${dir1}/${name}.out1 0.04 ${block} ${step} ${ess_min_cov_block} ${ess_lower_x_percent_discard} ${dir1} ${ess_paf_clean_match_base_number_lower}";
	python3 ${code_dir}/EsS2.py ${input_fastq_dir}/$file ${dir1}/${name}.out1 0.04 ${block} ${step} ${ess_min_cov_block} ${ess_lower_x_percent_discard} ${dir1} ${ess_paf_clean_match_base_number_lower};
	for file2 in ${dir1}/${name}.out2.*; do
		file2=`basename $file2`;
		number=`echo ${file2}| rev`; number=(${number//./ }); number=`echo ${number[0]}| rev`;
		name2=`echo ${file2}| sed -e "s/out2.${number}/out3.${number}/"`;
		$mafft --op 0 --ep 1 --thread 8 --threadit 0 ${dir1}/$file2 > ${dir1}/$name2;
	done;
	for file3 in ${dir1}/${name}.out3.*; do
		file3=`basename $file3`;
		python3 ${code_dir}/LGS3.py ${dir1}/${file3} ${dir1} ${dir1} ${block};
	done;
	for file7 in ${dir1}/${name}.out7.*; do
		file7=`basename $file7`;
		number=`echo ${file7}| rev`; number=(${number//./ }); number=`echo ${number[0]}| rev`;
		name8=`echo ${file7}| sed -e "s/out7.${number}/out8.${number}/"`;
		echo "Running command: $mafft --op 0 --ep 1 --thread 8 --threadit 0 ${dir1}/$file7 > ${dir1}/${name8}";
		$mafft --op 0 --ep 1 --thread 8 --threadit 0 ${dir1}/$file7 > ${dir1}/${name8};
	done;
	for file8 in ${dir1}/${name}.out8.*; do
		file8=`basename $file8`;
		number=`echo ${file8}| rev`; number=(${number//./ }); number=`echo ${number[0]}| rev`;
		name4=`echo ${file8}| sed -e "s/out8.${number}/out4.${number}/"`;
		python3 ${code_dir}/SCM2.py ${input_fastq_dir}/$file ${dir1}/$file8 0.6 0.5 0.4 ${dir1}/${name4};
	done;

	rm ${dir1}/${name}*out2*
	rm ${dir1}/${name}*out3*
	rm ${dir1}/${name}*out7*
	
	cnt_out4=`ls -1 ${dir1}/${name}.out4.*| wc -l`;
	python3 ${code_dir}/RCS2.py ${dir1}/${name} ${cnt_out4} ${block} ${step} $mafft;
	rm ${dir1}/${name}*out8*
	rm ${dir1}/*.concat.*
	rm ${dir1}/*.out4.*
	rm ${dir1}/*.out4a.*
	rm ${dir1}/*.tail_head.*
done

time2=`date +%s`


# Step2. hetero-class.
input_sam2nd_dir=${dir2}/sam2nd
mkdir ${input_sam2nd_dir}
for cons in ${dir1}/*.cs; do
	name=`basename ${cons}`; name=`echo ${name}|rev|cut -c 4-|rev`;
        $minimap2 -a ${cons} ${input_fastq_dir}/${name} > ${input_sam2nd_dir}/${name}_map_on_1st_cs.sam;	
	python3 ${code_dir}/ReadClassify6.py ${input_sam2nd_dir}/${name}_map_on_1st_cs.sam 0.1 500 100 100 8 ${input_fastq_dir}/${name} ${dir2} ${input_fastq2nd_dir} ${n_sigma} ${abso};
done
time3=`date +%s`


# Step3. 2nd-CS.
cd ${dir2}
nakami=`ls ${input_fastq2nd_dir}`
if [ -z "${nakami}" ]; then echo '3rd step finished. No hetero region detected.'; else echo 'Hetero region detected.';
for file in ${input_fastq2nd_dir}/*fastq*; do
	file=`basename $file`;
	name=`echo ${file}`;
	$minimap2 -x ava-${lr} ${input_fastq2nd_dir}/$file ${input_fastq2nd_dir}/$file > ${dir2}/${file}.out1;
	echo fastq file  : $file;
	#echo region name : $name;
	echo "Runnig command: python3 ${code_dir}/EsS2.py ${input_fastq2nd_dir}/$file ${dir2}/${file}.out1 0.04 ${block} ${step} ${ess_min_cov_block} ${ess_lower_x_percent_discard} ${dir2} ${ess_paf_clean_match_base_number_lower}";
	python3 ${code_dir}/EsS2.py ${input_fastq2nd_dir}/$file ${dir2}/${file}.out1 0.04 ${block} ${step} ${ess_min_cov_block} ${ess_lower_x_percent_discard} ${dir2} ${ess_paf_clean_match_base_number_lower};
	for file2 in ${dir2}/${file}.out2.*; do
		file2=`basename $file2`;
		number=`echo ${file2}| rev`; number=(${number//./ }); number=`echo ${number[0]}| rev`;
		name2=`echo ${file2}| sed -e "s/out2.${number}/out3.${number}/"`;
		$mafft --op 0 --ep 1 --thread 8 --threadit 0 ${dir2}/$file2 > ${dir2}/$name2;
	done;
	for file3 in ${dir2}/${file}.out3.*; do
		file3=`basename $file3`;
		python3 ${code_dir}/LGS3.py ${dir2}/${file3} ${dir2} ${dir2} ${block};
	done;
	for file7 in ${dir2}/${file}.out7.*; do
		file7=`basename $file7`;
		number=`echo ${file7}| rev`; number=(${number//./ }); number=`echo ${number[0]}| rev`;
		name8=`echo ${file7}| sed -e "s/out7.${number}/out8.${number}/"`;
		echo "Running command: $mafft --op 0 --ep 1 --thread 8 --threadit 0 ${dir2}/$file7 > ${dir2}/${name8}";
		$mafft --op 0 --ep 1 --thread 8 --threadit 0 ${dir2}/$file7 > ${dir2}/${name8};
	done;
	for file8 in ${dir2}/${file}.out8.*; do
		file8=`basename $file8`;
		number=`echo ${file8}| rev`; number=(${number//./ }); number=`echo ${number[0]}| rev`;
		name4=`echo ${file8}| sed -e "s/out8.${number}/out4.${number}/"`;
		python3 ${code_dir}/SCM2.py ${input_fastq2nd_dir}/$file ${dir2}/$file8 0.6 0.5 0.4 ${dir2}/${name4};
	done;

	rm ${dir2}/${file}*out2*
	rm ${dir2}/${file}*out3*
	rm ${dir2}/${file}*out7*
	
	cnt_out4=`ls -1 ${dir2}/${name}.out4.*| wc -l`;
	python3 ${code_dir}/RCS2.py ${dir2}/${file} ${cnt_out4} ${block} ${step} $mafft;
	rm ${dir2}/${file}*out8*
	rm ${dir2}/*.concat.*
	rm ${dir2}/*.out4.*
	rm ${dir2}/*.out4a.*
	rm ${dir2}/*.tail_head.*
done
fi
time4=`date +%s`

mv ${dir1}/*.cs ${CONSENSUS}
if [ -e ${dir2}/*.cs ]; then mv ${dir2}/*.cs ${CONSENSUS}; fi

time_firstCS=$((time2-time1))
time_readSep=$((time3-time2))
time_secondCS=$((time4-time3))
time_total=$((time4-time1))
echo time to finish 1st CS: $time_firstCS
echo time to finish read classification: $time_readSep
echo time to finish 2nd CS: $time_secondCS
echo total time: $time_total

