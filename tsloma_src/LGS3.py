# Written by python. 23.01.10
# LineGraphScoring.pyを改変しここではヘテロ分類しないようにする.

print("---LGS.py---")

import os
import re
import gc
import sys
import numpy as np
import matplotlib.pyplot as plt
print("imported modules")

# sys.argv[1] - mafft output (out3)
# sys.argv[2] - 結果図格納場所・分類後fastaファイル保存ディレクトリ

def preprocess():
    f = open(sys.argv[1])	# e.g. ST3GAL3_novel_0_ATL7.fastq.out3.2
    l = f.readlines()
    name_num, names, reads = [], [], []
    for i in range(len(l)):
        if l[i][0] == '>':
            name_num.append(i)
            names.append(l[i][1:-1])

    for i in range(len(name_num)):
        if i == len(name_num)-1:
            seq = ''
            for j in range(name_num[i]+1, len(l)):
                line = l[j].replace('\n', '')
                seq += line
            reads.append(seq)       
        else:
            seq = ''
            for j in range(name_num[i]+1, name_num[i+1]):
                line = l[j].replace('\n', '')
                seq += line
            reads.append(seq)   #リードのアラインメントされた配列(-込み)を文字列で順番にリストに格納
    f.close()
    del l
    gc.collect()

    return names, reads

def end_due_to_small_read_number(input_file, out_dir, read_num, names, reads):
    name = os.path.basename(input_file)
    if read_num == 1:		# ﾘｰﾄﾞ本数 = 1 --> out4 (sub consensus)
        new_name = name.replace("out3", "out4")
        with open(f'{out_dir}/{new_name}', mode='w') as f:
            s = names[0] + "\n" + reads[0] + "\n"
            f.write(s)
    else:
        new_name = name.replace("out3", "out8")
        with open(f'{out_dir}/{new_name}', mode='w') as f:
            s = ''
            for i in range(read_num):
                s += f'>{names[i]}\n'
                read_i = reads[i]
                s += f'{read_i}\n'
            f.write(s)
    print(name+" was saved in "+out_dir)
    return 0

def syukei_atgc(r, i):
	# 多重アラインメントi番目の{A,T,G,C,-}数をカウントする
	l = [r[k][i] for k in range(len(r))]
	d = {"a":l.count("a")/len(r), "t":l.count("t")/len(r), "g":l.count("g")/len(r), "c":l.count("c")/len(r), "-":l.count("-")/len(r)}
	# 割合の辞書を返す
	return d

def ErrorMatrix(reads):
	# 行数(ﾘｰﾄﾞ数)*列数(ｱﾗｲﾝﾒﾝﾄ長)の行列 EM_ji: j番目ﾘｰﾄﾞのi番目塩基のi番目塩基全体における割合
	EM = [[0 for i in range(len(reads[0]))] for i in range(len(reads))]
	for i in range(len(reads[0])):
		prop_i = syukei_atgc(reads, i)	#i番目の塩基割合 prop_i = {"a": 0.1, "t": 0.1, "g": 0.4, "c": 0.2, "-": 0.2}
		for j in range(len(reads)):
			EM[j][i] = prop_i[reads[j][i]]
	return EM

def sum_EM(em):
	# EM ﾘｰﾄﾞごとの和
	l = [0]*len(em)
	for j in range(len(em)):
		l_j = sum(em[j])
		l[j] = l_j
	return l

def SD(l):
	# リストlの標準偏差を計算する.
	l_sq = [i**2 for i in l]
	var = sum(l_sq)/len(l) - (sum(l)/len(l))**2
	return var**0.5

def sd_deviated(l, n):
	# sd*n個以上ずれているlistの値を検出する
	idx_detected = []
	ave_l = sum(l)/len(l)
	sd_l = SD(l)
	for i in range(len(l)):
		if ave_l-l[i] > sd_l*n:	# if low score
			idx_detected.append(i)
	return idx_detected

def print_deviated_reads(idx, names):
	print("The following reads are excluded due to more than 2SDs deviation.")
	for i in idx:
		print(names[i])
	return 0

def save_non_error_reads_as_fasta(input_file, out_dir, names, reads, idx):
	file_name = os.path.basename(input_file)
	new_file_name = file_name.replace("out3", "out7")
	content = ""
	for i in range(len(names)):
		if i in idx:
			continue
		else:
			content += ">" + names[i] + "\n" + reads[i].replace("-","") + "\n"
	with open(f'{out_dir}/{new_file_name}', mode="w") as f:
		f.write(content)
	print(new_file_name+" was saved in "+out_dir)
	return 0


if __name__ == "__main__":
	divide = 10	# ｴﾗｰﾘｰﾄﾞ検出を行うかどうかの境界本数.	
	names, reads = preprocess()
	
	if len(names) < divide:
		print("The number of this file < "+str(divide)+". Error reads are NOT explored.")
		end_due_to_small_read_number(sys.argv[1], sys.argv[2], len(names), names, reads)

	else:
		print("The number of this file < "+str(divide)+". Error reads are being explored.")
		EM = ErrorMatrix(reads)
		sumEM = sum_EM(EM)
		idx_detected = sd_deviated(sumEM, 2)
		print("#reads: ", len(reads))
		print("length of aligned length: ", len(reads[0]))

		print_deviated_reads(idx_detected, names)			# ｴﾗｰﾘｰﾄﾞ名を出力
		save_non_error_reads_as_fasta(sys.argv[1], sys.argv[2], names, reads, idx_detected)	# 非ｴﾗｰﾘｰﾄﾞをfastaで保存


