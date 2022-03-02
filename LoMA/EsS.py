"""//////////////////////////////////////////////////////////////////////////////////////////
                        リードを重ねて共通部分を切り出すコード(EsS)

ーーーコマンドライン引数の説明ーーー
sys.argv[1] - リードのfastq
sys.argv[2] - minimap2によるアラインメントの結果paf
sys.argv[3] - アラインメントされない端の許容範囲
sys.argv[4] - ブロックサイズ
sys.argv[5] - ステップサイズ
sys.argv[6] - FASTAファイルに入れるためのblockにかかっている長さの最低ライン
sys.argv[7] - アラインメント長下位 x % の切り捨て
sys.argv[8] - output


To run, python3 EsS.py <fq> <paf> <hashikko> <block size> <step size> <length restriction>

*** 細かいコーディングについては「研究ノート１」にメモ
//////////////////////////////////////////////////////////////////////////////////////////"""

print("---EsS.py---")

import re
import collections
import sys
import copy
import random
import numpy as np
import os
import time
startTime = time.time()

######ここからデータの成形######

def getSeqDiv_from_paf_format(line):
	if "\tdv:f:" in line:
		m = re.search(r'\tdv:f:.*\t', line)
		m = m.group()
		seq_div = float(m[6:-1])
	elif "\tde:f:" in line:
		m = re.search(r'\tde:f:.*\t', line)
		m = m.group()
		seq_div = float(m[6:-1])
	else:
		seq_div = "*"
	return seq_div


		



'''（１）'''

# Import reads from FASTQ file.
file_fastq = open(sys.argv[1])
list_fastq = file_fastq.readlines()
fq_seq = []          #リード配列
fq_name = []         #リード名
for i in range(len(list_fastq)):
    if i % 4 == 1:
        fq_seq.append(list_fastq[i])
    if i % 4 == 0:
        fq_name.append(list_fastq[i])
file_fastq.close()

for i in range(len(fq_seq)):
    fq_seq[i] = fq_seq[i].rstrip('\n')
for i in range(len(fq_name)):
    fq_name[i] = fq_name[i].rstrip('\n')


fq_name0 = []
for i in range(len(fq_name)):
    s = fq_name[i].split(' ')
    fq_name0.append(s[0][1:])


fq_name = fq_name0      #fq_name: All reads' names as list.
numFq = len(fq_name)    #numFq: Number of reads, in FASTQ file.

#for i in enumerate(fq_name):
#    print(i)
print('fastq内リード数:', numFq, '本')


'''（２）'''

# Import the result of all-versus-all alignments from PAF.
PAF_data_INIT = np.array([])


list_paf_name1 = []     #リード１の名前
length_1 = []
list_paf_start1 = []    #リード１のスタート位置
list_paf_end1 = []      #リード１の終了位置
list_paf_relative = []  #リード１とリード２が＋か-（+ = 1, - = -1）
list_paf_name2 = []     #リード２の名前
length_2 = []
list_paf_start2 = []    #リード２のスタート位置
list_paf_end2 = []      #リード２の終了位置
list_paf_matchBase = []	#アラインメントのマッチ塩基数
list_paf_alnLen = []	#アラインメント長
list_paf_seqDiv = []

with open(sys.argv[2]) as file_paf:
	for line in file_paf:
		seq_div = getSeqDiv_from_paf_format(line)
		list_paf_seqDiv.append(seq_div)

		line_split = re.split('\t', line)
			
		for i in range(len(line_split)):
			if i % 17 == 0 and line_split[i] != '':
				list_paf_name1.append(line_split[i])
			if i % 17 == 1:
				length_1.append(int(line_split[i]))
			if i % 17 == 2:
				list_paf_start1.append(int(line_split[i]))
			if i % 17 == 3:
				list_paf_end1.append(int(line_split[i]))
			if i % 17 == 4:
				list_paf_relative.append(line_split[i])
			if i % 17 == 5:
				list_paf_name2.append(line_split[i])
			if i % 17 == 6:
				length_2.append(int(line_split[i]))
			if i % 17 == 7:
				list_paf_start2.append(int(line_split[i]))
			if i % 17 == 8:
				list_paf_end2.append(int(line_split[i]))
			if i % 17 == 9:
				list_paf_matchBase.append(int(line_split[i]))
			if i % 17 == 10:
				list_paf_alnLen.append(int(line_split[i]))

#PAF_data_INIT = np.array(list_paf_name1)
#PAF_data_INIT = 


#For relative, replace, + with 1, - with -1.
for i in range(len(list_paf_relative)):
    if list_paf_relative[i] == '+':
        list_paf_relative[i] = 1
    else:
        list_paf_relative[i] = -1





'''（３）端のアラインメント不成の長さのフィルター　＜フィルター１＞'''

hashikko = float(sys.argv[3])   #Define both ends of read,to specify dangling reads.
#9/10 変更点　端の定義を割合ではなく絶対値にする。割合だと長さ10kbのリードは5%でも500bの不斉を許容することになるから。
#hashikko = int(sys.argv[3])

list_ox = []    #o for inclusion, x for exclusion, as list

#Hollow inside cases (,in which the middles of two reads are not aligned although both ends are filled) will be considered in 2nd filter.
for i in range(len(length_1)):

    if list_paf_relative[i] == 1:
        if (list_paf_start1[i] / length_1[i]) >= hashikko and (list_paf_start2[i] / length_2[i]) >= hashikko:
        #if list_paf_start1[i] >= hashikko and list_paf_start2[i] >= hashikko: #9/10
            list_ox.append('x')
        elif ((length_1[i] - list_paf_end1[i]) / length_1[i]) >= hashikko and ((length_2[i] - list_paf_end2[i]) / length_2[i]) >= hashikko:
        #elif length_1[i] - list_paf_end1[i] >= hashikko and length_2[i] - list_paf_end2[i] >= hashikko: #9/10
            list_ox.append('x') 
        else:
            list_ox.append('o')

    elif list_paf_relative[i] == -1:
        if (list_paf_start1[i] / length_1[i]) >= hashikko and ((length_2[i] - list_paf_end2[i]) / length_2[i]) >= hashikko:
        #if list_paf_start1[i] >= hashikko and length_2[i] - list_paf_end2[i] >= hashikko: #9/10
            list_ox.append('x')
        elif ((length_1[i] - list_paf_end1[i]) / length_1[i]) >= hashikko and (list_paf_start2[i] / length_2[i]) >= hashikko:
        #elif length_1[i] - list_paf_end1[i] >= hashikko and list_paf_start2[i] >= hashikko: #9/10
            list_ox.append('x')
        else:
            list_ox.append('o')



'''（４）PAFデータの選択　＜フィルター１＞'''

list2_paf_name1 = []
list2_paf_start1 = []
list2_paf_end1 = []
list2_paf_relative = []
list2_paf_name2 = []
list2_paf_start2 = []
list2_paf_end2 = []
list2_paf_matchBase = []
list2_paf_alnLen = []
list2_paf_seqDiv = []
length2_1 = []
length2_2 = []

for i in range(len(list_ox)):
    if list_ox[i] == 'o':
        list2_paf_name1.append(list_paf_name1[i])
        list2_paf_start1.append(list_paf_start1[i])
        list2_paf_end1.append(list_paf_end1[i])
        list2_paf_relative.append(list_paf_relative[i])
        list2_paf_name2.append(list_paf_name2[i])
        list2_paf_start2.append(list_paf_start2[i])
        list2_paf_end2.append(list_paf_end2[i])
        list2_paf_matchBase.append(list_paf_matchBase[i])
        list2_paf_alnLen.append(list_paf_alnLen[i])
        list2_paf_seqDiv.append(list_paf_seqDiv[i])
        length2_1.append(length_1[i])
        length2_2.append(length_2[i])
        
#list_paf系をlist2_paf系で書き換え
list_paf_name1 = list2_paf_name1
list_paf_start1 = list2_paf_start1
list_paf_end1 = list2_paf_end1
list_paf_relative = list2_paf_relative
list_paf_name2 = list2_paf_name2
list_paf_start2 = list2_paf_start2
list_paf_end2 = list2_paf_end2
list_paf_matchBase = list2_paf_matchBase
list_paf_alnLen = list2_paf_alnLen
list_paf_seqDiv = list2_paf_seqDiv
length_1 = length2_1
length_2 = length2_2



'''（５）２箇所以上／２回以上アラインメントされるリードペアの除外　＜フィルター２＞'''
#Hollow indide cases are included here.

#dic_pairs: shows how many times each pair of reads appears in PAF.
dic_pairs = {}
for i in range(len(list_paf_name1)):
    k1 = fq_name.index(list_paf_name1[i])
    k2 = fq_name.index(list_paf_name2[i])

    if k1 == k2:    #Send off with no argument.
        dic_pairs[(list_paf_name1[i], list_paf_name2[i])] = 2
        
    elif k1 > k2:   #Send off if more than one.
        pair_name = (list_paf_name2[i], list_paf_name1[i])
        if not pair_name in dic_pairs:
            dic_pairs[pair_name] = 1
        else:
            dic_pairs[pair_name] += 1

    elif k1 < k2:   #Send off if more than one.
        pair_name = (list_paf_name1[i], list_paf_name2[i])
        if not pair_name in dic_pairs:
            dic_pairs[pair_name] = 1
        else:
            dic_pairs[pair_name] += 1

set_pair_twice = set({})
for key in dic_pairs:
    if dic_pairs[key] >= 2:
        set_pair_twice |= {key}

#Get numbers (in PAF) of strange alignments, and omit the numbers from PAF系.
rem_pos = []
for i in range(len(list_paf_name1)):
    if (list_paf_name1[i], list_paf_name2[i]) in set_pair_twice or (list_paf_name2[i], list_paf_name1[i]) in set_pair_twice:
        rem_pos.append(i)

paf_name1 = [list_paf_name1[i] for i in range(len(list_paf_name1)) if i not in rem_pos]
paf_start1 = [list_paf_start1[i] for i in range(len(list_paf_start1)) if i not in rem_pos]
paf_end1 = [list_paf_end1[i] for i in range(len(list_paf_end1)) if i not in rem_pos]
paf_name2 = [list_paf_name2[i] for i in range(len(list_paf_name2)) if i not in rem_pos]
paf_start2 = [list_paf_start2[i] for i in range(len(list_paf_start2)) if i not in rem_pos]
paf_end2 = [list_paf_end2[i] for i in range(len(list_paf_end2)) if i not in rem_pos]
paf_relative = [list_paf_relative[i] for i in range(len(list_paf_relative)) if i not in rem_pos]
paf_matchBase = [list_paf_matchBase[i] for i in range(len(list_paf_matchBase)) if i not in rem_pos]
paf_alnLen = [list_paf_alnLen[i] for i in range(len(list_paf_alnLen)) if i not in rem_pos]
paf_seqDiv = [list_paf_seqDiv[i] for i in range(len(list_paf_seqDiv)) if i not in rem_pos]
length_1 = [length_1[i] for i in range(len(length_1)) if i not in rem_pos]
length_2 = [length_2[i] for i in range(len(length_2)) if i not in rem_pos]

numPaf = len(paf_name1)     #Number of renewed PAF data.





#12/11 追加 アラインメント長順に並び替えてアラインメント長の下位 X ％は使用しないようにする
X = float(sys.argv[7])

lenAln_l = [(i, paf_end1[i] - paf_start1[i]) for i in range(numPaf)] 
lenAln_sorted = sorted(lenAln_l, key=lambda x: x[1])
lower_l = []
for i in range(round(numPaf * X)):
    lower_l.append(lenAln_sorted[i][0])

print("下位", X*100, "%点アラインメント長フィルタリング", lenAln_sorted[round(numPaf * X)-1])

paf_name1 = [paf_name1[i] for i in range(numPaf) if i not in lower_l]
paf_start1 = [paf_start1[i] for i in range(numPaf) if i not in lower_l]
paf_end1 = [paf_end1[i] for i in range(numPaf) if i not in lower_l]
paf_name2 = [paf_name2[i] for i in range(numPaf) if i not in lower_l]
paf_start2 = [paf_start2[i] for i in range(numPaf) if i not in lower_l]
paf_end2 = [paf_end2[i] for i in range(numPaf) if i not in lower_l]
paf_relative = [paf_relative[i] for i in range(numPaf) if i not in lower_l]
paf_matchBase = [paf_matchBase[i] for i in range(numPaf) if i not in lower_l]
paf_alnLen = [paf_alnLen[i] for i in range(numPaf) if i not in lower_l]
paf_seqDiv = [paf_seqDiv[i] for i in range(numPaf) if i not in lower_l]
length_1 = [length_1[i] for i in range(numPaf) if i not in lower_l]
length_2 = [length_2[i] for i in range(numPaf) if i not in lower_l]


numPaf = len(paf_name1)     # Renew number of PAF data.




# 2021 05 18 minimap2のava alignment結果=pafのペアワイズアラインメントの信頼度の高いものだけを使用する.
# 関数定義
# match: 一致塩基数 --> default >= 500
# div: sequence divergence --> default < 0.1
def Cleansing_by_MatchBase_SeqDiv(match, div, paf_data):
	idx = []
	for i in range(len(paf_data[0])):
		if int(paf_data[9][i]) < int(match):		# match base
			idx.append(i)
		if float(paf_data[11][i]) >= float(div):	# sequence divergence
			idx.append(i)
	
	new_paf = np.delete(paf_data, idx, axis=1)
	return new_paf




######ここまではデータの成型######




# ここまでにフィルターされたpafデータを行列にして格納
# PAF_data
# [ ["read1-1", "read2-1", "raed3-1"],
#   [100,       500,       200      ],
#   [300,       1000,      1500     ],
#   ["read1-2", "read2-2", "read3-2"],
#   [900,       400,       0        ],
#   [1100,      900,       1300     ],
#   [1,         -1,        1        ],
#   [3000,      2000,      5000     ],
#   [2000,      1200,      1400     ],
#   [length_1],
#   [length_2],
#   [paf_matchBase],
#   [paf_alnLen],
#   [paf_seqDiv]		    ]
PAF_data = np.array([paf_name1, paf_start1, paf_end1, paf_name2, paf_start2, paf_end2, paf_relative, length_1, length_2, paf_matchBase, paf_alnLen, paf_seqDiv])


PAF_data = Cleansing_by_MatchBase_SeqDiv(500, 0.1, PAF_data)
numPaf = len(PAF_data[0])


##### 関数定義 #####

def IdxSortPaf_by_aln_len(paf_data):
	aln_len_l = np.array([int(paf_data[2][i])-int(paf_data[1][i]) for i in range(len(paf_data[0]))])
	paf_data_sort = paf_data[:, aln_len_l.argsort()]
	paf_data_sort = np.flip(paf_data_sort, axis=1)
	return paf_data_sort


######ここからは重ね合わせ######
'''（６）+ - を揃える'''

def beginWith(l, paf_data):
    plmi = [0 for i in range(numFq)]
    plmi[l] = 1

    for k in range(numPaf):
        if paf_data[0][k] == fq_name[l]:
            i = fq_name.index(paf_data[3][k])
            plmi[i] = int(paf_data[6][k])

        elif paf_data[3][k] == fq_name[l]:
            i = fq_name.index(paf_data[0][k])
            plmi[i] = int(paf_data[6][k])

    return plmi


def plusMinusOn(l, paf_data):
    lstup = []
    for i in range(numFq):
        plmi_i = beginWith(i, paf_data)
        lstup.append(plmi_i)

    plmi_l = lstup[l]
    plusMinusB = [0 for i in range(numFq)]
    plusMinusA = copy.copy(plmi_l)

    while plusMinusB != plusMinusA:
        plusMinusB = copy.copy(plusMinusA)
        plusMinusA = copy.copy(plusMinusB)

        for i in range(numFq):
            plmi_i = lstup[i]
            cnt = 0
            relation_i = 0
            for j in range(numFq):
                if plmi_i[j] != 0 and plusMinusA[j] != 0:
                    cnt += 1
                    relation_i = int(plmi_i[j] / plusMinusA[j])

            if cnt != 0:
                for m in range(numFq):
                    if plusMinusA[m] == 0 and plmi_i[m] != 0:
                        plusMinusA[m] = plmi_i[m] * relation_i

    return plusMinusA



PAF_data_sort = IdxSortPaf_by_aln_len(PAF_data)		# PAF_data をアラインメント長順に並び替え(降順)

paf_name1 = PAF_data_sort[0]
paf_start1 = np.array(PAF_data_sort[1], dtype=int)
paf_end1 = np.array(PAF_data_sort[2], dtype=int)
paf_name2 = PAF_data_sort[3]
paf_start2 = np.array(PAF_data_sort[4], dtype=int)
paf_end2 = np.array(PAF_data_sort[5], dtype=int)
paf_relative = np.array(PAF_data_sort[6], dtype=int)
length_1 = np.array(PAF_data_sort[7], dtype=int)
length_2 = np.array(PAF_data_sort[8], dtype=int)
#paf_matchBase = np.array(PAF_data_sort[9], dtype=int)
#paf_alnLen = np.array(PAF_data_sort[10], dtype=int)
#paf_seqDiv = np.array(PAF_data_sort[11], dtype=int)


#Ranking by length ---Descending order---
fq_length = [len(i) for i in fq_seq]
indexes_order = [i[0] for i in sorted(enumerate(fq_length), key=lambda x: x[1])]
indexes_order = indexes_order[::-1]
print('indexes_order', indexes_order)

# 長いリードを ± 基準にする.

plmi_l_10 = []
for i in range(10):     #長いリード上から10本
    plmi_l_10.append(plusMinusOn(indexes_order[i], PAF_data_sort))

cnt_same_plmi = [0 for i in range(10)]
for i in range(10):
    p_i = plmi_l_10[i]
    p_i_rev = [-p_i[k] for k in range(len(p_i))]
    for j in range(10):
        p_j = plmi_l_10[j]
        if p_i == p_j or p_i_rev == p_j:
            cnt_same_plmi[i] += 1

max_idx = 0
max_num = 0
for i in range(10):
    if max_num < cnt_same_plmi[i]:
        max_idx = i
        max_num = cnt_same_plmi[i]

base = indexes_order[max_idx]
print('base+-', base)


plusMinus = plusMinusOn(base, PAF_data_sort)
print("plusMinus\n", plusMinus)
#ここまで・plusMinus作る過程



'''（７）重ね合わせ'''

#The parameter used below, hashikko has been already defined above, in the 1st filter.
def yoke(k):
    lst = [paf_relative[k]]

    i = fq_name.index(paf_name1[k])
    j = fq_name.index(paf_name2[k])
    lst.append((i, j))

    range1 = (0, length_1[k])

    if lst[0] == 1:
        range2 = (paf_start1[k] - paf_start2[k], paf_start1[k] - paf_start2[k] + length_2[k])
        lst.append([range1, range2])

    elif lst[0] == -1:

        if ( (length_1[k] - paf_end1[k]) / length_1[k] ) < hashikko:
            range2 = (paf_start1[k] - (length_2[k] - paf_end2[k]), paf_start1[k] + paf_end2[k])
            lst.append([range1, range2])

        elif ( paf_start2[k] / length_2[k] ) < hashikko:
            range2 = (paf_end1[k] - (length_2[k] - paf_start2[k]), paf_end1[k] + paf_start2[k])
            lst.append([range1, range2])

    return lst      #[ ±1, (i,j), [(0,**), (**,**)] ] 


def combine(l, rtn1, rtn2):
    trio = [( 0, len(fq_seq[l]) )]

    if rtn1[1][0] == l and rtn2[1][0] == l:     ##場合分け 1,5,9,13
        i_range = rtn1[2][1]
        j_range = rtn2[2][1]
        trio.append(i_range)
        trio.append(j_range)
        return trio

    elif rtn1[0] == 1 and rtn2[0] == 1:         ##場合分け 2,3,4

        if rtn1[1][0] == l and rtn2[1][1] == l:     #(場合分け 2)
            i_range = rtn1[2][1]
            j_range = (-rtn2[2][1][0], -rtn2[2][1][0] + rtn2[2][0][1])
            trio.append(i_range)
            trio.append(j_range)
            return trio

        elif rtn1[1][1] == l and rtn2[1][0] == l:   #(場合分け 3)
            i_range = (-rtn1[2][1][0], -rtn1[2][1][0] + rtn1[2][0][1])
            j_range = rtn2[2][1]
            trio.append(i_range)
            trio.append(j_range)
            return trio

        elif rtn1[1][1] == l and rtn2[1][1] == l:   #(場合分け 4)
            i_range = (-rtn1[2][1][0], -rtn1[2][1][0] + rtn1[2][0][1])
            j_range = (-rtn2[2][1][0], -rtn2[2][1][0] + rtn2[2][0][1])
            trio.append(i_range)
            trio.append(j_range)
            return trio

    elif rtn1[0] == 1 and rtn2[0] == -1:        ##場合分け 6,7,8

        if rtn1[1][0] == l and rtn2[1][1] == l:     #(場合分け 6)
            i_range = rtn1[2][1]
            j_range = (rtn2[2][1][1] - rtn2[2][0][1], rtn2[2][1][1])
            trio.append(i_range)
            trio.append(j_range)
            return trio

        elif rtn1[1][1] == l and rtn2[1][0] == l:   #(場合分け 7)
            i_range = (-rtn1[2][1][0], -rtn1[2][1][0] + rtn1[2][0][1])
            j_range = rtn2[2][1]
            trio.append(i_range)
            trio.append(j_range)
            return trio

        elif rtn1[1][1] == l and rtn2[1][1] == l:   #(場合分け 8)
            i_range = (-rtn1[2][1][0], -rtn1[2][1][0] + rtn1[2][0][1])
            j_range = (-rtn2[2][0][1] + rtn2[2][1][1], rtn2[2][1][1])
            trio.append(i_range)
            trio.append(j_range)
            return trio

    elif rtn1[0] == -1 and rtn2[0] == 1:        ##場合分け 10,11,12

        if rtn1[1][0] == l and rtn2[1][1] == l:     #(場合分け 10)
            i_range = rtn1[2][1]
            j_range = (-rtn2[2][1][0], -rtn2[2][1][0] + rtn2[2][0][1])
            trio.append(i_range)
            trio.append(j_range)
            return trio

        elif rtn1[1][1] == l and rtn2[1][0] == l:   #(場合分け 11)
            i_range = (rtn1[2][1][1] - rtn1[2][0][1], rtn1[2][1][1])
            j_range = rtn2[2][1]
            trio.append(i_range)
            trio.append(j_range)
            return trio

        elif rtn1[1][1] == l and rtn2[1][1] == l:   #(場合分け 12)
            i_range = (rtn1[2][1][1] - rtn1[2][0][1], rtn1[2][1][1])
            j_range = (-rtn2[2][1][0], -rtn2[2][1][0] + rtn2[2][0][1])
            trio.append(i_range)
            trio.append(j_range)
            return trio

    elif rtn1[0] == -1 and rtn2[0] == -1:       ##場合分け 14,15,16

        if rtn1[1][0] == l and rtn2[1][1] == l:     #(場合分け 14)
            i_range = rtn1[2][1]
            j_range = (rtn2[2][1][1] - rtn2[2][0][1], rtn2[2][1][1])
            trio.append(i_range)
            trio.append(j_range)
            return trio

        elif rtn1[1][1] == l and rtn2[1][0] == l:   #(場合分け 15)
            i_range = (rtn1[2][1][1] - rtn1[2][0][1], rtn1[2][1][1])
            j_range = rtn2[2][1]
            trio.append(i_range)
            trio.append(j_range)
            return trio

        elif rtn1[1][1] == l and rtn2[1][1] == l:   #(場合分け 16)
            i_range = (rtn1[2][1][1] - rtn1[2][0][1], rtn1[2][1][1])
            j_range = (rtn2[2][1][1] - rtn2[2][0][1], rtn2[2][1][1])
            trio.append(i_range)
            trio.append(j_range)
            return trio


def pileOn(l):
    pool = []
    for k in range(numPaf):
        if fq_name[l] == paf_name1[k] or fq_name[l] == paf_name2[k]:
            pool.append(yoke(k))

    pile = [() for tup in range(numFq)]
    for m in range(len(pool)):
        if pool[m][1][0] == l:
            i = pool[m][1][1]
        elif pool[m][1][1] == l:
            i = pool[m][1][0]

        for n in range(m+1, len(pool)):
            if pool[n][1][0] == l:
                j = pool[n][1][1]
            elif pool[n][1][1] == l:
                j = pool[n][1][0]

            trio = combine(l, pool[m], pool[n])
            if pile[l] == ():
                pile[l] = trio[0]
            if pile[i] == ():
                pile[i] = trio[1]
            if pile[j] == ():
                pile[j] = trio[2]

    return pile


def flip(pile, l1, l2):
    #フリップ
    if plusMinus[l1] != plusMinus[l2]:
        for l in range(numFq):
            if pile[l] != ():
                pile[l] = list(pile[l])

                tmp = pile[l][0]
                pile[l][0] = pile[l][1]
                pile[l][1] = tmp
                pile[l][0] *= -1
                pile[l][1] *= -1

                pile[l] = tuple(pile[l])

    return pile


def collate(pile1, pile2, l1, l2):
    pile2 = flip(pile2, l1, l2)

    cnt = 0
    ij = (0, 0)
    for i in range(numFq):
        for j in range(numFq):
            if i < j:
                if pile1[i] != () and pile1[j] != () and pile2[i] != () and pile2[j] != ():
                    cnt += 1
                    ij = (i, j)

    if cnt != 0:
        i = ij[0]
        diff = pile2[i][0] - pile1[i][0]

        for l in range(numFq):
            if pile1[l] != () and pile2[l] !=():
                continue
            elif pile1[l] != () and pile2[l] == ():
                continue
            elif pile1[l] == () and pile2[l] != ():
                pile1[l] = (pile2[l][0] - diff,  pile2[l][1] - diff)
            elif pile1[l] == () and pile2[l] == ():
                continue

    pile = pile1

    return pile 


def placementBy(l1, l2):  #Based on the longest reads, l1, l2, call this function.

    pile1 = pileOn(l1)  #(l_rndm1)
    pile2 = pileOn(l2)  #(l_rndm2)
    plcmnt_tmp = collate(pile1, pile2, l1, l2)
    plcmnt = []
    while plcmnt != plcmnt_tmp:
        plcmnt = plcmnt_tmp
        for l in range(numFq):
            #plcmnt_tmp = collate(plcmnt_tmp, pileOn(l), l_rndm1, l)
            plcmnt_tmp = collate(plcmnt_tmp, pileOn(l), l1, l)

    return plcmnt


#Ranking by length ---Descending order---
fq_length = [len(i) for i in fq_seq]
indexes_order = [i[0] for i in sorted(enumerate(fq_length), key=lambda x: x[1])]
indexes_order = indexes_order[::-1]
#print('ranking by length', len(fq_seq[indexes_order[0]]), len(fq_seq[indexes_order[1]]))


# Added on June 27, after meeting with Prof., --> use the longest reads as the basis.

nonZero = []
for i in range(numFq):
    if plusMinus[i] != 0:
        nonZero.append(i)

print('plusMinus  nonZero数:', len(nonZero))

# 一番基礎になる長いリード２本を選択（plusMinusが0なら次に長いリードを選択）
ind1, ind2 = indexes_order[0], indexes_order[1]
print('ind1', ind1, 'ind2', ind2)
ind_cnt = 2
while plusMinus[ind1] == 0:
    ind1 = indexes_order[ind_cnt]
    ind_cnt += 1
while plusMinus[ind2] == 0:
    ind2 = indexes_order[ind_cnt]
    ind_cnt += 1

print(ind_cnt-2, '回基準が取り直された。基準は', ind_cnt-2, '番目')
placement = placementBy(ind1, ind2)

# もし選択された２本がなすplacement（アセンブリ）が空になる場合は、さらに次のリードを選択
while placement == [() for i in range(numFq)]:
    ind1_ = indexes_order[ind_cnt]
    placement = placementBy(ind1_, ind2)
    if placement != [() for i in range(numFq)]:
        ind1 = ind1_
        break
    if placement == [() for i in range(numFq)]:
        ind2_ = indexes_order[ind_cnt]
        placement = placementBy(ind1, ind2_)
        if placement != [() for i in range(numFq)]:
            ind2 = ind2_
            break

    ind_cnt += 2
    ind1 += 2
    ind2 += 2

    if placement == [() for i in range(numFq)]:
        placement = placementBy(ind1, ind2)

print("placement\n", placement)

######ここまで重ね合わせ######


# 10/2 追加 アセンブリの途中で途切れる箇所が存在するか確認





'''plusMinusを参照して（基底に対して）-鎖は逆転・相補をとる'''
# Reverse complement
for l in range(numFq):
    if plusMinus[l] * plusMinus[ind1] == -1:
        fq_seq[l] = fq_seq[l][::-1]
        fq_seq[l] = fq_seq[l].replace('A', 'X')
        fq_seq[l] = fq_seq[l].replace('T', 'A')
        fq_seq[l] = fq_seq[l].replace('X', 'T')
        fq_seq[l] = fq_seq[l].replace('C', 'X')
        fq_seq[l] = fq_seq[l].replace('G', 'C')
        fq_seq[l] = fq_seq[l].replace('X', 'G')



'''（８）可視化'''
import matplotlib.pyplot as plt

nonZero = []
for i in range(numFq):
#for i in range(0, 40):
    if placement[i] != ():
        nonZero.append(i)
        x = np.arange(placement[i][0], placement[i][1], 1)
        y = x ** 0 + i - 1
        plt.plot(x, y)


processTime = time.time() - startTime
print('Run time:', processTime, '[sec]')

#figure name
path_figure = sys.argv[8]
fastq_path = sys.argv[1]
#m = re.search(r'(NA|GM).*\.sam', fastq_path)
#pos_name = m.group()
pos_name = os.path.basename(fastq_path)
if "_hap" in pos_name:
	pos_name = pos_name[:-18]
else:
	pos_name = pos_name[:-10]

plt.savefig(f'{path_figure}/{pos_name}.png')



'''（９）Process to make FASTA files to apply to MAFFT'''
#Use both plusMinus and placement below.

blockSize = int(sys.argv[4])
stepSize = int(sys.argv[5])
rmv = float(sys.argv[6])

vacancies = []      #上のnonZero使うなら必要なし！
for i in range(numFq):
    if placement[i] == ():
        vacancies.append(i)

maximum = max(placement[i][1] for i in range(numFq) if not i in vacancies)
minimum = min(placement[i][0] for i in range(numFq) if not i in vacancies)
zenRange = maximum - minimum
print('min:', minimum, 'max:', maximum)


#リード平均長・分散の計算 9/7 追加


#Make FASTA files throughout the whole range of the reads' layout.
#Divide by block once a stepsize, then output truncated sub-sequences.

# 10/1 ブロック内リード１本の時でもfastaファイル生成

cnt = 0
#cnt_i = 0
cnt_reads_l = []
#block_num_eaten = []

path_file = sys.argv[8]
fastq_path = sys.argv[1]
#m = re.search(r'(NA|GM).*\.sam', fastq_path)
#pos_name = m.group()
pos_name = os.path.basename(fastq_path)
if "_hap" in pos_name:
	pos_name = pos_name[:-18]
else:
	pos_name = pos_name[:-10]

while minimum + stepSize*cnt + blockSize-1 < maximum:
    cnt_reads = 0

    file_name = f'{pos_name}_{minimum+stepSize*cnt}-{minimum+stepSize*cnt+blockSize-1}_out2.fa'
    s = ''

    blockRange = (minimum + stepSize*cnt, minimum + stepSize*cnt + blockSize-1)
    for i in range(numFq):
        if not i in vacancies:
            if blockRange[0] <= placement[i][1] and placement[i][0] <= blockRange[1]:   #かかる条件
                if placement[i][0] <= blockRange[0]:            # リード始まりがブロック左端より左
                    if blockRange[1] <= placement[i][1]:                # (1) リード終わりがブロック右端より右
                        start_i = blockRange[0] - placement[i][0]   # (0-origin index)
                        end_i = start_i + blockSize - 1
                        
                    else:                                               # (2) リード終わりがブロック右端より左
                        start_i = blockRange[0] - placement[i][0]
                        end_i = placement[i][1] - placement[i][0] - 1

                                                                # リード始まりがブロック左端より右
                elif blockRange[1] <= placement[i][1]:                  # (3) リード終わりがブロック右端より右
                    start_i = 0
                    end_i = blockRange[1] - placement[i][0] - 1

                else:                                                   # (4) リード終わりがブロック右端より左
                    start_i = 0
                    end_i = placement[i][1] - placement[i][0]

                readSize = end_i - start_i + 1
                if (readSize / blockSize) >= rmv:
                    s += f'>{fq_name[i]}\t{(start_i, end_i)}\t{plusMinus[i]}\n'
                    seq = fq_seq[i]
                    s += f'{seq[start_i:end_i+1]}\n'

                    cnt_reads += 1

    if cnt_reads == 1:  #次ステップでmafftのエラーを避けるためにリード数１本のブロックは名前を変更しておく(subconsensus終了状態にする) 
        file_name = f'{pos_name}_{minimum+stepSize*cnt}-{minimum+stepSize*cnt+blockSize-1}_out4.txt'
        name = f'>{pos_name}_{minimum+stepSize*cnt}-{minimum+stepSize*cnt+blockSize-1}\t1'    #'1'はリード数
        m_seq = re.search(r'\n(A|T|G|C)*\n', s)
        seq = m_seq.group()[1:-1]
        seq = seq.lower()
        s = ''
        s += name + '\n'
        s += seq + '\n'     #一応改行


    #if cnt_reads <= int(sys.argv[7]):     # sys.argv[7] - ブロックに含まれるリード数下限
    #if cnt_reads < b:                     # 9/7変更点 sys.argv[7] を本数ではなく、割合に変更?
        #block_num_eaten.append(cnt)
        #cnt += 1
        #continue

    #print('file_name', file_name)
    
    with open(f'{path_file}/{file_name}', mode='w') as f:
        f.write(s)

    cnt += 1
    #cnt_i += 1

    cnt_reads_l.append(cnt_reads)


