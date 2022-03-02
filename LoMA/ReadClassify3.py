"""
コンセンサス配列に対してリードをアラインメントした。
このスクリプトでアラインメントCIGARから「長いクリップ」と「長い挿入・欠失」を抽出する。
"""

print("---ReadClassify.py---")

import re
import sys
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.special import comb

# sys.argv[1] - SAM file
# sys.argv[2] - FLOAT 端％（default=0.1）
# sys.argv[3] - INT 端長さ（default=500）
# sys.argv[4] - INT 挿入・欠失　検出下限の塩基長（default=100）
# sys.argv[5] - INT 異常位置の変動幅（default=100）
# sys.argv[6] - INT 異常位置とみなし検出する下限リード本数（default=8）
# sys.argv[7] - FASTQ file
# sys.argv[8] - png 保存ディレクトリ
# sys.argv[9] - NA18943_chr1_123456-198765.sam.fastq.zi_hapj 保存ディレクトリ

#グローバル変数
hashi_per = float(sys.argv[2])
hashi_len = int(sys.argv[3])
indel_len = int(sys.argv[4])


'''____________関数定義___________'''

def CIGARist(cig):
    cig_MID = []
    for i in range(len(cig)):
        if cig[i] in ['M', 'I', 'D', 'S', 'H']:
            cig_MID.append(cig[i])
    cig_num = re.split(r'[MIDSH]', cig)
    cig_num = cig_num[:-1]
    cig_num = [int(c) for c in cig_num]

    return (cig_num, cig_MID)   #input: CIGAR配列 --> output: ([CIGAR長], [CIGAR記号])



def clipFinder(sam_i, length):
    clip_pos = {'H':[0, 0], 'S':[0, 0]}     #戻り値: 位置
    clip_len_per = [0, 0]
    clip_len_LR = [0, 0]    # あとで付け足し: 左端・右端合わせた長さだけでなく、その側の長さも見る

    if sam_i[3][0] == 'H' or sam_i[3][0] == 'S':
        clip_len_per[0] += sam_i[2][0]
        clip_len_LR[0] += sam_i[2][0]               # 付け足し
    if sam_i[3][-1] == 'H' or sam_i[3][-1] == 'S':
        clip_len_per[0] += sam_i[2][-1]
        clip_len_LR[1] += sam_i[2][-1]              # 付け足し

    clip_len_per[1] = clip_len_per[0] / length
    #if clip_len_per[0] > hashi_len and clip_len_per[1] > hashi_per:
    # 異常なclip --> 位置を記録
    if clip_len_LR[0] > hashi_len and clip_len_LR[0] / length > hashi_per:
        if sam_i[3][0] == 'H':
            clip_pos['H'][0] = sam_i[0]
        elif sam_i[3][0] == 'S':
            clip_pos['S'][0] = sam_i[0]
    if clip_len_LR[1] > hashi_len and clip_len_LR[1] / length > hashi_per:
        if sam_i[3][-1] ==  'H':
            clip_pos['H'][1] = sam_i[0] + sum(sam_i[2][j] for j in range(len(sam_i[2])) if sam_i[3][j] in ('M', 'D'))
        elif sam_i[3][-1] == 'S':
            clip_pos['S'][1] = sam_i[0] + sum(sam_i[2][j] for j in range(len(sam_i[2])) if sam_i[3][j] in ('M', 'D'))

    return clip_pos



def indelFinder(sam_i):
    # indel_len 変数 bp以上の挿入欠失を検出
    cig_num = sam_i[2]
    cig_MID = sam_i[3]
    indel_pos = sam_i[0]
    indel_pos_d = {'ins':[], 'del':[]}

    for i in range(len(cig_MID)): 
        if cig_MID[i] == 'I':
            if cig_num[i] >= indel_len:
                indel_pos_d['ins'].append(indel_pos)

        elif cig_MID[i] == 'D':
            if cig_num[i] >= indel_len:
                indel_pos_d['del'].append(indel_pos)

        if cig_MID[i] != 'S' and cig_MID[i] != 'H':
            if cig_MID[i] == 'M' or cig_MID[i] == 'D':  #'M' or 'D' ならリファレンスの位置が進む
                indel_pos += cig_num[i]

    return indel_pos_d



def StrangerWhereabout(clip_, indel_, low_, length_):
    #戻り値：異常位置のリスト　(存在しない場合は空のリストを戻す)
    pos_cnt = clip_ + indel_
    strangeWhere_l = []
    for i in range(len(pos_cnt)):
        if low_ <= pos_cnt[i]:
            strangeWhere_l.append(i)
    if strangeWhere_l == []:
        return []
    else:
        #はじめと最後の3000塩基は切り落とす
        strangeWhere_l = [i for i in strangeWhere_l if 3000 <= i <= (length_ - 3000)]
        if strangeWhere_l == []:
            return []
        else:
            #同じ異常はひとまとめにする
            cnt = []
            pnt = strangeWhere_l[0]
            for i in range(1, len(strangeWhere_l)):
                if not pnt <= strangeWhere_l[i] < pnt + 100:
                    cnt.append(i-1)
                pnt = strangeWhere_l[i]
        
            strangeWhere_integrated = []
            pnt = 0
            for i in cnt:
                m = (strangeWhere_l[pnt] + strangeWhere_l[i]) // 2
                strangeWhere_integrated.append(m)
                pnt = i + 1
            m = (strangeWhere_l[pnt] + strangeWhere_l[-1]) // 2
            strangeWhere_integrated.append(m)
            return strangeWhere_integrated



def Range(sam_):
    #戻り値：リードの範囲, clipの部分はそのままアラインメントに含めてカウント
    range_dict = {}
    for i in sam_:
        if sam_[i][3][0] == 'H' or sam_[i][3][0] == 'S':
            range_i = [sam_[i][0] - sam_[i][2][0]]
        else:
            range_i = [sam_[i][0]]
        ending_pos_i = sam_[i][0]
        for j in range(1, len(sam_[i][3])):
            if sam_[i][3][j] in ['M', 'D', 'H', 'S']:
                ending_pos_i += sam_[i][2][j]
        range_i.append(ending_pos_i)
        range_dict[i] = range_i
    return range_dict



def StrangerDetector(range_, pos_):
    #戻り値：異常なリード名
    stranger = []
    for i in range_:
        if range_[i][0] <= pos_ <= range_[i][1]:
            stranger.append(i)
    return stranger



def classifier(pos_i, z_i, clip_, indel_):
    # 戻り値：異常リードのリスト
    z_hap1 = []
    # ±500bp はその異常とみなす
    for read in clip_:
        # clip で z_i_hap1
        if clip_[read]['H'][0] != 0:
            if clip_[read]['H'][0] - 500 <= pos_i <= clip_[read]['H'][0] + 500:
                z_hap1.append(read)
        if clip_[read]['H'][1] != 0:
            if clip_[read]['H'][1] - 500 <= pos_i <= clip_[read]['H'][1] + 500:
                z_hap1.append(read)
        if clip_[read]['S'][0] != 0:
            if clip_[read]['S'][0] - 500 <= pos_i <= clip_[read]['S'][0] + 500:
                z_hap1.append(read)
        if clip_[read]['S'][1] != 0:
            if clip_[read]['S'][1] - 500 <= pos_i <= clip_[read]['S'][1] + 500:
                z_hap1.append(read)

    for read in indel_:
        # indel で z_i_hap1
        for j in indel_[read]['ins']:
            if j - 500 <= pos_i <= j + 500:
                z_hap1.append(read)
        for j in indel_[read]['del']:
            if j - 500 <= pos_i <= j + 500:
                z_hap1.append(read)

    return set(z_hap1)


# 2021 / 05 / 24
def ConsecutiveIntegrate(integ, consec):
	if integ == []:
		return []
	tmp = [1]
	for i in range(1, len(integ)):
		if integ[i]-consec < integ[i-1]:
			tmp.append(1)
		else:
			tmp.append(0)
	#print(tmp)
	pos_tmp_l = []
	pos_tmp = [integ[0]]
	for i in range(1, len(integ)):
		if tmp[i] == 1:
			pos_tmp.append(integ[i])
		else:
			pos_tmp_l.append(pos_tmp)
			pos_tmp = [integ[i]]
	pos_tmp_l.append(pos_tmp)
	pos_integrated = []
	for i in range(len(pos_tmp_l)):
		if len(pos_tmp_l[i]) != 1:
			average_i = round(sum(pos_tmp_l[i])/len(pos_tmp_l[i]))
		else:
			average_i = pos_tmp_l[i][0]
		pos_integrated.append(average_i)
	return pos_integrated


def Binomial(coverage, pr):
	cdf = []
	P = 0.5 ** coverage
	for i in range(0, coverage+1):
		if i == 0:
			cdf.append(P)
		else:
			p = comb(coverage, i) * P
			cdf.append(cdf[i-1]+p)
	#print("CDF", cdf)
	
	k = 0
	for i in range(len(cdf)):
		if cdf[i] >= pr:
			k = i
			break
	return k			# coverage本のうち確率prを超える最低本数を二項分布から計算して返す.

def Binomial2(coverage, n):
	mean = coverage/2
	SD = np.sqrt(coverage)/2
	sigma = n * SD
	sigma_range = [round(mean-sigma, 2), round(mean+sigma, 2)]
	return sigma_range

'''_____________リード辞書づくり_____________'''

# 698452f5-6fbf-4c3a-9378-86c3c7d7e313    2048    NA18943_chr1_79077314-79117314_80       37854   1       2885H29M2D4M1I10M6I37M3D3M1I51M1I2M2D59M1I8M2D62M1D11M3681H     *       0       0       GTGGCTCACGCCTATAATCCCTGCACTTTGAGGCCCAAGGCAGGAAAGCAGGATCACCTACGGTCAGGAGTTCAAGACCAGCCTGGCCATGGGCGAAACCCCATCTCTACTAAAAATACAAAAAGTAGCCAGGTATGGCGGCATCGCCTGTAATCCCAGCTACTCGGGAGGCTGAGGTAAGAGAATCACTTGAACCCAGGAGGCTGGAGGTTGGTGAGCTAAGATTGCACCATTGCACTCCAGCCTGGGTGACAGAGTGGGACTCCCTATCAAAAAAAAAACAAAG  =.06:=<?BA/9;(BJILHF;:<B<>B@<+48=<?>(-53//2%&('''''.344*=;+69..'''=(%+2,*&)8.:1.--/)()')''&''$%-1.+,00,33**,=9:::;;9&5====;2+65::75:+))6:9/5*((*((%CB@Z]NPMN17==?967603C<-(&$3,4.,20//0-0*,:JHM;7>/@:7<66)*$$)'$'(&&$#$((*)&(.+**/5:8<B5258@LMMK>=>FGDLQA3>@8,0:;<;46+,9@:=<;@<GHFGEB843*+///%  NM:i:59 ms:i:238        AS:i:238        nn:i:0  tp:A:P  cm:i:5  s1:i:45 s2:i:44 de:f:0.1713     SA:Z:NA18943_chr1_79077314-79117314_80,26484,-,771S324M6D5757S,57,51;NA18943_chr1_79077314-79117314_80,14010,+,1753S308M4791S,60,52;NA18943_chr1_79077314-79117314_80,17099,+,2575S299M8D3978S,10,67;   rl:i:0


f = open(sys.argv[1])
sam_dict = {}   # key:リード名, value:リスト [アラインメント開始点, リード塩基配列, [cigar num], [cigar MIDSH]]
for line in f:
    if line[0] == '@':
        continue
    else:
        line_split = line.split('\t')
        name = line_split[0]
        if name in sam_dict:   #複数の位置にアラインメントされているリードは一番上のアラインメントを採用
            continue

        cigar_num_mid = CIGARist(line_split[5])
        
        if cigar_num_mid[0] == [] or cigar_num_mid[1] == []:
            continue

        sam_dict[name] = [int(line_split[3])]      #アラインメント開始点:コンセンサス配列における座標
        sam_dict[name].append(line_split[9])       #リード塩基配列:Hard clip配列は含まない, M/I/Sの和に等しい  
        sam_dict[name].append(cigar_num_mid[0])    #cigar number:リストで順番に格納
        sam_dict[name].append(cigar_num_mid[1])    #cigar MIDSH:リストで順番に格納, cigar numberと対応

        #print('cig 長さ', sum(sam_dict[name][2][i] for i in range(len(sam_dict[name][2])) if sam_dict[name][3][i] != 'D'), 'seq 長さ', len(sam_dict[name][1]))
f.close()

#リファレンス(コンセンサス配列)の長さをSAMのヘッダーから抽出
f = open(sys.argv[1])
len_ref = 0
for line in f:
    if line[0:3] == '@SQ':
        header_split = line.split('\t')
        for i in header_split:
            if i[0:2] == 'LN':
                m = re.search('LN:[0-9]*', i)
                len_ref = m.group()
                len_ref = len_ref[3:]
                if len_ref[-1] == '\n' or len_ref[-1] == '\t':
                    len_ref = len_ref[:-1]
                len_ref = int(len_ref)

f.close()
#print('sam_dict', sam_dict)



'''___________クリップ・挿入欠失の記録___________'''

#total_num_read = 0
#end2end_num_read = 0
clip_dict = {}      #key: リード名, value: {'H': [len, len], 'S': [len, len]}
indel_dict = {}     #key: リード名, value: {'ins': [pos, ...], 'del': [pos, ...]}

for i in sam_dict:
    #total_num_read += 1

    #リード本来の長さ
    len_read = sum(sam_dict[i][2][j] for j in range(len(sam_dict[i][3])) if sam_dict[i][3][j] in ['M', 'I', 'S', 'H'])

    #クリップ
    clip_i = clipFinder(sam_dict[i], len_read)
    clip_dict[i] = clip_i

    #挿入・欠失
    indel_i = indelFinder(sam_dict[i])
    indel_dict[i] = indel_i



'''_______________異常な位置を検知_______________'''

''' ここでは、異常な位置を「クリップの始／終が集積する位置」or「挿入・欠失が集積する位置」とする。
変動幅(fluctuation) <-- パラメータ、の範囲内の異常は同一とみなす。
一定数 <-- パラメータ、以上の本数のリードがあれば、その位置に異常があるとみなす。'''

fluc = int(sys.argv[5])
low_limit = int(sys.argv[6])

#print("clip_dict", clip_dict)
#print("indel_dict", indel_dict)

# 位置をリスト化
clip_pos_l = []
indel_pos_l = []
for name in sam_dict:
    # clip 追加
    if clip_dict[name]['H'][0] != 0:
        clip_pos_l.append(clip_dict[name]['H'][0])
    elif clip_dict[name]['S'][0] != 0:
        clip_pos_l.append(clip_dict[name]['S'][0])
    if clip_dict[name]['H'][1] != 0:
        clip_pos_l.append(clip_dict[name]['H'][1])
    elif clip_dict[name]['S'][1] != 0:
        clip_pos_l.append(clip_dict[name]['S'][1])

    # indel 追加
    for i in indel_dict[name]['ins']:
        indel_pos_l.append(i)
    for i in indel_dict[name]['del']:
        indel_pos_l.append(i)

# カウント
clip_cnt_l = []
indel_cnt_l = []
cnt = -1
while cnt * 1 < len_ref:
    cnt += 1
    cover_range = [cnt*1, cnt*1 + fluc]

    cnt_clip_in_range = 0
    for i in clip_pos_l:
        if cover_range[0] <= i < cover_range[1]:
            cnt_clip_in_range  += 1
    clip_cnt_l.append(cnt_clip_in_range)

    cnt_indel_in_range = 0
    for i in indel_pos_l:
        if cover_range[0] <= i < cover_range[1]:
            cnt_indel_in_range += 1
    indel_cnt_l.append(cnt_indel_in_range)


y = np.array(clip_cnt_l) + np.array(indel_cnt_l)
#print("y", list(y))
plt.plot(y)
#plt.show()
plt.savefig(sys.argv[8]+'/'+os.path.basename(sys.argv[7])+'.png')

clip_cnt_l = np.array(clip_cnt_l)       #各位置から+fluc塩基の範囲のclipの数
indel_cnt_l = np.array(indel_cnt_l)     #各位置から+fluc塩基の範囲のindelの数

strangeWhere_integrated = StrangerWhereabout(clip_cnt_l, indel_cnt_l, low_limit, len_ref)     #low_limit本以上の異常リードが集積している位置を特定しリスト化
# さらに統合
strangeWhere_integrated = ConsecutiveIntegrate(strangeWhere_integrated, 700)


print(sys.argv[7])
print('異常位置候補:', len(strangeWhere_integrated), 'コ', strangeWhere_integrated)



# 検出する最低値を定める
max_cnt = np.amax(y)
print("max_cnt", max_cnt)
average_cnt = np.average(y)
print("average_cnt", average_cnt)


'''______________平均coverageによる二項確率の検出下限______________'''

range_dict = Range(sam_dict)    #リードの被覆範囲を辞書化（クリップも含める）
all_read_name = [i for i in sam_dict]
# 5/24
print("全体の#read", len(all_read_name))
coverage = round(sum(range_dict[i][1]-range_dict[i][0] for i in range_dict)/len(y))
print("全体の平均coverage", coverage)


'''______________異常位置のリード______________'''

range_dict = Range(sam_dict)    #リードの被覆範囲を辞書化（クリップも含める）
all_read_name = [i for i in sam_dict]




Z_hap1 = []     # 異常位置 Z1, Z2, ... , Zn (n個の異常位置) のclip or indel リード名を集める
Z_hap2 = []     # 異常位置 Z1, Z2, ... , Zn (n個の異常位置) にあって clip, indel のないリード名を集める
notZ = []       # Z_hap1 にも Z_hap2 にもないリード名を集める
for i in strangeWhere_integrated:
    z_i = StrangerDetector(range_dict, i)    #Ziを被覆している全リード
    z_i_hap1 = list(classifier(i, z_i, clip_dict, indel_dict))
    z_i_hap2 = list(set(z_i) - set(z_i_hap1))
    notz_i = list(set(all_read_name) - set(z_i))

    print("z_i_hap1", len(z_i_hap1), z_i_hap1)
    print("z_i_hap2", len(z_i_hap2), z_i_hap2)
    
    # 二項分布のn-sigma区間に入る本数なら最終的な異常位置とする.
    #k_minimum = Binomial(len(z_i), float(sys.argv[])/2)
    sigma_range = Binomial2(len(z_i), 3)
    print("異常候補領域", i, "でのリード総本数:", len(z_i))
    print("異常候補領域", i, "でのn-sigma区間", sigma_range)
    #if k_minimum > len(z_i_hap1) or k_minimum > len(z_i_hap2):
    if len(z_i_hap1) < sigma_range[0] or sigma_range[1] < len(z_i_hap1) or len(z_i_hap2) < sigma_range[0] or sigma_range[1] < len(z_i_hap2):
        print("異常候補", i, "は、二項分布の3σ区間外であり、異常とは認められませんでした。")
        continue
    #elif len(z_i) / :
    #    print("異常候補", i, "は、")

    # 追加
    Z_hap1.append(z_i_hap1)
    Z_hap2.append(z_i_hap2)
    notZ.append(notz_i)



'''____________ファイル作成_____________'''

fastq = sys.argv[7]
fastq_name = os.path.basename(fastq)

#fastq_l = open(fastq).readlines()
ff = open(fastq)
fastq_l = ff.readlines()
for i in range(len(Z_hap1)):
    with open(f'{sys.argv[9]}/{fastq_name}.z{i+1}_hap1', mode='w') as f:
        for j in Z_hap1[i]:
            for k in range(len(fastq_l)):
                if j in fastq_l[k]:
                    s = fastq_l[k] + fastq_l[k+1] + fastq_l[k+2] + fastq_l[k+3]
                    f.write(s)
        for j in notZ[i]:
            for k in range(len(fastq_l)):
                if j in fastq_l[k]:
                    s = fastq_l[k] + fastq_l[k+1] + fastq_l[k+2] + fastq_l[k+3]
                    f.write(s)
ff.close()

#fastq_l = open(fastq).readlines()
ff = open(fastq)
fastq_l = ff.readlines()
for i in range(len(Z_hap2)):
    with open(f'{sys.argv[9]}/{fastq_name}.z{i+1}_hap2', mode='w') as f:
        for j in Z_hap2[i]:
            for k in range(len(fastq_l)):
                if j in fastq_l[k]:
                    s = fastq_l[k] + fastq_l[k+1] + fastq_l[k+2] + fastq_l[k+3]
                    f.write(s)
        for j in notZ[i]:
            for k in range(len(fastq_l)):
                if j in fastq_l[k]:
                    s = fastq_l[k] + fastq_l[k+1] + fastq_l[k+2] + fastq_l[k+3]
                    f.write(s)
ff.close()
