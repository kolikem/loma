'''///////////////////////////////////////////////////////////////////////

　MAFFTの多重アラインメント結果から部分コンセンサス配列を作成する行程

////////////////////////////////////////////////////////////////////////'''

print("---ConsMaker.py---")

import re
import gc
import sys
import os
import numpy as np
import collections

#sys.argv[1] - fastqファイルパス
#sys.argv[2] - mafft結果.txt名前
#sys.argv[3] - 何割が'*'でないときからconsに含めるか
#sys.argv[4] - QS閾値割合
#sys.argv[5] - QS閾値リード数
#sys.argv[6] - sub-consensus配列ファイル（新規作成）



'''////////// アラインメントされたリード集をリスト構造で作成 //////////'''
f = open(sys.argv[2])
l = f.readlines() 

name_num = []
names = []
for i in range(len(l)):
    if l[i][0] == '>':
        name_num.append(i)
        names.append(l[i][1:-1])

reads = []
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


'''////////// 変数定義 //////////'''
readSu = len(reads)     #リード数
readN = len(reads[0])   #アラインメントリード長

print('mafftファイル内リード数:', readSu)

'''////////// fastqファイルからmafft結果.txtに入っているリードのQSを抜き出す //////////'''
f_txt = open(sys.argv[2])
l_txt = f_txt.readlines()
name_txt = []
for i in range(len(l_txt)):
    if l_txt[i][0] == '>':
        name = re.findall('>.*?\t', l_txt[i])
        name = name[0][1:-1:]
        name_txt.append(name)

f_fq = open(sys.argv[1])
l_fq = f_fq.readlines()
qss = {}
for i in range(len(l_fq)):
    if i % 4 == 0 and l_fq[i][0] == '@':
        name = re.findall('@.*? ', l_fq[i])
        name = name[0][1:-1:]
        if name in name_txt:
            qss[name] = l_fq[i+3][:-1]
        else:
            continue

qs = []
for i in range(len(name_txt)):
    qs_i = qss[name_txt[i]]
    qs.append(qs_i)

f_txt.close()
f_fq.close()
del l_txt
gc.collect()
del l_fq
gc.collect()
del name_txt
gc.collect()
del qss
gc.collect()

#print('Quality scores: \n', qs)


'''////////// QSのキリトリ //////////'''
f_txt = open(sys.argv[2])
l_txt = f_txt.readlines()
plmi = []
kaishi = []
for i in range(len(l_txt)):
    if l_txt[i][0] == '>':
        plmi_i = re.findall('\t..?\n', l_txt[i])
        plmi_i = plmi_i[0][1:-1:]
        plmi.append(plmi_i)
        kaishi_i = re.findall('\t.*\t', l_txt[i])
        kaishi_i = kaishi_i[0][2:-2:]
        kaishi.append(kaishi_i) 
f_txt.close()
del l_txt
gc.collect()

for i in range(readSu):
    kaishi[i] = kaishi[i].split(', ')
    kaishi[i][0] = int(kaishi[i][0])
    kaishi[i][1] = int(kaishi[i][1])

for i in range(readSu):
    if int(plmi[i]) == -1:
        qs_i = qs[i][::-1]
        #qs_i = qs_i[1:]
    else:
        qs_i = qs[i] #+[0:-2]

    qs_truncated_i = qs_i[kaishi[i][0]:kaishi[i][1]+1]
    qs[i] = qs_truncated_i


'''////////// 各リードはどこからどこまでか //////////'''
starts = []
ends = []
for i in range(readSu):
    start_i = 0
    for j in range(readN):
        if reads[i][j] == '-':
            start_i += 1
        else:
            break
    starts.append(start_i)  #startsリストにstart位置(0-origin)を格納

    end_i = len(reads[i])-1
    for j in range(readN):
        if reads[i][-j-1] == '-':
            end_i -= 1
        else:
            break
    ends.append(end_i)      #endsリストにend位置(0-origin)を格納


'''////////// 各リードでstart-end外の「-」-->「*」に変換 //////////'''
for i in range(readSu):
    read_i = list(reads[i])
    for j in range(0, starts[i]):
        read_i[j] = '*'
    for j in range(0, readN-ends[i]-1):
        read_i[-j-1] = '*'

    read_i = ''.join(read_i)
    reads[i] = read_i

# * は始まっていないギャップ
# - は始まっているギャップ


'''////////// 切り落とす場所 //////////'''
i = 0
cons_start = 0
while i == cons_start:
    d = {'*':0, '-':0, 'a':0, 'g':0, 'c':0, 't':0}
    for k in range(readSu):
        d[reads[k][i]] += 1
    i += 1
    if 1 - d['*']/readSu < float(sys.argv[3]):
        cons_start += 1

j = readN - 1
cons_end = readN - 1
while j == cons_end:
    d = {'*':0, '-':0, 'a':0, 'g':0, 'c':0, 't':0}
    for k in range(readSu):
        d[reads[k][j]] += 1
    j -= 1
    if 1 - d['*']/readSu < float(sys.argv[3]):
        cons_end -= 1


'''/////////// cons_start, cons_endまでQSを進めておく //////////'''
qs_pos = [0 for i in range(readSu)]     #0-origin index
for i in range(cons_start):
    for j in range(readSu):
        if reads[j][i] != '*' and reads[j][i] != '-':
            qs_pos[j] += 1



#9/17 変なリードを取り除く  ・・・・# の行　追加


'''////////// 部分コンセンサス配列作成 //////////'''
consensus = ''
print('ブロック名', sys.argv[2])#
irr_num_l = [0 for i in range(readSu)]#
for i in range(cons_start, cons_end+1):     #前後XX塩基分は切り捨てる
    #print('i', i)

    cnt_type = {'a':0, 't':0, 'g':0, 'c':0, '-':0, '*':0}
    for j in range(readSu):
        cnt_type[reads[j][i]] += 1

    if cnt_type['a'] == 0 and cnt_type['t'] == 0 and cnt_type['g'] == 0 and cnt_type['c'] == 0:
        continue

    qs_sum = {'a':0, 't':0, 'g':0, 'c':0}
    for j in range(readSu):
        if reads[j][i] in ('a', 't', 'g', 'c'):
            qs_sum[reads[j][i]] += ord(qs[j][qs_pos[j]]) - 33   # ASCII - 33
            #print('i', i, '塩基', reads[j][i], 'QS', ord(qs[j][qs_pos[j]])-33)
            qs_pos[j] += 1
    #print('i', i, 'qs_sum', qs_sum)
    #print('qs_pos', qs_pos)

    threshold_QS = float(sys.argv[4])
    threshold_su = float(sys.argv[5])
    if max(qs_sum['a'], qs_sum['t'], qs_sum['g'], qs_sum['c']) / (qs_sum['a'] + qs_sum['t'] + qs_sum['g'] + qs_sum['c']) >= threshold_QS and cnt_type['a'] + cnt_type['t'] + cnt_type['g'] + cnt_type['c'] >= threshold_su * readSu:
        if qs_sum['a'] >= qs_sum['t'] and qs_sum['a'] >= qs_sum['g'] and qs_sum['a'] >= qs_sum['c']:
            #consensus += 'a'
            con = 'a'
        elif qs_sum['t'] >= qs_sum['a'] and qs_sum['t'] >= qs_sum['g'] and qs_sum['t'] >= qs_sum['c']:
            #consensus += 't'
            con = 't'
        elif qs_sum['g'] >= qs_sum['a'] and qs_sum['g'] >= qs_sum['t'] and qs_sum['g'] >= qs_sum['c']:
            #consensus += 'g'
            con = 'g'
        else:
            #consensus += 'c'
            con = 'c'
    #else:
    #    consensus += 'N'
    else:
        con = '-'

    if con != '-':
        consensus += con

    for j in range(readSu):#
        if con != '-' and reads[j][i] != con:#
            irr_num_l[j] += 1#



#print('ブロック範囲　コンセンサス長', cons_end-cons_start, len(consensus))
for i in range(readSu):#
    irr_num_l[i] = 100 - int(100 * irr_num_l[i] / len(consensus))#
print('リード名\n', names)
print('リードの採用率\n', irr_num_l)#
    
  


#sub-consensus配列の最後は改行しておく 
consensus += '\n'

f = open(sys.argv[6], mode='w')
#out8(mafftファイル)の名前から_out8の直前までを抜き取る
out8_name = sys.argv[2]
out4_name = os.path.basename(out8_name)
out4_name = out4_name[:-5]

s = f'>{out4_name}\t{readSu}\n'     #10/2 リード数付け足し(tab delimited)
s += consensus
f.write(s)
f.close()


'''///// Write out sub-consensus.txt files for minimap2 /////'''
'''
for i in range(len(subCons)):
    if i < 10:
        name = sys.argv[6] + f'/subCons_{i}.fa' #0{i}.fa'
    else:
        name = sys.argv[6] + f'/subCons_{i}.fa'
    f = open(name, mode='w')
    s = f'>{i}\n'
    s += subCons[i]
    f.write(s)
    f.close()
'''

