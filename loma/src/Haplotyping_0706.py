'''///////////////////////////////////////////////////////////
        Haplotyping 分離　可能・不可能　判断
////////////////////////////////////////////////////////////'''

print("---Haplotyping.py---")

import sys
import glob
import re
import os

#sys.argv[1] - out4 & out8 (MAFFT２周目)ディレクトリ
#sys.argv[2] - name : NA18943_chr1_79552999-79572999

direc = sys.argv[1]
if direc[-1] == '/':
    direc = direc[:-1]
l = glob.glob(direc + '/' + sys.argv[2] + '*out8*')
#print('l', l)


#ヘテロ分類ブロックが存在することを確認
cnt = 0
num_h1 = []
num_h2 = []
for i in range(len(l)):
    if '_h1_' in l[i]:
        cnt += 1
        num_h1.append(i)
    elif '_h2_' in l[i]:
        num_h2.append(i)
#print('cnt', cnt)
#print(num_h1, num_h2)

if cnt == 0:
    print(sys.argv[2], '# Hetero classified block: 0')
    #print('分離可能')

elif cnt == 1:
    print(sys.argv[2], '# Hetero classified block: 1')
    #print('分離可能')
    out8 = l[num_h1[0]]
    out4 = out8.replace('_h1_out8', '_h1_out4')
    out4_replaced = out4.replace('_h1_out4', '_h1_typed_out4')
    #print('out4', out4)
    #print('out4_replaced', out4_replaced)
    os.rename(out4, out4_replaced)
    out8 = l[num_h2[0]]
    out4 = out8.replace('_h2_out8', '_h2_out4')
    out4_replaced = out4.replace('_h2_out4', '_h2_typed_out4')
    #print('out4', out4)
    #print('out4_replaced', out4_replaced)
    os.rename(out4, out4_replaced)



else:
    print(sys.argv[2], '# Hetero classified block:', cnt)
    #名前を順序づけ
    h1_name_l = [[l[i], i] for i in num_h1]
    h2_name_l = [[l[i], i] for i in num_h2]
    
    for i in range(len(h1_name_l)):
        m = h1_name_l[i][0]
        m = re.search(r'_-?[0-9]*--?[0-9]*_h1_out8', m)
        m = m.group()
        m = re.search(r'_-?[0-9]*-', m)
        m = m.group()
        m = m[1:-1]
        m = int(m)
        h1_name_l[i].append(m)
    #print('h1_name_l', h1_name_l)

    for i in range(len(h2_name_l)):
        m = h2_name_l[i][0]
        m = re.search(r'_-?[0-9]*--?[0-9]*_h2_out8', m)
        m = m.group()
        m = re.search(r'_-?[0-9]*-', m)
        m = m.group()
        m = m[1:-1]
        m = int(m)
        h2_name_l[i].append(m)
    #print('h2_name_l', h2_name_l)

    h1_name_sorted = sorted(h1_name_l, key=lambda x: x[2])
    #print('h1_name_sorted', h1_name_sorted)

    h2_name_sorted = sorted(h2_name_l, key=lambda x: x[2])
    #print('h2_name_sorted', h2_name_sorted)

    #順序づけされた順番でリード名を列挙
    h1_reads = []   #sortedの順にリード名を格納
    for i in range(len(h1_name_sorted)):
        h1s = []
        with open(h1_name_sorted[i][0]) as h1_out8:
            for line in h1_out8:
                if line[0] == '>':
                    line = re.search(r'>.*\t\(', line)
                    line = line.group()
                    line = line[1:-2]
                    h1s.append(line)
        h1_reads.append(h1s)
        
    h2_reads = []
    for i in range(len(h2_name_sorted)):
        h2s = []
        with open(h2_name_sorted[i][0]) as h2_out8:
            for line in h2_out8:
                if line[0] == '>':
                    line = re.search(r'>.*\t\(', line)
                    line = line.group()
                    line = line[1:-2]
                    h2s.append(line)
        h2_reads.append(h2s)

    #print('h1_reads', h1_reads)
    #print('h2_reads', h2_reads)

    #つながりを探す
    haplotyping1 = [1]
    for i in range(1, len(h1_name_sorted)):
        h1_h1 = set(h1_reads[i-1]).intersection(set(h1_reads[i]))   #タイプ1
        h1_h2 = set(h1_reads[i-1]).intersection(set(h2_reads[i]))   #タイプ2
        h2_h1 = set(h2_reads[i-1]).intersection(set(h1_reads[i]))   #タイプ2
        h2_h2 = set(h2_reads[i-1]).intersection(set(h2_reads[i]))   #タイプ1

        if haplotyping1[i-1] == 1:
            if len(h1_h1) != 0 and len(h1_h1) >= len(h1_h2):
                haplotyping1.append(1)
            elif len(h1_h2) != 0 and len(h1_h2) > len(h1_h1):
                haplotyping1.append(2)
            elif len(h1_h1) == 0 and len(h1_h2) == 0:
                print('1 分離不可能')
                break
        elif haplotyping1[i-1] == 2:
            if len(h2_h1) != 0 and len(h2_h1) >= len(h2_h2):
                haplotyping1.append(1)
            elif len(h2_h2) != 0 and len(h2_h2) > len(h2_h1):
                haplotyping1.append(2)
            elif len(h2_h1) == 0 and len(h2_h2) == 0:
                print('2 分離不可能')
                break

    haplotyping2 = []
    for k in haplotyping1:
        if k == 1:
            haplotyping2.append(2)
        elif k == 2:
            haplotyping2.append(1)

    #print('haplotyping1', haplotyping1)
    #print('分離可能')
    for i in range(len(haplotyping1)):
        #print('i', i)
        if i == 0:
            out8 = h1_name_sorted[0][0]
            m = out8[:-5]
            out4 = m + '4.txt'
            out4_replaced = out4.replace('_h1_', '_h1_typed_')
            #print('out8', out8)
            #print('out4', out4)
            #print('out4_replaced', out4_replaced)
            os.rename(out4, out4_replaced)

            out8 = h2_name_sorted[0][0]
            m = out8[:-5]
            out4 = m + '4.txt'
            out4_replaced = out4.replace('_h2_', '_h2_typed_')
            #print('out8', out8)
            #print('out4', out4)
            #print('out4_replaced', out4_replaced)
            os.rename(out4, out4_replaced)

        else:
            if haplotyping1[i] == 1:
                out8 = h1_name_sorted[i][0]
                m = out8[:-5]
                out4 = m + '4.txt'
                out4_replaced = out4.replace('_h1_', '_h1_typed_')
                #print('out8', out8)
                #print('out4', out4)
                #print('out4_replaced', out4_replaced)
                os.rename(out4, out4_replaced)

                out8 = h2_name_sorted[i][0]
                m = out8[:-5]
                out4 = m + '4.txt'
                out4_replaced = out4.replace('_h2_', '_h2_typed_')
                #print('out8', out8)
                #print('out4', out4)
                #print('out4_replaced', out4_replaced)
                os.rename(out4, out4_replaced)

            elif haplotyping1[i] == 2:
                out8 = h2_name_sorted[i][0]
                m = out8[:-5]
                out4 = m + '4.txt'
                out4_replaced = out4.replace('_h2_', '_h1_typed_')
                #print('out8', out8)
                #print('out4', out4)
                #print('out4_replaced', out4_replaced)
                os.rename(out4, out4_replaced)

                out8 = h1_name_sorted[i][0]
                m = out8[:-5]
                out4 = m + '4.txt'
                out4_replaced = out4.replace('_h1_', '_h2_typed_')
                #print('out8', out8)
                #print('out4', out4)
                #print('out4_replaced', out4_replaced)
                os.rename(out4, out4_replaced)

