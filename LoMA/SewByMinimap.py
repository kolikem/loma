'''//////////////////////////////////////////
minimapアラインメント(paf)から
一本のコンセンサス配列につなぎ合わせる
//////////////////////////////////////////'''

print("---SewByMinimap.py---")

import sys
import os
import glob
import re

# sys.argv[1] - {name}_two_N.txt パス
# sys.argv[2] - {name}_two_minimap_N.paf パス
# sys.argv[3] - {name}_N.txt パス


with open(sys.argv[1]) as two:
    name = []
    seq = []
    for line in two:
        if line[0] == '>':
            m = re.search(r'>.*_[0-9]*\t', line)
            if m == None:
                m = re.search(r'>.*\t', line)
            nam = m.group()
            nam = nam[1:-1]
            name.append(nam)
        else:
            seq.append(line[:-1])
#print('name', name)
#print('seq', seq)

sewnSeq = ''
f = open(sys.argv[2])
for line in f:
    l = line.split('\t')
    if l[0] == name[0] and l[5] == name[1] and l[4] == '+':
        if int(l[3])/int(l[1]) > 0.95: #and 2 * (int(l[1])-int(l[3])) < int(l[6])-int(l[8]):
            sewnSeq += seq[0][:int(l[3])+1]
            sewnSeq += seq[1][int(l[8])+1:int(l[6])]
            break
    elif l[0] == name[1] and l[5] == name[0] and l[4] == '+':
        if int(l[8])/int(l[6]) > 0.95: #and 2 * (int(l[6])-int(l[8])) < int(l[1])-int(l[3]):
            print('範囲チェック, l[8], l[3], l[1]', l[8], l[3], l[1])
            sewnSeq += seq[0][:int(l[8])+1]
            sewnSeq += seq[1][int(l[3])+1:int(l[1])]
            break
f.close()
print('sewnSeq length', len(sewnSeq))


with open(sys.argv[1]) as two:
    for line in two:
        if line[0] == '>':
            m = re.search(r'_[0-9]*\t', line)
            if m != None:
                num = m.group()
                num = int(num[1:-1])

            m1 = re.search(r'[a-zA-Z]*[0-9]*_chr[0-9XY]*_*[0-9]*-[0-9]*', line)
            if m1 != None:
                nam = m1.group()
            
            m2 = re.search(r'\t[0-9]*\t.*\n', line)
            if m2 != None:
                tag = m2.group()
                tag = tag[1:-1]
                continue

            m3 = re.search(r'\t[0-9]*\n', line)
            if m3 != None:
                add = m3.group()
                add = add[1:-1]
num += 1
num = str(num)

with open(sys.argv[3], 'w') as out:
    #名前記入
    out.write(f'>{nam}_{num}\t{tag}\t{add}\n')
    #配列記入
    out.write(sewnSeq + '\n')

