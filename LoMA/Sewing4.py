'''//////////////////////////////////////////
        ２本の部分コンセンサスを統合
//////////////////////////////////////////'''

import sys
import re
import glob
import os

# sys.argv[1] - mafft結果パス

seq_l = []
nam = []
with open(sys.argv[1]) as f:
    cnt = 0
    seq1 = ''
    seq2 = ''
    for line in f:
        if line[0] == '>':
            cnt += 1
            nam.append(line)

        if cnt == 1 and line[0] != '>':
            seq1 += line[:-1]

        elif cnt == 2 and line[0] != '>':
            seq2 += line[:-1]



seq_left, seq_right = seq1, seq2
cnt = [0, 0, 0]
#cnt[1] - 重なりの始め
for i in range(len(seq_left)):
    if seq_left[i] != '-' and seq_right[i] != '-':
        fst = i
        break
cnt[1] = fst
for i in range(len(seq_left)):
    if seq_left[i] != '-' and seq_right[i] != '-':
        if i + 100 <= len(seq_left):
            dens = 0
            for j in range(i, i+100):
                if seq_left[j] != '-' and seq_right[j] != '-':
                    dens += 1
            if dens > 70:
                cnt[1] = i
                break

#cnr[2] - 重なりの終わり
for i in range(len(seq_left)):
    if seq_left[-1-i] != '-' and seq_right[-1-i] != '-':
        fst = len(seq_left) - i
        break
cnt[2] = fst
for i in range(len(seq_left)):
    if seq_left[-1-i] != '-' and seq_right[-1-i] != '-':
        if len(seq_left) - i - 100 >= 0:
            dens = 0
            for j in range(len(seq_left)-i-100, len(seq_left)-i):
                if seq_left[j] != '-' and seq_right[j]  != '-':
                    dens += 1
            if dens > 70:
                cnt[2] = len(seq_left) - i
                break
#print(cnt)

#重なり方が期待通りかどうか確かめる
dens = 0
for i in range(cnt[1], cnt[2]):
    if seq_left[i] != '-' and seq_right[i] != '-':
        dens += 1

tail = 0
for i in range(cnt[2], len(seq_left)):
    if seq_left[i] == '-' and seq_right[i] != '-':
        tail += 1


if dens / (cnt[2] - cnt[1] + 1) < 0.7 or tail / (len(seq_left) - cnt[2] + 1) < 0.7:
    print(1)
else:
    print(0)
