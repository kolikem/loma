# Written by Python on 230117.
# Concatenate sub-consensus sequences to a CS.

import sys
import os
import subprocess
import re

# sys.argv[1] - input swap file path:       /home/kikemoto/mafft_para_0625/output1-15/1/dir1/NA18943_chr10_93121869-93161869_55.swap
# sys.argv[2] - output tail-head file path: /home/kikemoto/mafft_para_0625/output1-15/1/dir1
# sys.argv[3] - block size (window size): 3000
# sys.argv[4] - step size: 1000
# sys.argv[5] - mafft path
# sys.argv[6] - consensus file path
# sys.argv[7] - consensus番号(新)

print("input parameters")
print("block size: ", sys.argv[3])
print("step size: ", sys.argv[4])
print("input file name: ", sys.argv[1])
print("output file location: ", sys.argv[2])
print("consensus file location: ", sys.argv[6])
print("consensus number: ", sys.argv[7])


with open(sys.argv[1], mode="r") as f:
	cnt = 0
	for line in f:
		cnt += 1
		if cnt == 1:
			left_name = line
		elif cnt == 2:
			if line[-1] == "\n":
				left_seq = line[:-1]
			else:
				left_seq = line
		elif cnt == 3:
			right_name = line
		elif cnt == 4:
			if line[-1] == "\n":
				right_seq = line[:-1]
			else:
				right_seq = line
		else:
			break

	# tail-headファイルを作成.
	tail_head = open(sys.argv[2]+"/"+os.path.basename(sys.argv[1])+".tail_head", mode="w")
	content = left_name + left_seq[int(sys.argv[4])-int(sys.argv[3]):]+"\n" + right_name + right_seq[:int(sys.argv[3])-int(sys.argv[4])]	# 前方配列の後2000bp, 後方配列の前2000bp.
	tail_head.write(content)
	tail_head.close()

cmd = sys.argv[5] + " " + sys.argv[2]+"/"+os.path.basename(sys.argv[1])+".tail_head" + ">" + sys.argv[2]+"/"+os.path.basename(sys.argv[1])+".tail_head" + ".mafft"
subprocess.check_call(cmd, shell=True)

# middle
with open(sys.argv[2]+"/"+os.path.basename(sys.argv[1])+".tail_head.mafft", mode="r") as f:
	cnt = 0
	for line in f:
		cnt += 1
		if cnt == 2:
			l = line[:-1]
		elif cnt == 4:
			r = line

nam = []
with open(sys.argv[2]+"/"+os.path.basename(sys.argv[1])+".tail_head.mafft", mode="r") as f:
    cnt = 0
    l = ''
    r = ''
    for line in f:
        if line[0] == '>':
            cnt += 1
            nam.append(line)

        if cnt == 1 and line[0] != '>':
            l += line[:-1]

        elif cnt == 2 and line[0] != '>':
            r += line[:-1]

cnt = [0, 0, 0]
	#cnt[1] - 重なりの始め
for i in range(len(l)):
	if l[i] != '-':
		cnt[1] = i
		break

	#cnr[2] - 重なりの終わり
for i in range(len(l)):
	if r[-1-i] != '-':
		cnt[2] = len(l)-i-1
		break


# 接続
left = left_seq[:int(sys.argv[4])-int(sys.argv[3])]
right = right_seq[int(sys.argv[3])-int(sys.argv[4]):]
middle_l = l[cnt[1]:cnt[2]+1]
middle_r = r[cnt[1]:cnt[2]+1]
middle = ''
for i in range(len(middle_l)):
	if middle_l[i] == '-':
		middle += middle_r[i]
	elif middle_r[i] == '-':
		middle += middle_l[i]
	else:
		middle += middle_l[i]

concat = left + middle + right

# concat配列格納
a = re.search(r'\t[0-9]*.*\n' , nam[0])
a = a.group()
a = a[1:-1]
b = re.search(r'\t[0-9]*\.*\n' , nam[1])
b = b.group()
b = b[1:-1]
ab = a + '\t' + b

name = os.path.basename(sys.argv[1])
name = name[:-5]
with open(sys.argv[6], mode='w') as f:
	s = f'>{name}\t{ab}\n'
	s += concat
	s += '\n'
	f.write(s)

cons = os.path.basename(sys.argv[6])
m = re.search(r'consensus_.*_', cons)
m = m.group()
cmd = "mv " + sys.argv[6] + " " + os.path.dirname(sys.argv[6]) + "/" + m + sys.argv[7]
subprocess.check_call(cmd, shell=True)


