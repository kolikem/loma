# Written by Python on 230117.
# Concatenate sub-consensus sequences to a CS.

import sys
import os
import subprocess
import re

# sys.argv[1] - common name: e.g) /Users/ikemotoko/Desktop/loma_v2/out/dir1/ST3GAL3_novel_0_ATL7
# sys.argv[2] - number of out4 (#blocks)
# sys.argv[3] - block size (window size)
# sys.argv[4] - step size
# sys.argv[5] - mafft path

print("input parameters")
print("block size:", sys.argv[3])
print("step size:", sys.argv[4])
print("number of blocks:", sys.argv[2])


def make_init_concat():
    cmd = "cp "+sys.argv[1]+".out4.1"+" "+sys.argv[1]+".concat.1"
    subprocess.check_call(cmd, shell=True)
    print("made an initial concat. file:", cmd)
    return 0

def make_tail_head_file(CNT, name, dir_name):
    # the tail of concat.{CNT} & the head of out4.{CNT+1} > out
    with open(dir_name+"/"+name+".concat."+str(CNT-1)) as f:
        for line in f:
            line = line.replace("\n","")
            if line[0] == ">":
                left_name = line
            else:
                left_seq = line

    with open(dir_name+"/"+name+".out4."+str(CNT)) as f:
        for line in f:
            line = line.replace("\n","")
            if line[0] == ">":
                right_name = line
            else:
                right_seq = line

    # make a new file
    OUT = dir_name+"/"+name+".tail_head."+str(CNT)
    with open(OUT, mode="w") as f:
        content = left_name + "\n" + left_seq[-(int(sys.argv[3])-int(sys.argv[4])):] + "\n"     # 後 blockSize-stepSize[nt]
        content += right_name + "\n" + right_seq[:int(sys.argv[3])-int(sys.argv[4])]            # 前 blockSize-stepSize[nt]
        f.write(content)
    print("made a new file:", OUT)
    return OUT, len(left_seq), len(right_seq), left_name.split("\t")[1], right_name.split("\t")[1]  #出力ファイル名, 配列長, ﾌﾞﾛｯｸｶﾊﾞﾚｯｼﾞ情報を返す.


def align_tail_head_file(FILE, mafft):
    cmd = mafft+" "+FILE+" > "+FILE+".mafft"
    subprocess.check_call(cmd, shell=True)
    print("ran mafft:", cmd)
    return FILE+".mafft"    #出力ファイル名を返す.



def find_middle_location(th_mafft, len_left_seq, len_right_seq):
    with open(th_mafft) as f:
        cnt = 0
        l_seq, r_seq, nam = '', '', []
        for line in f:
            line = line.replace("\n","")
            if line[0] == '>':
                cnt += 1
                nam.append(line)

            if cnt == 1 and line[0] != '>':
                l_seq += line
            elif cnt == 2 and line[0] != '>':
                r_seq += line

    loc_left = [len_left_seq-(int(sys.argv[3])-int(sys.argv[4])), len_left_seq+1]   # left_seqの何塩基目から何塩基目までがMIDDLEか  (1-origin)
    loc_right = [0, int(sys.argv[3])-int(sys.argv[4])+1]                            # right_seqの何塩基目から何塩基目までがMIDDLEか (1-origin)
    for i in range(len(l_seq)):
        if l_seq[i] != "-":
            loc_left[0] += 1
            if r_seq[i] != "-":
                break
    for i in range(len(l_seq)):
        if l_seq[-1-i] != "-":
            loc_left[1] -= 1
            if r_seq[-1-i] != "-":
                break
    for i in range(len(l_seq)):
        if r_seq[i] != "-":
            loc_right[0] += 1
            if l_seq[i] != "-":
                break
    for i in range(len(l_seq)):
        if r_seq[-1-i] != "-":
            loc_right[1] -= 1
            if l_seq[-1-i] != "-":
                break

    return loc_left, loc_right


def Concatenate(loc_left, loc_right, CNT, name, dir_name):
    with open(dir_name+"/"+name+".concat."+str(CNT-1)) as f:
        for line in f:
            line = line.replace("\n","")
            if line[0] == ">":
                left_name = line
            else:
                left_seq = line

    with open(dir_name+"/"+name+".out4."+str(CNT)) as f:
        for line in f:
            line = line.replace("\n","")
            if line[0] == ">":
                right_name = line
            else:
                right_seq = line

    LEFT = left_seq[:loc_left[0]-1]                 # 次世代concatのLEFT配列
    MIDDLE = left_seq[loc_left[0]-1:loc_left[1]]     # 次世代concatのMIDDLE配列 (left配列を配列を採用する)
    RIGHT = right_seq[loc_right[1]:]

    return LEFT+MIDDLE+RIGHT



num_out4 = int(sys.argv[2])
dir_name = os.path.dirname(sys.argv[1])
base_name = os.path.basename(sys.argv[1])
make_init_concat()
CNT = 1     # out4を右でどこまで繋いだか. CNT = sys.argv[2] で終了(.cs).
for i in range(num_out4-1):
    CNT += 1
    OUT, len_left_seq, len_right_seq, cov_l, cov_r = make_tail_head_file(CNT, base_name, dir_name)
    OUT_mafft = align_tail_head_file(OUT, sys.argv[5])
    loc_left, loc_right = find_middle_location(OUT_mafft, len_left_seq, len_right_seq)
    concat_CNT_seq = Concatenate(loc_left, loc_right, CNT, base_name, dir_name)
    if i < num_out4-2:
        name_line = ">"+base_name+".concat."+str(CNT)+"\t"+cov_l+","+cov_r
        with open(dir_name+"/"+base_name+".concat."+str(CNT), mode="w") as f:
            content = name_line + "\n" + concat_CNT_seq
            f.write(content)
    else:
        name_line = ">"+base_name+".cs"+"\t"+cov_l+","+cov_r
        with open(dir_name+"/"+base_name+".cs", mode="w") as f:
            content = name_line + "\n" + concat_CNT_seq
            f.write(content)

print("concatenation finished.")
print("a consensus sequence (.cs) finished.")

