'''///////////////////////////////////////////////////////////////////////

　MAFFTの多重アラインメント結果からコンセンサス配列を作成する行程

////////////////////////////////////////////////////////////////////////'''

print("---CM.py---")

import re
import gc
import sys
import os
import numpy as np
import collections

#sys.argv[1] - fastqファイルパス
#sys.argv[2] - mafft結果
#sys.argv[3] - 何割が'*'でないときからconsに含めるか
#sys.argv[4] - QS閾値割合
#sys.argv[5] - QS閾値リード数
#sys.argv[6] - sub-consensus配列ファイル（新規作成）


def get_aln_reads():
    #アラインメントされたリード(mafft)をリスト構造で作成
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

    return reads, names


def basic_info(reads):
    readSu = len(reads)
    readN = len(reads[0])
    return readSu, readN


def get_QS():
    #fastqファイルからmafft結果.txtに入っているリードのQSを抽出
    f_txt = open(sys.argv[2])
    l_txt = f_txt.readlines()
    plmi_l = []
    name_txt = []
    for i in range(len(l_txt)):
        if l_txt[i][0] == '>':
            plmi = re.findall('\t..?\n', l_txt[i])[0][1:-1:]
            plmi_l.append(int(plmi))
            name = re.split("\t| ", l_txt[i])[0][1:]
            name_txt.append(name)

    f_fq = open(sys.argv[1])
    l_fq = f_fq.readlines()
    qss = {}
    for i in range(len(l_fq)):
        if i % 4 == 0 and l_fq[i][0] == '@':
            name = re.split("\t| ", l_fq[i])[0]
            name = name.replace("\n","")
            name = name[1:]
            if name in name_txt:
                qss[name] = l_fq[i+3].replace("\n","")

    print("\nname_txt\n", name_txt)
    print("\nplmi_l\n", plmi_l)
    print("\nqss\n", qss)
    # reverse QS string in case of -1.
    for i in range(len(name_txt)):
        if plmi_l[i] == -1:
            qss[name_txt[i]] = qss[name_txt[i]][::-1]
    print("\nqss reversed\n", qss)

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
        #~~~~~~~~~~~~~~~QSのリバースを取ることまだしてない！！！~~~~~~~~
    return qs


def start_end(readSu, readN, reads):
    #各リードはどこからどこまでか
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
    return starts, ends


def convert_to_asterisk(starts, ends, readSu, readN, reads):
    #各リードでstart-end外の「-」-->「*」に変換
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
    return reads


def get_start_and_end_point_of_cs(readSu, readN, reads):
    #切り落とす場所
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

    return cons_start, cons_end


def move_QS(readSu, cons_start, reads):
    #cons_start, cons_endまでQSを進めておく
    qs_pos = [0 for i in range(readSu)]     #0-origin index
    for i in range(cons_start):
        for j in range(readSu):
            if reads[j][i] != '*' and reads[j][i] != '-':
                qs_pos[j] += 1
    return qs_pos


def make_cs(cons_start, cons_end, readSu, reads, qs, qs_pos, names):
    # CS配列作成
    consensus = ''
    base_info = '#number of bases'+"\t\t\t\t"+"proportion of bases"+"\t\t\t"+"sum of QS"+"\n"
    base_info += "#\tA\tG\tC\tT\t-\tA\tG\tC\tT\t-\tA\tG\tC\tT\n"
    base_cnt = 0
    irr_num_l = [0 for i in range(readSu)]#
    for i in range(cons_start, cons_end+1):     #前後XX塩基分は切り捨てる
        base_cnt += 1
        cnt_type = {'a':0, 't':0, 'g':0, 'c':0, '-':0, '*':0}
        base_info_i = ''
        for j in range(readSu):
            cnt_type[reads[j][i]] += 1

        if cnt_type['a'] == 0 and cnt_type['t'] == 0 and cnt_type['g'] == 0 and cnt_type['c'] == 0:
            continue

        qs_sum = {'a':0, 't':0, 'g':0, 'c':0}
        for j in range(readSu):
            if reads[j][i] in ('a', 't', 'g', 'c'):
                qs_sum[reads[j][i]] += ord(qs[j][qs_pos[j]]) - 33   # ASCII - 33
                qs_pos[j] += 1

        numEffectiveReads = cnt_type["a"]+cnt_type["g"]+cnt_type["c"]+cnt_type["t"]+cnt_type["-"]
        base_info_i += str(base_cnt)+"\t"+str(cnt_type["a"])+"\t"+str(cnt_type["g"])+"\t"+str(cnt_type["c"])+"\t"+str(cnt_type["t"])+"\t"+str(cnt_type["-"])+"\t"
        base_info_i += str(round(cnt_type["a"]/numEffectiveReads, 3)) + "\t" + str(round(cnt_type["g"]/numEffectiveReads, 3)) + "\t" + str(round(cnt_type["c"]/numEffectiveReads, 3)) + "\t" + str(round(cnt_type["t"]/numEffectiveReads, 3)) + "\t" + str(round(cnt_type["-"]/numEffectiveReads, 3)) + "\t"
        base_info_i += str(qs_sum["a"])+"\t"+str(qs_sum["g"])+"\t"+str(qs_sum["c"])+"\t"+str(qs_sum["t"])
        base_info += base_info_i + "\n"

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
    #        consensus += 'N'
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
    print('Read name:\n', names)
    print('Adoption proportion:\n', irr_num_l)#
    

    #sub-consensus配列の最後は改行しておく 
    consensus += '\n'

    f = open(sys.argv[6], mode='w')
    out8_name = sys.argv[2]
    out4_name = os.path.basename(out8_name)
    out4_name = out4_name.replace("out8", "cs")

    s = f'>{out4_name}\t{readSu}\n'     #read number (tab delimited)
    s += consensus
    f.write(s)
    f.close()

    with open(sys.argv[6].replace("cs","out4a"), mode="w") as f:
        f.write(base_info)
    print("made "+sys.argv[6]+" annotation file.")


def main():
    reads, names = get_aln_reads()
    readSu, readN = basic_info(reads)
    print("readSu:", readSu)
    print("readN:", readN)
    qs = get_QS()
    starts, ends = start_end(readSu, readN, reads)
    reads = convert_to_asterisk(starts, ends, readSu, readN, reads)
    cons_start, cons_end = get_start_and_end_point_of_cs(readSu, readN, reads)
    qs_pos = move_QS(readSu, cons_start, reads)
    make_cs(cons_start, cons_end, readSu, reads, qs, qs_pos, names)
    return print("main finished.")

if __name__ == "__main__":
    main()

