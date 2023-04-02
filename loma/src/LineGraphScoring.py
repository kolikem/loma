'''//////////////////////////////////////////////////////////

            MAFFT後　各ブロック内リードの分類
    アセンブリエラー・マッピングエラー　ヘテロ接合

//////////////////////////////////////////////////////////'''

print("---LineGraphScoring.py---")

import os
import re
import gc
import sys
import numpy as np
import matplotlib.pyplot as plt
#from scipy.special import comb		<-- Segmentation fault (コアダンプ)原因
print("imported modules")

#name = name[1:-9]       #e.g. 'NA18943_chrX_126288870-126308870_6663-8162'
# sys.argv[1] - mafft output (out3)
# sys.argv[2] - 結果図格納場所
# sys.argv[3] - 分類後fastaファイル保存ディレクトリ
# sys.argv[4] - ブロックサイズ

'''///// 分配確率 /////'''      #10/14 保留
'''
def Prob(n, l, m):
    if l >= m:
        summ = 0
        for i in range(m+1):
            summ += comb(n, i)
        return summ / 2**n
    else:
        summ = 0
        for i in range(l+1):
            summ += comb(n, i)
        return summ / 2**n
'''


'''///// データ下準備 /////'''
f = open(sys.argv[1])
l = f.readlines()
#print("l", l)
name_num = []
names = []
for i in range(len(l)):
    if l[i][0] == '>':
        name_num.append(i)
        names.append(l[i][1:-1])
print("names, name_num", names, name_num)
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




# リード数×リード数の行列を作成（y_l格納用）
num = len(names)                            #行列の行数＝列数＝リード数
print('リード数', num)

# divide: ブロックのリード分類を行うかどうかの境界本数
divide = 10          # ここはパラメータ

# divide 本未満なら終了
if num < divide:
    name = sys.argv[1]
    name = os.path.basename(name)
    name = name[:-9]	# e.g. NA19240_chr15_75493250-75514000_45000-47999

    with open(f'{sys.argv[3]}/{name}_out7.fa', mode='w') as f:
        s = ''
        for i in range(num):
            s += f'>{names[i]}\n'
            read_i = reads[i].replace('-', '')
            s += f'{read_i}\n'
        f.write(s)

# divide 本以上なら分類
elif num >= divide:
    '''///// スコア行列 /////'''
    SM = [[] for i in range(num)]               # SM: Score Matrix
    for i in range(num):
        SM[i] = [[] for j in range(num)]        # SM: 各要素はリスト[]
    
    for i in range(num):
        name1, seq1 = names[i], reads[i]
        for j in range(i, num):
            name2, seq2 = names[j], reads[j]

            #２本のリードの前後の'-'連続は省く
            leftmost = 0
            rightmost = len(reads[i]) - 1
            l_notgap = 0
            while l_notgap < 5:
                l_notgap = 0
                for k in range(leftmost, leftmost+8):
                    if reads[i][k] != '-' and reads[j][k] != '-':
                        l_notgap += 1
                leftmost += 1
            r_notgap = 0
            while r_notgap < 5:
                r_notgap = 0
                for k in range(rightmost-8, rightmost):
                    if reads[i][k] != '-' and reads[j][k] != '-':
                        r_notgap += 1
                rightmost += -1

            seq1 = reads[i][leftmost:rightmost]
            seq2 = reads[j][leftmost:rightmost]
            #print('i', i, 'j', j, 'leftmost:rightmost', leftmost, len(reads[i])-rightmost)

            if i != j:
                y_l = [0]
                cnt = 0
                for k in range(len(seq1)):
                    if seq1[k] == '-' and seq2[k] == '-':
                        cnt += 0
                    elif seq1[k] == '-' and seq2[k] in ['a', 't', 'g', 'c']:
                        cnt += -1
                    elif seq2[k] == '-' and seq1[k] in ['a', 't', 'g', 'c']:
                        cnt += -1
                    elif seq1[k] == seq2[k]:
                        cnt += 1
                    elif seq1[k] != seq2[k]:
                        cnt += -1
                    y_l.append(cnt)
    
                SM[i][j] = y_l
    
                #x = np.array([i for i in range(-1, len(seq1))])
                #y = np.array(y_l)
                #plt.figure()
                #plt.plot(x, y)
                #plt.title(f'{name1} and {name2} scoring')
                #plt.xlabel('position')
                #plt.ylabel('pair score')
                #plt.savefig(f'{sys.argv[2]}/{i}and{j}.png')
    
    #SM: Score Matrix
    SM = np.array(SM)
    #print('ScoreMatrix\n', SM)
    
    
    

    '''///// 判定（ヘテロ有無/アセンブリエラー・マッピングエラー） /////'''
    #length = len(SM[0][1])
    TM = np.zeros((num, num)) #Typing Matrix: リード数×リード数, リード関係を記録
    haba = 200
    
    # ---関数---
    # 同タイプ or 異タイプ判定
    # 初めにEnd scoreをみて、これが小さすぎる場合 --> アセンブリ・マッピングエラーと判断 --> 0
    # 次に全体のscoreをみて、同タイプの場合 --> 1, 異タイプの場合 --> 2
    def Judge(i, j, scores, length):
        if scores[-1] < 150:
            return 0
        else:
            cnt = 0
            cnt_decrease = 0
            while haba*(cnt) < length:
                if haba*(cnt+1) < length:
                    diff = scores[haba*(cnt+1)] - scores[haba*cnt]
                else:
                    diff = scores[-1] - scores[haba*cnt]
    
                if diff < 0:
                    cnt_decrease += 1
                cnt += 1
            if cnt_decrease >= round(int(sys.argv[4])/1000):    # haba 塩基区切りで、端の差が負である小ブロックが複数個(block size / 1000)あれば、そのブロックは異タイプと判断
                return 2
            else:                   #それ以外であれば、そのブロックは同タイプと判断
                return 1
    
    
    
    # ---関数---
    # リード断片の分類　ハプロタイプ1 / ハプロタイプ2 / エラー
    def Classifier(tm):
        read_groups = []                        #リード数（行/列数）だけ各々仲間を集める
        for i in range(num):                                #|
            group_i = [i]                                   #|
            for j in range(i, num):                         #|
                if i != j:                                  #|
                    if tm[i][j] == 1:                       #|
                        group_i.append(j)                   #|
            for j in range(0, i):                           #|
                if tm[j][i] == 1:                           #|
                    group_i.append(j)                       #|
            read_groups.append(set(group_i))                #|
                                                            #\/
        meta_groups = []                        #縮退させて、グループを同定する（アセンブリエラー・マッピングエラーは少数と仮定）
        for i in range(num):
            if not read_groups[i] in meta_groups:
                meta_groups.append(read_groups[i])
    
        meta_groups2 = []
        for group in meta_groups:
            group_copied = group
            for group2 in meta_groups:
                if group2 != group:
                    for r in group2:
                        if r in group:
                            for r in group2:
                                group_copied.add(r)
            meta_groups2.append(group_copied) 
    
        meta_groups = []
        for group in meta_groups2:
            if not group in meta_groups:
                meta_groups.append(group)
                                    
        return meta_groups
    
    
    
    
    '''///// タイプ行列 /////'''
    for i in range(num):
        for j in range(i, num):
            if i != j:
                TM[i][j] = Judge(i, j, SM[i][j], len(SM[i][j]))
    
    meta_groups = Classifier(TM)
    print('read groups\n', meta_groups)
    
    # この後、分類したリード数の割合でヘテロ合格 or 不合格を判定 30% ~ で採用
    errors = set({})
    haplos = []
    for group in meta_groups:
        if len(group) / num >= 0.3:
            haplos.append(group)
        else:
            for i in group:
                errors.add(i)

    #リード本数の確率を計算して閾値ならホモ判定・以上ならヘテロ判定
    #保留　10/14

    if haplos == []:
        haplos.append(tuple(errors))   #メインの方が検出できないときはerrorsリードをhaplosとする
        errors = set({})

    if len(haplos) > 2:
        print('エラー　推定ヘテロタイプが３種類以上あります')
    
    elif len(haplos) == 1:
        haplo1 = haplos[0]
        print('ホモ接合（推定）\n', 'リード番号\t', haplo1)
        print('エラーリード\t', errors)
    
    elif len(haplos) == 2:
        haplo1 = haplos[0]
        haplo2 = haplos[1]
        print('ヘテロ接合（推定）\n', 'リード番号\t', 'haplo 1', haplo1, '\t', 'haplo 2', haplo2)
        print('エラーリード\t', errors)
    
    
    
    #リードをチーム分けした結果をfasta形式で保存
    name = sys.argv[1]
    print('name', name)
    name = os.path.basename(name)
    name = name[:-9]
    
    if len(haplos) == 1:
        if len(haplo1) == 1:        #リードが１本しかない時はout4（サブコンセンサスoutput）に飛ぶ
            with open(f'{sys.argv[3]}/{name}_out4.txt', mode='w') as f:
                s = ''
                for i in haplo1:
                    s += f'>{name}\t1\n'
                    read_i = reads[i].replace('-', '')
                    s += f'{read_i}\n'
                f.write(s)
        else:
            with open(f'{sys.argv[3]}/{name}_out7.fa', mode='w') as f:
                s = ''
                for i in haplo1:
                    s += f'>{names[i]}\n'
                    read_i = reads[i].replace('-', '')
                    s += f'{read_i}\n'
                f.write(s)
    
    elif len(haplos) == 2:
        if len(haplo1) == 1:        #リードが１本しかない時はout4（サブコンセンサスoutput）に飛ぶ
            with open(f'{sys.argv[3]}/{name}_h1_out4.txt', mode='w') as f:
                s = ''
                for i in haplo1:
                    s += f'>{name}\t1\n'
                    read_i = reads[i].replace('-', '')
                    s += f'{read_i}\n'
                f.write(s)
        else:
            with open(f'{sys.argv[3]}/{name}_h1_out7.fa', mode='w') as f:
                s = ''
                for i in haplo1:
                    s += f'>{names[i]}\n'
                    read_i = reads[i].replace('-', '')
                    s += f'{read_i}\n'
                f.write(s)
    
        if len(haplo2) == 1:
            with open(f'{sys.argv[3]}/{name}_h2_out4.txt', mode='w') as f:
                s = ''
                for i in haplo2:
                    s += f'>{name}\t1\n'
                    read_i = reads[i].replace('-', '')
                    s += f'{read_i}\n'
                f.write(s)
        else:
            with open(f'{sys.argv[3]}/{name}_h2_out7.fa', mode='w') as f:
                s = ''
                for i in haplo2:
                    s += f'>{names[i]}\n'
                    read_i = reads[i].replace('-', '')
                    s += f'{read_i}\n'
                f.write(s)

