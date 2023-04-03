'''/////////////////////////////////////////////////////////////////////////////////
        閾値以下の参考リード数のみで構成されているサブコンセンサスは切る
/////////////////////////////////////////////////////////////////////////////////'''

import sys
import glob
import re
import os

# sys.argv[1] - 領域名
# sys.argv[2] - out4 ディレクトリ
# sys.argv[3] - forward 実行 or backward 実行; 'F' or 'B'
# sys.argv[4] - __本未満の参考リード数で端を切る
# sys.argv[5] - ヘテロ分類ブロックあり：１、なし：０

thresh = int(sys.argv[4])

sys2 = sys.argv[2]
if sys2[-1] == '/':
    sys2 = sys2[:-1]

path_4 = sys2 + '/' + sys.argv[1]

# glob 関数: 条件にマッチするファイル名を列挙する
file_l = glob.glob(path_4 + '*out4_*')
#print(file_l)

# 閾値１０(thresh本)リードで端のブロックを切り落とす
N = len(file_l)


##### ヘテロ分類ブロックない場合 #####

if int(sys.argv[5]) == 0:
    # Forward
    if sys.argv[3] == 'F':
        cnt_forward = 0
        stop_i = 0
        for i in range(N):
            stop_j = 0
            for j in range(N):
                if not f'out4_{i}.txt' in file_l[j]:
                    continue
                else:
                    stop_j = 1
                    f = open(file_l[j])
                    for line in f:
                        if line[0] == '>':
                            num = re.search('\t[0-9]*\n', line)
                            num = num.group()
                            num = int(num[1:-1])
                            if num < thresh:
                                cnt_forward += 1
                            elif num >= thresh:
                                stop_i = 1
                            break
                    f.close()
                if stop_j == 1:
                    break
            if stop_i == 1:
                break
        print(cnt_forward)

    # Backward
    elif sys.argv[3] == 'B':
        cnt_backward = N - 1
        stop_i = 0
        for i in range(N):
            stop_j = 0
            for j in range(N):
                if not f'out4_{N-1-i}.txt' in file_l[j]:
                    continue
                else:
                    stop_j = 1
                    f = open(file_l[j])
                    for line in f:
                        if line[0] == '>':
                            num = re.search('\t[0-9]*\n', line)
                            num = num.group()
                            num = int(num[1:-1])
                            if num < thresh:
                                cnt_backward -= 1
                            elif num >= thresh:
                                stop_i = 1
                            break
                    f.close()
                if stop_j == 1:
                    break
            if stop_i == 1:
                break
        print(cnt_backward)




##### ヘテロ分類ブロックある場合 #####

elif int(sys.argv[5]) == 1:
    f_sorted = []
    for i in range(len(file_l)):
        f = [file_l[i]]
        if '_typed_' in file_l[i]:
            f.append('o')   # ヘテロ分類ブロックの場合: f[1] = 'o'
        else:
            f.append('x')   # ヘテロ分類ブロックの場合: f[1] = 'x'
        m = re.search(r'out4_[0-9]*', file_l[i])
        num = m.group()
        num = int(num[5:])  # f[2] : ブロック番号
        f.append(num)
        f_sorted.append(f)
    f_sorted = sorted(f_sorted, key=lambda x: x[2])
    #print('f_sorted', f_sorted)


    # Forward
    if sys.argv[3] == 'F':
        stop = 0
        for i in range(len(f_sorted)):
            if f_sorted[i][1] == 'x':   #ヘテロ分類ではないブロック
                with open(f_sorted[i][0]) as f:
                    for line in f:
                        if line[0] == '>':
                            num = re.search('\t[0-9]*\n', line)
                            num = num.group()
                            num = int(num[1:-1])
                            if num >= thresh:
                                print(f_sorted[i][2])       #このブロック番号をprintしてコマンド変数に引き継ぐ
                                stop += 1
                                break
                            else:
                                break
            elif f_sorted[i][1] == 'o': #ヘテロ分類であるブロック
                continue
            
            if stop == 1:
                break


    # Backward
    elif sys.argv[3] == 'B':
        stop = 0
        for i in range(len(f_sorted)):
            if f_sorted[N-1-i][1] == 'x':   #ヘテロ分類ではないブロック
                with open(f_sorted[N-1-i][0]) as f:
                    for line in f:
                        if line[0] == '>':
                            num = re.search('\t[0-9]*\n', line)
                            num = num.group()
                            num = int(num[1:-1])
                            if num >= thresh:
                                print(f_sorted[N-1-i][2])       #このブロック番号をprintしてコマンド変数に引き継ぐ
                                stop += 1
                                break
                            else:
                                break
            elif f_sorted[N-1-i][1] == 'o': #ヘテロ分類であるブロック
                continue
            
            if stop == 1:
                break


# コマンドにこの実行結果（cnt_forward と cnt_backward）を引き渡して変数にする。
# cnt_FB=`python3 CutOut4.py ....`
