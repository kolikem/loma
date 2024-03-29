'''///////////////////////////////////////////////////////////////////////////////
    ・out4 (ConsMaker.py output) を座標順に並び替える
    ・ヘテロチェックを行い
        ヘテロブロックが存在しない場合 --> 名前を連番に変更
        ヘテロブロックが存在する場合   --> ヘテロブロックが連続していたら配列決定
///////////////////////////////////////////////////////////////////////////////'''

print("---sort_out4.py---")

# sys.argv[1] - out4.txtのdirectory
# sys.argv[2] - NA?????_chr??_?????-?????
# sys.argv[3] - EsS.pyのステップサイズ


import glob
import sys
import re
import os

name = sys.argv[2]
out4_path = glob.glob(sys.argv[1]+f'/*{name}*_out4.txt')
# 例 /Users/ikemotoko/Desktop/100points_consensus/20kbp_chr1_101119569-101159569_2/NA18943_chr1_101119569-101159569_-1337-1662_out4.txt 
# 例 /Users/ikemotoko/Desktop/10p/chr1_844/NA18943_chr1_84467176-84487176_27384-28883_h1_out4.txt

het_num = 0
for out4 in out4_path:
    if '_h1_' in out4:
        het_num += 1
print('ヘテロブロック数', het_num)



# ヘテロブロック存在しなければ普通に名前を変更
if het_num == 0:
    number_l = []
    number_r = []
    for out4 in out4_path:
        num = os.path.basename(out4)
        print("file name:", num)
        num = num[:-9]
        num = num.split("_")[-1]    # "-1337-1662", "-10000--5000"
        if num[0] == "-":
            if '--' in num:
                num = num.split("-")
                numL = "-"+num[1]
                numR = "-"+num[3]
            else:
                num = num.split("-")
                numL = "-"+num[1]
                numR = num[2]
        else:
            num = num.split("-")
            numL = num[0]
            numR = num[1]

        print("numL", numL)
        number_l.append(int(numL))
        number_r.append(int(numR))

    number_l_sorted = sorted(number_l)
    number_r_sorted = sorted(number_r)


    #ブロック連続が途切れるところを記録（右側ブロック番号）
    stepSize = int(sys.argv[3])
    block_gap_l = []
    for i in range(len(number_l_sorted)-1):
        if number_l_sorted[i+1] !=  number_l_sorted[i] + stepSize:
            block_gap_l.append(i+1)

    block_gap_l.append(len(number_l_sorted))

    #名前の書き換え

    for i in range(len(number_l_sorted)):
        for out4 in out4_path:
            if f'_{number_l_sorted[i]}-{number_r_sorted[i]}_' in out4:
                os.rename(out4, out4[:-4]+'_'+str(i)+'.txt')
                break         

