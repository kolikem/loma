'''///////////////////////////////////////////////////////////////////////////////
    ・out4 (ConsMaker.py output) を座標順に並び替える
    ・ヘテロチェックを行い
        ヘテロブロックが存在しない場合 --> 名前を連番に変更
        ヘテロブロックが存在する場合   --> ヘテロブロックが連続していたら配列決定
///////////////////////////////////////////////////////////////////////////////'''

print("---sort_out4_haplo.py---")

# sys.argv[1] - out4.txtのdirectory
# sys.argv[2] - NA?????_chr??_?????-?????
# sys.argv[3] - EsS.pyのステップサイズ


import glob
import sys
import re
import os

name = sys.argv[2]
out4_path = glob.glob(sys.argv[1]+f'/*{name}*_out4.txt')
# 例 /Users/ikemotoko/Desktop/100p/---/NA18943_chr1_101119569-101159569_-1337-1662_out4.txt 
# 例 /Users/ikemotoko/Desktop/100p/---/NA18943_chr1_84467176-84487176_27384-28883_h1_out4.txt

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
        #m = re.search(r'/(NA|GM).*chr.*-.*out4', out4)
        #num = m.group()
        num = os.path.basename(out4) #
        num = num[:-9] #
        #num = num[1:-5]
        m = re.search(r'\d*-\d*_.*-.*', num)
        num = m.group()
        m = re.search(r'-.*_.*-.*', num)
        num = m.group()
        m = re.search(r'_.*', num)
        num = m.group()
        num = num[1:]

        m = re.search(r'.*\d-', num)
        numL = m.group()
        numL = numL[:-1]

        if '--' in num:
            m = re.search(r'--.*', num)
            numR = m.group()
            numR = numR[1:]
        else:
            m = re.search(r'\d-.*', num)
            numR = m.group()
            numR = numR[2:]

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

#ヘテロブロックが存在する場合
elif het_num != 0:
    number_l = []
    number_r = []
    for out4 in out4_path:
        #print('out4', out4)
        m = re.search('/' + sys.argv[2] + r'.*-.*out4', out4)
        num = m.group()
        num = num[1:]
        
        if '_typed_' in out4:
            m = re.search(r'-?[0-9]*--?[0-9]*_h._typed_', num)
            num = m.group()
            num = num[:-10]

        else:
            m = re.search(r'-?[0-9]*--?[0-9]*_out4', num)
            num = m.group()
            num = num[:-5]

        m = re.search(r'-?[0-9]*-', num)
        numL = m.group()
        numL = numL[:-1]

        if '--' in num:
            m = re.search(r'--[0-9]*', num)
            numR = m.group()
            numR = numR[1:]
        else:
            m = re.search(r'[0-9]-[0-9]*', num)
            numR = m.group()
            numR = numR[2:]

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


    #print('number_l_sorted', number_l_sorted)
    #print('out4_path', out4_path)

    #名前の書き換え
    cnt = [0]
    for i in range(1, len(number_l_sorted)):
        if number_l_sorted[i-1] == number_l_sorted[i]:
            cnt.append(cnt[i-1])
        else:
            cnt.append(cnt[i-1]+1)

    for i in range(len(number_l_sorted)):
        num_l = number_l_sorted[i]
        num_r = number_r_sorted[i]
        for out4 in out4_path:
            if f'_{num_l}-{num_r}_' in out4:
                os.rename(out4, out4[:-4] + '_' + str(cnt[i]) + '.txt')
                out4_path.remove(out4)
                break 

