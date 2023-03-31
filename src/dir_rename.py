'''h1/h2ディレクトリ内の名前の統一'''

print("---dir_rename.py---")

import sys
import glob
import re
import os

# sys.argv[1] - ディレクトリ名
# sys.argv[2] - h1 --> 1, h2 --> 2

direc = sys.argv[1]
if direc[-1] == '/':
    direc = direc[:-1]

file_dir_name = glob.glob(direc + '/' + f'*_h{sys.argv[2]}_typed_*out4*')

for name in file_dir_name:
    name0 = name.split('/')[0]
    name1 = name.split('/')[1]

    mae = re.search(r'.*_h' + sys.argv[2] + '_', name1)
    mae = mae.group()
    mae = mae[:-3]
    #print('mae', mae)

    ushiro = re.search(r'out4_.*', name1)
    ushiro = ushiro.group()
    #print('ushiro', ushiro)

    new_name = mae + ushiro
    #print(new_name)

    os.rename(name, name0 + '/' + new_name)
