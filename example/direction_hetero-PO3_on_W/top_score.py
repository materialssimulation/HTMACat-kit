import numpy as np

num_top = 3 # 取前几个构型



fcontent = ''
with open('score_log.txt', 'r') as f:
    fcontent = f.read()
fcontent = fcontent.split('\n')



scores_l = []
for i,line in enumerate(fcontent):
    if ' | ' in line:
        score = float( line.split('=')[1].split()[0] )
        scores_l.append(score)
print('Configurations with top scores:')
for i in range(num_top):
    print( np.argsort(scores_l)[::-1][i],  scores_l[np.argsort(scores_l)[::-1][i]])



if 1: # 如需生成计算文件夹
    import os, shutil
    for i in range(num_top):
        sidx = str( np.argsort(scores_l)[::-1][i] )
        folder_name = sidx + '_'
        dir_path = os.path.join(os.getcwd(), folder_name)
        if not os.path.isdir(dir_path):
            vasp_files = [f for f in os.listdir('.') if f.endswith(sidx+'.vasp')]
            os.mkdir(dir_path)
            shutil.copy(vasp_files[0], os.path.join(dir_path, 'POSCAR'))
