# get a list of commands to run when adding PM

import os

curDir = os.getcwd()
folder = 'processing'
filePath = os.path.join(curDir, folder)
print(filePath)
fileList = sort(os.listdir(filePath))
f = open(os.path.join(curDir, 'addAFstuff.txt'), 'w')
for fi in fileList:
    line = 'python addPM.py AF {}\n'.format(os.path.join(filePath,fi))
    f.write(line)
f.close()
