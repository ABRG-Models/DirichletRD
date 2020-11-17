import sys
import time
from math import *
import numpy as np
import os
for j in range(0,99):
    logpathin = 'logsStaticDn/Dn' + str(j)
    logpathout1 = 'logsStaticDn/Images'
    infilename1 = logpathin + '/nn0.png'
    outfilename1 = logpathout1 +  '/' + str(j) + 'nnField.png'
    command1 = 'cp ' + infilename1 + ' ' + outfilename1
    os.system(command1)
    print(command1)
print("script finished")
