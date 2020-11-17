import sys
import time
from math import *
import numpy as np
import os
for j in range(0,8):
    valTemp = 3.6 -  j*0.4
    Dn = exp(valTemp)
    Dc = 0.3*Dn
    for k in range(1,2):
        Chi = (0.0 + k*0.4)*Dn
        logpath = "./logs/Dn" + str(j);
        os.system("cp centres.data " + logpath);
        print(logpath)
        print(command)
        os.system(command)
        time.sleep(10)
print("script finished");
