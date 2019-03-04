import sys
import time
from math import *
import numpy as np
import os
for j in range(0,2):
    valTemp = 3.0 - j*0.1
    Dn = exp(valTemp)
    for k in range(0,2):
        Chi = (k+1)*1.0
        logpath = "logs/Dn" + str(j) + "Chi" + str(k)
        print(logpath)
        command="sim/build/process " + logpath +" 0.00001 " + str(Dn) +" 5000 " + str(Chi) +" 1000"
        os.system(command)
print("script finished")
