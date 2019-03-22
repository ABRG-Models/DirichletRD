import sys
import time
from math import *
import numpy as np
import os
import concurrent.futures
for j in range(0,3):
    valTemp = 3.0 - j*2.0
    Dn = exp(valTemp)
    for k in range(0,3):
        Chi = (k+1)*1.0
        logpath = "logs/Dn" + str(j) + "Chi" + str(k)
        # print(logpath)
        command="sim/build/process " + logpath +" 0.00001 " + str(Dn)  + " " +  str(Chi) + " 50 " + " 10 " + " 1 &"
        # print(command)
        #  os.system(command)
        with concurrent.futures.ProcessPoolExecutor(max_workers=32) as executor:
          executor.submit(os.system(command))
        #  executor.submit(command)
print("sript finished")

