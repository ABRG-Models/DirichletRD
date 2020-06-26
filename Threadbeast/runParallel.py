import sys
import time
from math import *
import numpy as np
import os
import concurrent.futures
os.system("cp centres.h logs")
for j in range(0,8):
    valTemp = 3.6 - j*0.4
    Dn = exp(valTemp)
    Dc = Dn*0.3
    for k in range(0,4):
        Chi = (0.4 + k*0.2)*Dn
        logpath = "logs/Dn" + str(j) + "Chi" + str(k)
        print(logpath)
        command1 = "mkdir " + logpath
        os.system(command1)
        command="sim/build/processNoDisplay " + logpath +" 0.0001 " + str(Dn)  + " " +  str(Chi) + " " + str(Dc) +  " 8 " + " 6.0 " + " 400000 " + " 0 > " + logpath + "/output &"
        with concurrent.futures.ProcessPoolExecutor(max_workers=32) as executor:
          executor.submit(os.system(command))
        #  executor.submit(command)
print("sript finished")

