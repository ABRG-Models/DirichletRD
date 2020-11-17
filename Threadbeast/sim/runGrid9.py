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
        os.system("mkdir " + logpath);
        for m in range(0,21):
            reg = "./logs/" + str(m)
            os.system("cp ./logs/" + str(m) + ".h5 " + logpath);
            command="$HOME/Neuroscience/DirichletRD/Threadbeast/sim/build/processCollect " + logpath + " 0.0001 " + str(Dn)  + " " +  str(Chi) + " " + str(Dc) + " 8 " + " 5.0 " +  " 200000" + " 199999 " + " 0 " + " 0 " + " 1 " + " 0  > " + logpath + "/output "
            print(logpath)
            print(command)
            os.system(command)
            time.sleep(10)
print("script finished");
