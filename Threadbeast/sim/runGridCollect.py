import sys
import time
from math import *
import numpy as np
import os
for j in range(0,8):
    valTemp = 3.6
    Dn =  exp(valTemp)
    Dc = 0.3*Dn
    logpath = "./logsStaticDn/Dn" + str(j);
    os.system("cp ./centres.inp " + logpath);
    for k in range(1,2):
        Chi = (0.0 + k*0.4)*Dn
        command="$HOME/Neuroscience/DirichletRD/Threadbeast/sim/build/processCollect " + logpath + " 0.0001 " + str(Dn)  + " " +  str(Chi) + " " + str(Dc) + " 8 " + " 5.0 " +  " 50 " + " 50 " + " 0 " + " 0 " + " 1 " + " 0  > " + logpath + "/output"
        print(logpath)
        print(command)
        os.system(command)
        time.sleep(10)
print("script finished");
