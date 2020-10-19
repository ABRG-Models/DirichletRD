rm  ./Morphlogs/*.data ./Morphlogs/*.Vect ./Morphlogs/*.txt ./Morphlogs/*.out
# arg1 logpath
# arg2 dt
# arg3 Dn
# arg4 Dchi
# arg5 Dc
# arg6 scale
# arg7 x extent
# arg8 number of time steps
# arg9 printing interval
# arg10 lstart true means start from previous run
# arg11 Lfixedseed if true we are using a fixed seed
# arg12 Lgraphics if true we use graphics
./build/processField /home/john/Neuroscience/DirichletRD/Threadbeast/sim/Morphlogs/ 0.0001 36.0 104.0 12.0 8 5.0 100 99 1 0 1

