rm ../logs/*.txt ../logs/*.out ../logs/*.data ../logs/*.Vect
#arg1 logpath
#arg2 dt timestep
#arg3 Dn
#arg4 Dchi
#arg5 Dc
#arg6 scale
#arg7 x extent
#arg8 numsteps
#arg9 Lcontinue, if true load initial field from file
#arg10 phaseShift, if true we shift all regions
./build/processCosines  ../logs 0.0001 36.0 18.0 10.8 8 5.0  10 0 0
