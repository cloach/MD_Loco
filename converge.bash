#! /bin/bash
rm out.temp
for i in `seq 1.2 0.01 10.0` 
do
java Loco atom_list4 1 $i > temp.dump
pairE=`grep "pairEnergyP" temp.dump | cut -d" " -f3`
embedE=`grep "ingEnergyP" temp.dump | cut -d" " -f3`
totalE=`grep "alEnergyP" temp.dump | cut -d" " -f3`
pairF=`grep "pairForceP" temp.dump | cut -d" " -f3`
#pairFS=`grep "1 0.0 0.0 0.0" temp.dump | cut -d" " -f5`
pairFS=`grep "magx1:" temp.dump | cut -d" " -f2`
embedF=`grep "ingForceP" temp.dump | cut -d" " -f3`
totalF=`grep "alForceP" temp.dump | cut -d" " -f3`
pairFM=`grep "PairForceVectorMag" temp.dump | cut -d" " -f3`
embedFM=`grep "totalEmbedVectorMag" temp.dump | cut -d" " -f3`
#echo $i $pairE $pairF $pairFM $pairFS >> out.temp 
echo $i $pairE $pairF $pairFM $pairFS $embedE $embedF $embedFM $totalE $totalF >> out.temp 
echo ""
done
