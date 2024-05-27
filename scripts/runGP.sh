#!/bin/bash
## This script is invoked as:
##doparallel "runGP.sh {1} $molecule $cwd $nojf $nchan" "$(seq $noj1 $nojf)"
##The first job is the tors scan
batchn=$(($1-1))
if [ $1 -eq 1 ];then
   batch=torsion
else
   batch=batch$batchn
fi
#Make batch directory
if [ -d ${batch} ]; then  rm -r ${batch} ; fi
mkdir ${batch}
#copy neccesary files.
cp amk.dat $2.xyz ${batch}
echo "tsdirll $3/tsdirLL_$2"  >> ${batch}/amk.dat
##copy frags only if sampling is vdw
copyfrags=$(awk 'BEGIN{cf=0};{if($1=="sampling" && $2=="vdW") cf=1};END{print cf}'  amk.dat)
if [ $copyfrags -eq 1 ]; then
   for i in $(seq 1 "$number_of_fragments"); do
      frag=( "${frag[@]}"  "$(awk -v i="$i" '{if($1=="fragment_'$i'"){print $2;exit}}' "$inputfile")" )
   done
   for i in $(seq 0 $((number_of_fragments-1))); do
      cp "${frag[$i]}".xyz ${batch}
   done
fi
##Launch the calcs
if [ $1 -eq 1 ]; then
   (cd ${batch} && tors.sh > scan_tors.log)
elif [ $5 -eq 0 ]; then
   (cd ${batch} && amk.sh amk.dat >amk.log)
else
   (cd ${batch} && LocateTs.py amk.dat $((batchn-1)) 1 $6 > amk.log)
fi
