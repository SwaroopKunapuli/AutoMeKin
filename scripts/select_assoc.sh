#!/bin/bash
source utils.sh
sharedir=${AMK}/share
###
elements=${sharedir}/elements

cwd=$PWD
inputfile=$1
##reading inputfile
read_input
##

concat_fr=${frag[0]}
for i in $(seq 1 "$((number_of_fragments-1))"); do
   concat_fr=${concat_fr}-${frag[i]}
done
echo "${concat_fr}"
assocdir=${cwd}/assoc_${concat_fr}


rm -f $assocdir/selected_min* fort.*

for i in $(seq 0 "$((number_of_fragments-1))"); do
   createMat.py "${frag[$i]}".xyz 1
   cp ConnMat ConnMat"${frag[$i]}"
done

n=0
nmin="$(ls $assocdir/assoc*.out | wc -l | awk '{print $1}')"
#echo "A total of $nmin minima have been optimized"
for i in $(ls $assocdir/assoc*.out )
do
   ((n=n+1))
   if [ "$program_opt" = "mopac" ]; then
      get_geom_mopac.sh $i | awk '{if(NF==4) print $0}' >mingeom0
   elif [ "$program_opt" = "qcore" ]; then
      awk '/Final structure/{flag=1; next} EOF{flag=0} flag' $i >mingeom0
   fi
   if [ $n -eq 1 ]; then
      met_label=$(awk 'NR==FNR{l[NR]=$1;tne=NR}
      NR>FNR{IGNORECASE = 1
         for(i=1;i<=tne;i++){
            if( $1 == l[i] && i==13) print FNR
            if( $1 == l[i] && i>=21 && i<=30) print FNR
            if( $1 == l[i] && i>=39 && i<=48) print FNR
            if( $1 == l[i] && i>=72 && i<=80) print FNR
            }
      }' $elements mingeom0)
      if [ -z "$met_label" ]; then met_label=0 ; fi
   fi
   

   currentline=0
   total_diff=0
   for j in $(seq 0 "$((number_of_fragments-1))"); do
      echo "${natomfr[$j]}" > mingeom"${frag[$j]}"
      echo ''>> mingeom"${frag[$j]}"
      futureline=$((currentline+${natomfr[$j]}))
      echo $futureline
      awk -v currentline=$currentline -v futureline=$futureline -v frag="${frag[$j]}" '{if(NR>currentline && NR<=futureline) print $0 >> "mingeom'$frag'"}' mingeom0
      currentline=$futureline
      creatMat.py mingeom"${frag[$j]}" 1
      cp ConnMat ConnMat"${frag[$j]}"p
      diff"${frag[$j]}"="$(paste ConnMat"${frag[$j]}" ConnMat"${frag[i]}"p | awk 'BEGIN{diff=0;natom='${natomfr[$j]}'} 
      {for(k=1;k<=natom;k++){d=$k-$(k+natom);diff+=sqrt(d*d)}
      }
      END{print diff}')"
      total_diff=$((total_diff+diff"${frag[$j]}"))
   done
  

   val=$(awk 'BEGIN{bo=0};{if('$met_label'==0) {print "0";exit}}
      /BOND ORDERS/{bo=1}
      {if(bo==1 && $1=='$met_label') {print $3;exit} }' $i | sed 's@(@@;s@)@@')
   if [ "$program_opt" = "mopac" ];then
      e=$(awk '/FINAL HEAT OF FORMATION/{print $6;exit}' $i)
   elif [ "$program_opt" = "qcore" ];then
      e=$(awk '/Energy=/{e0=$2};END{print e0}' $i)
   fi
   echo $i $val $e>> $assocdir/selected_min_$total_diff
done

for i in $(seq 0 10)
do
   if [ ! -f $assocdir/selected_min_$i ]; then
      continue 
   else
      if [ $i -eq 0 ]; then
         echo "Structures found with no changes in the geometries of "${frag[@]}
      elif [ $i -eq 1 ]; then
         echo "Structures found with 1 change in the geometries of "${frag[@]}
      else
         echo "Structures found with $i changes in the geometries of "${frag[@]}
      fi
      awk 'BEGIN{min=10^10}
      {name[NR]=$1
      val[NR]=$2
      e[NR]=$3
      if($3<min) min=$3
      }
      END{i=1
      while(i<=NR){
         en=e[i]-min
         point=2^val[i]-0.1*en
         print name[i],point
         i++
         }
      }' $assocdir/selected_min_$i > $assocdir/selected_min.out
      break
   fi
done

smin=$(awk 'BEGIN{max=-10^10}
{if($2>max) {max=$2;name=$1 }}
END{print name }' $assocdir/selected_min.out)
sminprint=$(echo $smin | sed 's@'$cwd'/@@')
echo "Selected minimum output file   = $sminprint"
echo "Selected minimum XYZ structure = ${molecule}.xyz" 
if [ "$program_opt" = "mopac" ];then
   get_geom_mopac.sh $smin > ${molecule}.xyz
elif [ "$program_opt" = "qcore" ];then
   awk '/Final structure/{flag=1; next} EOF{flag=0} flag{++n;a[n]=$0};END{print n"\n";for(i=1;i<=n;i++) print a[i]}' ${smin} > ${molecule}.xyz
fi
