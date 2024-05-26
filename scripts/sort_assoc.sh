inputfile=$1
#specA=$(awk '{if($1=="fragmentA") print $2}' $inputfile)
#specB=$(awk '{if($1=="fragmentB") print $2}' $inputfile)
#specC=$(awk '{if($1=="fragmentC") print $2}' $inputfile)
#assocdir=${PWD}/assoc_${specA}_${specB}_${specC}


concat_fr=${frag[0]}
for i in $(seq 1 "$((number_of_fragments-1))"); do
   concat_fr=${concat_fr}_${frag[i]}
done
echo "${concat_fr}"
assocdir=${cwd}/assoc_${concat_fr}

# We start with the TSs
if [ -f tmp_oren ]; then rm tmp_oren ; fi
if [ -f tmp_nonoren ]; then rm tmp_nonoren ; fi

#echo "non-ordered energies"
for i in $(ls $assocdir/assoc*out)
do
  if [ "$2" = "mopac" ]; then
     awk '/FINAL HEAT OF FORMATION =/{e=$6};END{print "'$i'",e}' $i >> tmp_nonoren
  elif [ "$2" = "qcore" ]; then
     awk '/Energy=/{e0=$2};END{print "'$i'",e0}' $i >> tmp_nonoren
  fi
done
file=tmp_nonoren
#echo "ordered energies"
awk '{a[++d]=$2}
END{
n = asort(a,b)
for (i=1; i<=n; i++){
 print b[i]
 }
}' $file > tmp_oren
#echo "final energies"
paste tmp_oren tmp_nonoren | awk '{oe[NR]=$1;flag[NR]=$2;e[NR]=$3}
END{
i=1
while(i<=NR){
  j=1
  while(j<=NR){
    if(oe[i]==e[j]) {print "Structure",i,flag[j],e[j];break}
    ++j
    }
  ++i
  }
}' > $assocdir/assoclist_sorted

