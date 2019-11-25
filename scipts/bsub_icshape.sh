

for file in $1/*fastq
# for file in $1/*trimmed.fastq
do

echo $file

bsub -q Z-HNODE -oo ${file}.out -eo ${file}.err python bsub_icshape.py ${file}
sleep 10

done
