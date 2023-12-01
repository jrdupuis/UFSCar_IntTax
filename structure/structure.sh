qsub_folder=qsubs
output_folder=StructureResults

mkdir $qsub_folder
mkdir $output_folder
mkdir err
mkdir out

for k in {1..4}; 
do 
for r in {1..5}; 
do
echo k_$k.rep_$r 
cat structure_qsub_header > $qsub_folder/k_$k.rep_$r.qsub
echo "structure -D $RANDOM -K $k -o $output_folder/$k.$r.output" >> $qsub_folder/k_$k.rep_$r.qsub
sleep 1
sbatch $qsub_folder/k_$k.rep_$r.qsub &
done
sleep 1
done

