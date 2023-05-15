for sample in R2826 R2835 R2839 R2845 R2855 R2857 R2869 R2874 R2894 R2895 ; do
    echo python Step_1_splam_prediction.py -f ../results/polyA/$sample/fasta/junction.fa -o ../results/polyA/$sample/junction_score.bed -m ../src/MODEL/SPLAM_v11/splam_14.pt
    python Step_1_splam_prediction.py -f ../results/polyA/$sample/fasta/junction.fa -o ../results/polyA/$sample/junction_score.bed -m ../src/MODEL/SPLAM_v11/splam_14.pt
done
