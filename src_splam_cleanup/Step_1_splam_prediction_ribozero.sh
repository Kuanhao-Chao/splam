for sample in R12258  R12260  R12263  R12265  R12266  R12277  R12278  R12280  R12285  R12287 ; do
    echo python Step_1_splam_prediction.py -f ../results/ribozero/$sample/fasta/junction.fa -o ../results/ribozero/$sample/junction_score.bed -m ../src/MODEL/SPLAM_v11/splam_14.pt
    python Step_1_splam_prediction.py -f ../results/ribozero/$sample/fasta/junction.fa -o ../results/ribozero/$sample/junction_score.bed -m ../src/MODEL/SPLAM_v11/splam_14.pt
done
