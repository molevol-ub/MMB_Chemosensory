makeblastdb -in DB/GR_db.fasta -dbtype prot
makeblastdb -in DB/OR_db.fasta -dbtype prot

ANNOTATIONFASTA=Annotations/GR_genomic_and_annotated_proteins_trimmed.fasta

blastp -query $ANNOTATIONFASTA -db DB/GR_db.fasta -outfmt "6 std qlen slen" -out "$ANNOTATIONFASTA"_GRblast.txt -num_threads 1 -max_target_seqs 5
blastp -query $ANNOTATIONFASTA -db DB/OR_db.fasta -outfmt "6 std qlen slen" -out "$ANNOTATIONFASTA"_ORblast.txt -num_threads 1 -max_target_seqs 5

makeblastdb -in DB/GR_Dmel_db.fasta -dbtype prot
blastp -query $ANNOTATIONFASTA -db DB/GR_Dmel_db.fasta -outfmt "6 std qlen slen" -out "$ANNOTATIONFASTA"_GRdmelblast.txt -num_threads 1 -max_target_seqs 5


hmmsearch -o "$ANNOTATIONFASTA"_hmmer.out --notextw --tblout "$ANNOTATIONFASTA"_hmmer.tblout --domtblout "$ANNOTATIONFASTA"_hmmer.domtblout DB/GR_db.hmm $ANNOTATIONFASTA


# Edit the input files in the script (see lines 13 to 30) and run:
perl get_annotation_classification.pl
