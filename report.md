## Report of bioinformatic analysis of: Aedes_Recife_R1.fastq
1. Dataset retrieved from thoth01 server at /home/data_for_classes/2022_MIP_280A4/final_project_datasets  
```
cp Aedes_Recife_R1.fastq ~/2022_MIP_280A4_final_project
```

2. Activate conda environment
```
conda activate bio_tools
```

3. Dataset analyzed with Fastqc v0.11  
```
fastqc Aedes_Recife_R1.fastq
```

4. Review output here: [Aedes_Recife_R1.fastq FastQC Report.pdf](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167081/Aedes_Recife_R1.fastq.FastQC.Report.pdf)

5. Adapter sequences and low quality bases trimmed using Cutadapt v3.5  
```
cutadapt \  
  -a AGATCGGAAGAG \  
  -q 30,30 \  
  --minimum-length 80 \  
  -o Aedes_Recife_R1_trimmed2.fastq \  
  Aedes_Recife_R1.fastq \  
  | tee cutadapt2.log
```  
Review output here: [cutadapt2.log](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167019/cutadapt2.log)

6. Dataset again analyzed with Fastqc v0.11
```
fastqc Aedes_Recife_R1_trimmed2.fastq
```

7. Review output here: [Aedes_Recife_R1_trimmed2.fastq FastQC Report.pdf](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167093/Aedes_Recife_R1_trimmed2.fastq.FastQC.Report.pdf)

8. Download *Aedes aegypti* genome from RefSeq database into directory on thoth01 server
```
curl -OL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz
```

9. Create index of *A. aegypti* genome to filter out host reads using Bowtie2 v2.4
```
bowtie2-build GCF_002204515.2_AaegL5.0_genomic.fna.gz aedes_aeg_index
```

10. Map reads against *A. aegypti* genome and retain reads that did not map in new fastq file 
```
bowtie2 -x aedes_aeg_index \
  --local \
  -U Aedes_Recife_R1_trimmed2.fastq \
  --no-unal \
  --threads 24 \
  -S aedes_recife_mapped_to_aegypti_genome.sam \
  --un unmapped_aedes_reads.fastq
```
Review output here: [Mapping_against_aedes_output.txt](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167014/Mapping_against_aedes_output.txt)

11. Assemble non-mapped reads using Spades v3.15
```
spades.py -o non_mapped_reads_assembly \
  -s unmapped_aedes_reads.fastq \
  -m 24 -t 18
```
Review output here: [contigs_from_unmapped_reads.txt](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167054/contigs_from_unmapped_reads.txt)

### To validate and further characterized assembled contigs:

12. Build new index from assembled contigs using Bowtie2 v2.4
```
bowtie2-build contigs.fasta contigs_index
```

13. Map previously unmapped reads against index built from contigs
```
bowtie2 -x contigs_index --local -U unmapped_aedes_reads.fastq --no-unal --threads 24 -S unmapped_reads_validat
ion.sam
```
Review output here: [Mapping_against_contigs_output.txt](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167190/Mapping_against_contigs_output.txt)
