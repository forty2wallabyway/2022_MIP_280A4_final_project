### A. Workflow of bioinformatic analysis of: Aedes_Recife_R1.fastq  

---

1. Retrieve sequence data from thoth01 server at /home/data_for_classes/2022_MIP_280A4/final_project_datasets:    
```
cp Aedes_Recife_R1.fastq ~/2022_MIP_280A4_final_project
```

2. Activate conda environment created with bio_tools.yaml:
```
conda activate bio_tools
```

3. Analyze sequence data with Fastqc v0.11:  
```
fastqc Aedes_Recife_R1.fastq
```

4. Identified increased variance (and decrease in overall trend) of base call quality beginning around base position 60. Also identified presence of 'Illumina Universal Adapter' content. Review output here: [Aedes_Recife_R1.fastq FastQC Report.pdf](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167081/Aedes_Recife_R1.fastq.FastQC.Report.pdf)

5. Trim adapter sequences and low quality bases using Cutadapt v3.5:  
```
cutadapt \  
  -a AGATCGGAAGAG \  
  -q 30,30 \  
  --minimum-length 80 \  
  -o Aedes_Recife_R1_trimmed2.fastq \  
  Aedes_Recife_R1.fastq \  
  | tee cutadapt2.log
```  
Identified 4,095,494 reads in total. Removing adapter sequences (5.9%) and low quality reads (12.7%) yielded a remaining 3,574,307 reads. Review output here: [cutadapt2.log](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167019/cutadapt2.log)

6. Analyze sequence data again with Fastqc v0.11:
```
fastqc Aedes_Recife_R1_trimmed2.fastq
```

7. Identified decreased variance (and increase in overall trend) of base call quality beginning around base position 85. Also identified absence of 'Illumina Universal Adapter' content. Review output here: [Aedes_Recife_R1_trimmed2.fastq FastQC Report.pdf](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167093/Aedes_Recife_R1_trimmed2.fastq.FastQC.Report.pdf)

8. Download *Aedes aegypti* genome from RefSeq database into working directory on thoth01 server:
```
curl -OL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz
```

9. Create index of *A. aegypti* genome to filter out host reads using Bowtie2 v2.4:
```
bowtie2-build GCF_002204515.2_AaegL5.0_genomic.fna.gz aedes_aeg_index
```

10. Map reads against *A. aegypti* genome and retain reads that did not map in new fastq file:
```
bowtie2 -x aedes_aeg_index \
  --local \
  -U Aedes_Recife_R1_trimmed2.fastq \
  --no-unal \
  --threads 24 \
  -S aedes_recife_mapped_to_aegypti_genome.sam \
  --un unmapped_aedes_reads.fastq
```
Identified 99.26% of reads mapped to *A. aegypti* genome, as expected. Review output here: [Mapping_against_aedes_output.txt](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167014/Mapping_against_aedes_output.txt)

11. Assemble remaining, non-mapped reads (0.74% of trimmed and host-filtered data) using Spades v3.15:
```
spades.py -o non_mapped_reads_assembly \
  -s unmapped_aedes_reads.fastq \
  -m 24 -t 18
```
Generated 135 contigs (nodes) varying in length from 229 to 6,781 bp. Review output here: [contigs_from_unmapped_reads.txt](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167054/contigs_from_unmapped_reads.txt)

---

### B. Validation and further characterization of assembled contigs:

12. Build new index from assembled contigs using Bowtie2 v2.4:
```
bowtie2-build contigs.fasta contigs_index
```

13. Map previously unmapped reads against index built from contigs:
```
bowtie2 -x contigs_index --local -U unmapped_aedes_reads.fastq --no-unal --threads 24 -S unmapped_reads_validat
ion.sam
```
Identified a lower-than-expected overall alignment rate of 66.69%. It is likely that reads that did not align were initially not assembled into contigs and therefore not included in the contigs_index build. These reads could be comprised of low quality or low coverage reads that could not be assembled using Spades. Review output here: [Mapping_against_contigs_output.txt](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10167190/Mapping_against_contigs_output.txt)

14. Using Geneious v2022.2.2, import unmapped_reads_validation.sam along with the contigs.fasta output from the Spades assembly. Sort contigs by 'The number of sequences in an alignment or set of sequences' filter. This will allow for identification of the **top 12** contigs with the greatest coverage depth. 

15. BLAST each of the 12 contigs identified in the previous step using the [blastn tool](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) on the NCBI website. 

Report of identified contigs, including average coverage depth, can be reviewed here: [final_contigs_table.pdf](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10169713/final_contigs_table.pdf)

---

### C. Assessment of pairwise identity and variant analysis of Yeltsovka tombus-like virus (YTLV) putative polyprotein from contig 2

16. Build index from Yeltsovka tombus-like virus isolate Koltsovo (MW251332.1) using Bowtie v2.4:
```
bowtie2-build YTLV.fasta YTLV_index
```
Note: *YTLV.fasta* was downloaded from NCBI and moved into the thoth01 working directory using secure file transfer protocol (sftp) 

17. Map previously unmapped reads against index built from YTLV.fasta:
```
bowtie2 -x YTLV_index --local -U unmapped_aedes_reads.fastq --no-unal --threads 24 -S YTLV_mapped_reads.sam
```

18. Using Geneious v2022.2.2, import YTLV_mapped_reads.sam along with YTLV.fasta. This step will provide a view of the reads in this dataset mapped directly against the closest NCBI assession available. It was identified that *genome* coverage in this case was relatively low when compared to reads mappeded against host-filtered contigs assembled from the same dataset.

19. While still using Geneious, conduct a MAFFT pairwise alignment (Align/Assemble > Pairwise Align > MAFFT Alignment) between YTLV.fasta and contig 2 from the Spades assembly. Options for Algorithm, Scoring matrix, and Gap open penalty can be left at default. Choose option "Automatically determine sequences' direction." 

20. Once again working from the shell terminal, locate YTLV_mapped_reads.sam and perform the following steps in order to produce a variant analysis output using Lofrequ v2.1:  

a. First convert the .sam file to a .bam file
```
samtools view -b YTLV_mapped_reads.sam > YTLV_mapped_reads.bam 
```

b. Then sort the .bam file
```
samtools sort -T tmp -O 'bam' YTLV_mapped_reads.bam > YTLV_mapped_reads.sorted.bam 
```

21. Run Lofreq v2.1 as follows:
```
lofreq call -f YTLV.fasta YTLV_mapped_reads.sorted.bam -o YTLV_variant_analysis.vcf
```

22. The .vcf file can be viewed in the terminal using the `less` or `cat` commands. 

---

### D. Assessment of pairwise identity and variant analysis of Phasi Charoen-like virus (PCLV) segment M from contigs 3

23. All of the same steps listed above may be followed to conduct pairwise identity and variant analysis of PCLV segment M (derived from contig 3). 

Note: For these analyses, *PCLV_segM.fasta* was created using the NCBI entry for assession NC_038261.1

---

### Thank you for reviewing this report! Please make sure to review the associate report of these analysis located in: 
