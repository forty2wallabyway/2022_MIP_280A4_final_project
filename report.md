## Report of bioinformatic analysis of: Aedes_Recife_R1.fastq
1. Dataset retrieved from thoth01 server at /home/data_for_classes/2022_MIP_280A4/final_project_datasets  
```
cp Aedes_Recife_R1.fastq ~/2022_MIP_280A4_final_project
```

2. Activate conda environment
```
conda activate bio_tools
```

3. Dataset analyzed with fastqc v0.11  
```
fastqc Aedes_Recife_R1.fastq
```

4. Quality report examined

5. Adapter sequences and low quality bases trimmed using cutadapt v3.5  
```
cutadapt \  
  -a AGATCGGAAGAG \  
  -q 30,30 \  
  --minimum-length 80 \  
  -o Aedes_Recife_R1_trimmed2.fastq \  
  Aedes_Recife_R1.fastq \  
  | tee cutadapt2.log
```  
<<cutadapt2.log>>

6. Dataset again analyzed with fastqc v0.11
```
fastqc Aedes_Recife_R1_trimmed2.fastq
```

7. Quality report of trimmed sequences examined

8. Download *Aedes aegypti* genome from RefSeq database into directory on thoth01 server
```
curl -OL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz
```

9. Create index of *A. aegypti* genome to filter out host reads using Bowtie2 v2.4
```
bowtie2-build GCF_002204515.2_AaegL5.0_genomic.fna.gz aedes_aeg_index
```
