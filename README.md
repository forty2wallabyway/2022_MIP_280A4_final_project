# Final Project for MIP 280A4 (Fall 2022)

### General Overview
This project aims to characterize metagenomic shotgun data from sequenced total RNA extracted from a pool of *Aedes aegypti* mosquitoes. Intially captured in Recifce, Brazil, these mosquitoes now represent an established colony within CSU's Center for Vectorborne Infectious Disease (CVID) annex. These mosquitoes were not previously experimentally infected or treated with any compounds. 

To accomplish the aim of this project, sequence reads will be quality assesed, trimmed, host-filtered, and assembled. Sequence reads will then be mapped against indices created from the assembled contigs to derive the **top 12** contigs with the greatest coverage depth. These contigs will be referenced against the NCBI database using the BLAST tool to identify sequence similarities. For assembled contigs of particular interest, further investigation will include genome annotation, pairwise alginment with exisiting NCBI sequences, phylogenetic analysis, and variant characterization.

### Data Source & Work Environment
- Sequence reads analyzed here were produced by an Illumina NextSeq500 using single-end, 150 basepair reads.  

- Bioinformatic analysis was conducted on CSU's department of MIP thoth01 server. Software installed in the conda environment used for this project can be found in [conda_list.log](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10168390/conda_list.log) 

- Specific tools (including versions) used for this analysis can be found in [bio_tools.yaml](https://raw.githubusercontent.com/stenglein-lab/MIP_280A4_Fall_2022/main/conda_environment/bio_tools.yaml)

### Analysis & Outputs
- A detailed workflow (including shell commands) of the bioinformatic analysis conducted for this project can be found in this repo's [report.md](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/blob/289c6ae8d4768d59acb9e8be637140568838ed4e/report.md) 

- A PDF report of the analysis, including biological interpretation, can be reviewed in [Final_Presentation_Sherman.pdf](https://github.com/forty2wallabyway/2022_MIP_280A4_final_project/files/10186521/Final_Presentation_Sherman.pdf)
