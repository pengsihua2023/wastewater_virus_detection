## Files in this repository
- ✅ The full list of Human and zoonotic viruses: Complete list of human and zoonotic viruses.md 
- ✅ A brief overview of this workflow: Workflow overview.md
- ✅ A list of 1,111 genomes from human and zoonotic viruses, representing 539 viral species, some of which include multiple distinct strains: complete_precise_human_animal_viral_summary.txt
- ✅ An example out log file: viral_enhanced_comp_40626015.out
## Evaluation of the example results
The sample for testing our workflow was download from the repository of CDC.
<img width="700" height="140" alt="image" src="https://github.com/user-attachments/assets/0506bad9-b15d-4fe3-9182-8e87bfbd65c2" />


### A summary of the example results
```
Kraken2 Classification Report Summary
87.47%	314 reads	unclassified
12.53%	45 reads	root
12.53%	45 reads	Viruses
8.08%	29 reads	unclassified
6.69%	24 reads	unclassified
6.69%	24 reads	unclassified
6.69%	24 reads	Pandoravirus
1.95%	7 reads	Pandoravirus
1.39%	5 reads	Pandoravirus
1.39%	5 reads	Pandoravirus
0.28%	1 reads	Pandoravirus
1.39%	5 reads	Acanthamoeba
4.46%	16 reads	Duplodnaviria
4.46%	16 reads	Heunggongvirae
4.46%	16 reads	Uroviricota
4.46%	16 reads	Caudoviricetes
3.34%	12 reads	unclassified
3.34%	12 reads	Enterobacteria
0.84%	3 reads	Pantevenvirales
0.84%	3 reads	Straboviridae
```

### Kraken2 Classification Results (Bacteriophages vs. Non-Bacteriophages)

| Category               | Entry           | Notes                                            |
| ---------------------- | --------------- | ------------------------------------------------ |
| **Bacteriophages**     | Duplodnaviria   | Includes dsDNA phages                            |
|                        | Heunggongvirae  | Phage-related viral clade                        |
|                        | Uroviricota     | Phylum of tailed phages                          |
|                        | Caudoviricetes  | Class of tailed phages (e.g., T4, λ)             |
|                        | Pantevenvirales | Viral order containing phages infecting bacteria |
|                        | Straboviridae   | Phage family                                     |
| **Non-Bacteriophages** | Pandoravirus    | Giant viruses infecting amoebae                  |
|                        | Acanthamoeba    | Amoeba host, not a virus                         |
|                        | Enterobacteria  | **Bacteria**, hosts of phages                        |

---
## Why were bacteriophages detected? That's so strange!
I used only human viral genome reference sequences to filter the reads, and the expected result should have been no phage sequences. However, the result was the opposite; most of the aligned sequences were phage sequences. My thoughts are:
1. Kraken2 is unreliable for classifying reads. I plan to replace it with Centrifuge in the workflow in the next step.
2. We had too few aligned reads (I set the threshold to 1000 reads before starting the assembly process), which was insufficient for assembly.
3. The number of human viruses we selected ([539](https://github.com/pengsihua2023/wastewater_viral_detection/blob/main/Complete%20list%20of%20human%20and%20zoonotic%20viruses.md)) is too small. Next, we plan to use [1,131 human viral species](https://github.com/pengsihua2023/wastewater_viral_detection/blob/main/To%20be%20used%20reference%20genome%20list/TableS1.xlsx) (including zoonotic viruses) to filter the reads, based on this article: [Bioinformatics, 2022, DOI: 10.1093/bioinformatics/btac275](https://github.com/pengsihua2023/wastewater_viral_detection/blob/main/To%20be%20used%20reference%20genome%20list/An%20atlas%20of%20human%20viruses.pdf).
4. Our reference genome database includes both DNA and RNA viruses. Therefore, our workflow can detect both DNA and RNA viruses. In this example, we used metagenomic paired-end sequencing data, which allows us to detect only DNA viruses. To detect RNA viruses, we would need metatranscriptomic paired-end sequencing data.
## It is puzzling that 87.47% of the putative human viral reads could not be classified by Kraken2 despite being filtered as viral sequences.
The reason why 87.47% of sequences were unclassified is:  
  
**Filtering retains non-target sequences:** Loose parameters or lack of negative screening allow non-human viruses (e.g., Pandoravirus), human ERVs, or contaminant sequences to pass through the filter. These sequences cannot be matched to the 539 viral taxIDs in Kraken2.  
  
**Sequence variation:** The target virus in the sample may contain variants that differ significantly from the reference sequences in the database, causing k-mer matching to fail.  

**Sequence quality:** Short sequences or low-quality reads are difficult to match with k-mers.  
  
**Kraken2 mechanism:** The k-mer matching and LCA algorithm require high specificity, which is not fully consistent with the global alignment used in filtering.   

Currently, we are using **BWA** as the alignment tool for virus screening; our next step is to try using **Bowtie2** instead.

## I used the reference genomes of 539 human viruses to build a Kraken2 database (to reduce computation time), but why does the Kraken2 classification result show a significant number of reads belonging to Pandoravirus?

- 1 **Classification tree contains Pandoravirus taxID**  
   Kraken2 uses the NCBI Taxonomy classification tree by default, which includes the Pandoravirus taxID. Even though the database only contains taxIDs for 539 human viruses, the classification tree may still introduce Pandoravirus, leading to the LCA algorithm misassigning some reads.

- 2 **BWA screening parameters are too lenient, retaining non-human virus sequences**  
   The default parameters of BWA-MEM (e.g., `-k 19`, `-T 30`) may allow non-target reads that are partially similar to the 539 viruses (e.g., Pandoravirus or Acanthamoeba-related sequences) to pass through the screening process.  

- 3 **Sample complexity and environmental contamination**  
   Metagenomic data may contain environmental viruses (e.g., Pandoravirus, which coexists with Acanthamoeba), and the report detected Acanthamoeba (1.39%, 5 sequences), supporting the possibility of sample contamination.

- 4 **k-mer matching false positives**  
   Pandoravirus may share conserved genes with human viruses (e.g., DNA replication sequences), causing the k-mers of reads retained after BWA screening to match Pandoravirus, triggering misclassification.

## Although Pandoravirus wasn't our target virus, it was detected. Is Pandoravirus truly present in this wastewater environment?
1. **Evidence Supporting Pandoravirus Presence**  
   - **Multiple Read Assignments**: The report lists Pandoravirus at different percentages (6.69%, 1.95%, 1.39%, 0.28%), totaling 42 reads. This consistency suggests that certain reads have been assigned to the Pandoravirus taxID (1349409).  
   - **Co-occurrence with Acanthamoeba**: The detection of Acanthamoeba (1.39%, 5 reads) is significant, as Pandoravirus is a giant virus that often coexists with Acanthamoeba as its host. This ecological relationship supports the possibility of Pandoravirus being present in the sample.  
   - **Environmental Context**: Since the data is metagenomic, likely obtained from a complex environment (e.g., soil, water, or human microbiome), and Pandoravirus inhabits such environments, its presence is reasonable.

2. **Potential Misclassification: Classification May Be Due to the Following Factors**  
   - **Conserved Sequence Similarity**: Pandoravirus may share conserved genes with other viruses (e.g., DNA polymerase), leading to k-mer matches.  
   - **Database/Tree Bias**: The Kraken2 database, built from 539 human viruses, may not fully represent the sample, and the presence of the Pandoravirus taxID in the classification tree could cause misassignments.  
   - **Screening Artifacts**: BWA may have retained non-target reads, which Kraken2 then misclassifies due to the tree structure.

3. **Conclusion**  
   Based solely on this report, it cannot be determined whether Pandoravirus is truly present in the sample. The 42 reads classified as Pandoravirus are more likely a combination of misclassification and potential low-level presence, influenced by the following factors:  
   - The NCBI Taxonomy tree includes the Pandoravirus taxID.  
   - BWA retained non-target reads (e.g., Acanthamoeba-related sequences) due to lenient parameters.  
   - k-mer false positives arise from conserved regions.
## Summary of this workflow
- This workflow functions correctly, but some steps failed to produce results, including viral_assembly, orf_prediction，diamond_analysis，and abundance_estimation, because the number of filtered reads was too small to trigger the sequence assembly process.
- Overall, this workflow is just a preliminary exploration and does not yet meet our requirements for wastewater metagenomic virus detection.
- Further work is needed to continuously optimize and improve the workflow.  
