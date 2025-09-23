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
