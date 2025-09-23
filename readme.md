## Files in this repository
- ✅ The full list of Human and zoonotic viruses: Complete list of human and zoonotic viruses.md 
- ✅ A brief overview of this workflow: Workflow overview.md
- ✅ A list of 1,111 genomes from human and zoonotic viruses, representing 539 viral species, some of which include multiple distinct strains: complete_precise_human_animal_viral_summary.txt
- ✅ An example out log file: viral_enhanced_comp_40626015.out
## Evaluation of the example results
### There are bacteriophages and bacteria in the results, how to explain it?

=== Kraken2 Classification Report Summary ===
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
---

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
|                        | Enterobacteria  | Bacteria, hosts of phages                        |

---



