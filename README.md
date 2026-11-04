# Childhood Colibactin Exposure

Data and code to reproduce the analyses and figures in:

> Marshall H. et al. *Common early-life exposure of the colorectal epithelium
> to the microbiome-derived mutagen colibactin.* 

## Data availability

Raw sequencing data are deposited on EGA (accession code EGAD00001015830).
Processed mutation calls, metadata and summary tables required to run the code are in `Data/`.

Within the folder `Data/`, 'Main_Outputs/' contain the main summary output data used in this manuscript.
Outputs from re-analysis of *Lee-Six et al., Nature, 2019* are found in 'Data/Lee_Six_Nature_2019/'.
Publically available data from *Robinson et al., Nature Communications, 2019*, *Robinson et al., Nature Genetics, 2021*, and *Cornish et al., Nature, 2024* that were used in analyses are found in their respective folders within 'Data/'.  

## Repository structure
| Folder | Contents |
|---|---|
| `Data/` | Contains summary data required to reproduce analyses in publication.  |
| `Signatures/` | Contains folders detailing extraction of signatures using SigProfiler and msigHDP, and outputs from attribution per crypt using different signature attribution methods. |
| `snv_phylogenetic_trees/` | Phylogenetic trees built from shared SNVs. |
| `Metagenomic_analysis/` | kraken2/metaphlan consensus outputs, with scripts and summary figures used in publication. |
| `Figure1/` | R markdowns associated with Figure 1 and manuscript sections: "Study participants and samples" and "Background mutation rates." |
| `Figure2/` | R markdowns associated with Figure 2 and manuscript section: "Prevalence and mutation burden of colibactin mutagenesis." |
| `Figure3/` | R markdowns associated with Figure 3. |
| `Figure4/` | R markdowns associated with Figure 4 and manuscript sections: "Timing of colibactin mutagenesis." |
| `Supplementary_Figures/` | R markdowns for supplementary figures.  |

## Licence

Code released under the [MIT] licence. See `LICENSE`.

## Contact

Henry Marshall — hm14@sanger.ac.uk
