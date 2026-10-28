# Childhood Colibactin Exposure

Data and code to reproduce the analyses and figures in:

> Marshall H. et al. *Common early-life exposure of the colorectal epithelium
> to the microbiome-derived mutagen colibactin.* 

## Data availability

Raw sequencing data are deposited on EGA (accession code EGAD00001015830)
Processed mutation calls, metadata and summary tables required to run the code are in `Data/`.

## Repository structure
| Folder | Contents |
|---|---|
| `Data/` | Contains summary data required to reproduce analyses in publication  |
| `Signatures/` | Contains folders detailing extraction of signatures using SigProfiler and msigHDP, and outputs from attribution per crypt using different signature attribution methods |
| `snv_phylogenetic_trees/` | Phylogenetic trees built from shared SNVs |
| `Metagenomic_analysis/` | kraken2/metaphlan consensus outputs, with scripts and summary figures used in publications |
| `Figure1/` – `Figure4/` | R markdowns generating each main figure with associated figure outputs |
| `Supplementary_Figures/` | R markdowns for supplementary figures with associated figure outputs |

## Reproducing the analyses

Each `FigureX/` folder contains self-contained R markdowns that read from
`Data/` using paths relative to the repository root.

## Licence

Code released under the [MIT] licence. See `LICENSE`.

## Contact

Henry Marshall — hm14@sanger.ac.uk
