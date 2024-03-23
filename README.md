# Modstrings: working with modified nucleotide sequences <img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/Modstrings/Modstrings.png" height="200" align="right">

<!-- badges: start -->
[![R-CMD-check-bioc-devel](https://github.com/FelixErnst/Modstrings/workflows/R-CMD-check-bioc-devel/badge.svg)](https://github.com/FelixErnst/Modstrings/actions/)
[![BioC Build](https://bioconductor.org/shields/build/devel/bioc/Modstrings.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/Modstrings/)
[![codecov](https://codecov.io/gh/FelixErnst/Modstrings/branch/master/graph/badge.svg)](https://codecov.io/gh/FelixErnst/Modstrings)
[![BioC Years](https://bioconductor.org/shields/years-in-bioc/Modstrings.svg)](https://doi.org/doi:10.18129/B9.bioc.Modstrings)
<!-- badges: end -->

RNA are usually in some form post-transcriptionally modified. Most prominent
examples are of course rRNA and tRNA, but in recent years mRNA was also
discovered to be post-transcriptionally modified.

In many resources, like the tRNAdb ([Juehling et al. 2009](#Literature)) or the
modomics database ([Boccaletto et al. 2018](#Literature)), a dictionary for
modified nucleotides was published. However, in the Bioconductor universe these
information were not directly accessible or representable, since they rely 
extensively on special characters in the RNA modification alphabet.

Therefore, the`Modstrings` package implements the `ModRNAString` class by
extending the `BString` class from the `Biostrings` ([Pages et
al.](#Literature)) package. It can store RNA sequences containing special
characters of the RNA modification alphabet and thus can store location and 
identity of modifications. Functions for conversion to a tabular format are 
implemented as well. A `ModDNAString` class is implemented analogously, which
is based on the DNA modification alphabet from the Hoffman lab ([Sood et
al.](#Literature)).

The implemented classes inherit most of the functions from the parental 
`BString` class and it derivatives, which allows them to behave like the 
normal `XString` classes within the bioconductor universe.

# Installation

The current version of the `Modstrings` package is available from 
Bioconductor.
 
```{r}
# Installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Modstrings")
# Load and attach the package
library("Modstrings")
```

# Literature

- Jühling, Frank; Mörl, Mario; Hartmann, Roland K.; Sprinzl, Mathias; Stadler,
Peter F.; Pütz, Joern (2009): "TRNAdb 2009: Compilation of tRNA Sequences and
tRNA Genes." Nucleic Acids Research 37 (suppl_1): D159–D162.
doi:[10.1093/nar/gkn772](https://doi.org/10.1093/nar/gkn772). 
- Boccaletto, Pietro; Machnicka, Magdalena A.; Purta, Elzbieta; Piatkowski,
Pawel; Baginski, Blazej; Wirecki, Tomasz K. et al. (2018): "MODOMICS: a database
of RNA modification pathways. 2017 update." Nucleic Acids Res. 46 (D1),
D303-D307. doi:[10.1093/nar/gkx1030](https://doi.org/10.1093/nar/gkx1030).
- Pagès, H.; Aboyoun, P.; Gentleman, R.; DebRoy, S. (2018). "Biostrings: 
Efficient manipulation of biological strings." R package version 2.50.1.
- Sood, Ankur Jai, Coby Viner, and Michael M. Hoffman. 2019. “DNAmod: The Dna 
Modification Database.” Journal of Cheminformatics 11 (1):30. 
doi:[10.1186/s13321-019-0349-4](https://doi.org/10.1186/s13321-019-0349-4).
