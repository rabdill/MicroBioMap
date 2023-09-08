# Microbiome Compendium


Our dataset includes over 170,000 samples of publicly available 16S rRNA
amplicon sequencing data, all processed using the same pipeline and reference
database.

The goal of the MicroBioMap package is simply to expose these data
to the broad community of Bioconductor and R users with 
the smallest fuss possible. 

See the [vignette](https://seandavi.github.io/MicroBioMap/articles/overview.html) for more details. 

## Installation

For the smoothest installation experience, use the [`BiocManager` Bioconductor 
package](https://bioconductor.org/packages/BiocManager).

```
BiocManager::install('seandavi/MicroBioMap')
```

## Usage

Load the compendium using:

```
library(MicroBioMap)
cpd <- getCompendium()
```

The resulting object is a `TreeSummarizedExperiment` object. Currently, the
"tree" part of the TreeSummarizedExperiment is not populated, but that is
on the roadmap.

After loading the compendium, you will have immediate access to nearly
170,000 microbiome samples. 