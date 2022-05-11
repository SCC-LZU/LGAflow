# LGAflow

LGAflow - A genome annotation pipeline for eukaryotes

LGAflow stands for Labtorary Genome Annotation workflow which is based on Nextflow framework.

> NOTE: LGAflow is still under development, there might be some bugs during running

- [LGAflow](#lgaflow)
  - [Quick Start](#quick-start)
    - [Install Nextflow](#install-nextflow)
    - [Install Miniconda](#install-miniconda)
    - [Run Test](#run-test)
  - [Usage](#usage)
    - [profile](#profile)
    - [inputs](#inputs)
    - [databases](#databases)
    - [computing options](#computing-options)
  - [Citations](#citations)

## Quick Start

LGAflow is written in [`Nextflow`](https://github.com/nextflow-io/nextflow?msclkid=8ee65eeaa9df11ec9e177e3bb2743c73) and can be used in any Linux distribution. LGAflow uses `conda` or Linux Container (eg. `Docker`, `Singularity`) to manage dependencies. `Nextflow` and either `conda`, `Docker`, `Singuilarity` need to be installed before running the pipeline.

### Install [Nextflow](https://nf-co.re/usage/installation)

1. Typical installation
    ``` bash
    # Make sure that Java v8+ is installed:
    java -version

    # Install Nextflow
    curl -fsSL get.nextflow.io | bash

    # Add Nextflow binary to your user's PATH:
    # mv nextflow ~/bin/

    # OR system-wide installation:
    # sudo mv nextflow /usr/local/bin
    ```
2. Using Bioconda
    ``` bash
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    conda install nextflow
    ```

3. Using pre-built binaiers
    ``` bash
    # Make sure that Java v8+ is installed:
    java -version

    # Get the executable from https://github.com/nextflow-io/nextflow/releases eg :
    wget https://github.com/nextflow-io/nextflow/releases/download/v21.10.6/nextflow-21.10.6-all -O nextflow

    # Add Nextflow binary to your user's PATH:
    # mv nextflow ~/bin/

    # OR system-wide installation:
    # sudo mv nextflow /usr/local/bin
    ```

### Install [Miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

### Run Test
Test data can be obtained from http://seafile.itslab.lzu.edu.cn/f/f886da69fd4b4aafad65/

## Usage
eg:
```
nextflow run dakehero/lgaflow -profile local,conda --cores 8 --memory 16.GB --genome ./test/c_elegans.fa --proteins ./test/proteins/ --reads ./test/RNA_Seq/ --species c_elegans --rm_species Caenorhabditis --augustus_species caenorhabditis --busco_db ./database/nematoda_odb10/ 
```

The above command shows how to perform genome annotation pipeline for *Caenorhabditis elegans* with LGAflow, several options should be specified.

### profile
- data
  - `test`: using test data to run test annotation.
- executer(choose one)
  - `local`: use local executer to run the pipeline on a single-node computer.
  - `slurm`: use slurm executer to run the pipeline on a HPC cluster.
  - `pbs`: use torque executer to run the pipeline on a HPC cluster.
- environment engine(choose one)
  - `conda`
  - `docker`
  - `singularity`
### inputs
- `genome`: genome assemably file in fasta format(eg: species.fa)
- `reads`: path of RNA-Seq fastq files (eg: /path/of/reads/)
- `proteins`: path of proteins file in fasta format(eg: /path/to/proteins)
### databases
- `busco_db`: path of OrthoDB for the specie
### computing options
- `cores`: max CPU cores per process
- `memory`: max memory per process
- `queueSize`max queue size for batch scheduler
- `useMamba` using `mamba` instead of `conda` to create conda-env
## Citations
- [1]“The UniVec Database.” https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/ (accessed Dec. 16, 2021).
- [2]M. Tarailo-Graovac and N. Chen, “Using RepeatMasker to identify repetitive elements in genomic sequences,” Curr Protoc Bioinformatics, vol. Chapter 4, p. Unit 4.10, Mar. 2009, doi: 10.1002/0471250953.bi0410s25.
- [3]J. Storer, R. Hubley, J. Rosen, T. J. Wheeler, and A. F. Smit, “The Dfam community resource of transposable element families, sequence models, and genome annotations,” Mobile DNA, vol. 12, no. 1, p. 2, Jan. 2021, doi: 10.1186/s13100-020-00230-y.
- [4]M. Stanke and S. Waack, “Gene prediction with a hidden Markov model and a new intron submodel,” Bioinformatics, vol. 19 Suppl 2, pp. ii215-225, Oct. 2003, doi: 10.1093/bioinformatics/btg1080.
- [5]G. S. C. Slater and E. Birney, “Automated generation of heuristics for biological sequence comparison,” BMC Bioinformatics, vol. 6, p. 31, Feb. 2005, doi: 10.1186/1471-2105-6-31.
- [6]F. A. Simão, R. M. Waterhouse, P. Ioannidis, E. V. Kriventseva, and E. M. Zdobnov, “BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs,” Bioinformatics, vol. 31, no. 19, pp. 3210–3212, Oct. 2015, doi: 10.1093/bioinformatics/btv351.
- [7]W. Shen, S. Le, Y. Li, and F. Hu, “SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation,” PLOS ONE, vol. 11, no. 10, p. e0163962, Oct. 2016, doi: 10.1371/journal.pone.0163962.
- [8]A. Sharma, R. Podolsky, J. Zhao, and R. A. McIndoe, “A modified hyperplane clustering algorithm allows for efficient and accurate clustering of extremely large datasets,” Bioinformatics, vol. 25, no. 9, pp. 1152–1157, May 2009, doi: 10.1093/bioinformatics/btp123.
- [9]C. L. Schoch et al., “NCBI Taxonomy: a comprehensive update on curation, resources and tools,” Database, vol. 2020, p. baaa062, Jan. 2020, doi: 10.1093/database/baaa062.
- [10]A. R. Quinlan and I. M. Hall, “BEDTools: a flexible suite of utilities for comparing genomic features,” Bioinformatics, vol. 26, no. 6, pp. 841–842, Mar. 2010, doi: 10.1093/bioinformatics/btq033.
- [11]W. H. Majoros, M. Pertea, and S. L. Salzberg, “TigrScan and GlimmerHMM: two open source ab initio eukaryotic gene-finders,” Bioinformatics, vol. 20, no. 16, pp. 2878–2879, Nov. 2004, doi: 10.1093/bioinformatics/bth315.
- [12]M. Lataretu and M. Hölzer, “RNAflow: An Effective and Simple RNA-Seq Differential Gene Expression Pipeline Using Nextflow,” Genes, vol. 11, no. 12, p. 1487, Dec. 2020, doi: 10.3390/genes11121487.
- [13]D. Kim, J. M. Paggi, C. Park, C. Bennett, and S. L. Salzberg, “Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype,” Nat Biotechnol, vol. 37, no. 8, pp. 907–915, Aug. 2019, doi: 10.1038/s41587-019-0201-4.
- [14]K. J. Hoff and M. Stanke, “Predicting Genes in Single Genomes with AUGUSTUS,” Current Protocols in Bioinformatics, vol. 65, no. 1, p. e57, 2019, doi: 10.1002/cpbi.57.
- [15]B. J. Haas et al., “Automated eukaryotic gene structure annotation using EVidenceModeler and the Program to Assemble Spliced Alignments,” Genome Biology, vol. 9, no. 1, p. R7, Jan. 2008, doi: 10.1186/gb-2008-9-1-r7.
- [16]M. G. Grabherr et al., “Trinity: reconstructing a full-length transcriptome without a genome from RNA-Seq data,” Nat Biotechnol, vol. 29, no. 7, pp. 644–652, May 2011, doi: 10.1038/nbt.1883.
- [17]S. Chen, Y. Zhou, Y. Chen, and J. Gu, “fastp: an ultra-fast all-in-one FASTQ preprocessor,” Bioinformatics, vol. 34, no. 17, pp. i884–i890, Sep. 2018, doi: 10.1093/bioinformatics/bty560.
- [18]B. Buchfink, C. Xie, and D. H. Huson, “Fast and sensitive protein alignment using DIAMOND,” Nat Methods, vol. 12, no. 1, Art. no. 1, Jan. 2015, doi: 10.1038/nmeth.3176.
