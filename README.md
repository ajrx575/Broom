README.md
================

- [Introduction](#introduction)
- [Getting Started](#getting-started)
  - [Requirements](#requirements)
  - [Installation](#installation)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [Notes](#notes)

## Introduction

*Sweeping away extraneous information in your genomic data files,
cleaning them up for downstream processing.*

Broom is a genomic data file pre-processing tool. This package is run
from the command line that will take a GFF file and a genome FASTA file
from Ensembl, RefSeq, and JGI and convert them into files with the
minimum amount of necessary information for important features, such
that files from either repository are structurally indistinguishable.
This tool is currently integrated with the variant annotator snpEff.

Despite being a widely used “standardized” file format, many GFF files
contain non-compliant features that have inappropriate characters in
different columns, poorly defined feature relationships, and excess
information in the attributes column that convolute parsing. It has also
been observed, mainly in files from RefSeq, that in set of files for a
given species, the identifiers used across files do not map to each
other making it increasingly difficult to understand which features and
sequences are supposed of the same family. The RefSeq identifier label
is also messy with weird characters and long descriptions that are also
difficult to parse. These issues can plague files from any repository.
Broom is used to clean up files from RefSeq, Ensmebl, and both of JGI’s
databases to make traceability, parsing, and downstream annotation easy.

## Getting Started

### Requirements

This package is stored in a Docker container, and requires Docker to be
installed to be run. You can install Docker Desktop or the Engine
following the instructions on the Docker website.

### Installation

Broom is being held in a Docker container. You can pull the image using
the following command in the Terminal:
`docker pull ar9481/python_broom:v2` or, more simply, execute a docker
run command as shown in the Usage section below.

Broom has been configured to run on both amd64 and arm64 architectures.

## Usage

This package is run via the command line and takes eight total
arguments, four of which are optional. Broom can be used to solely to
clean up and reconstitute your files for use with other tools, or the
files that are created in the process can be passed along to snpEff for
variant annotation with definition of the optional arguments. Command
line examples for the two use cases are shown below, though quoting of
strings passed as arguments is not necessary and is done here to clarify
the arguments.

``` bash
docker run -v /path/to/data/files/locally:/data \
ar9481/broom:v2 \
--gff /data/'GFF file name' \ 
--genome_fasta /data/'FASTA file name' \ 
--repository 'name of repository of origin' \
--genome_name '(nick)name of species to which the files belong'
```

``` bash
docker run -v /path/to/data/files/locally:/data \
ar9481/broom:v2 \
--gff /data/'GFF file name' \ 
--genome_fasta /data/'FASTA file name' \ 
--repository 'name of repository of origin' \
--genome_name '(nick)name of species to which the files belong'
--vcf /data/path/to/VCF/file \
--config /data/path/to/snpEff/config/file \
--excl 'IDs of chromosomes to not include in GFF' \
-a
```

## Input

1.  **–gff (required)** :The input for this argument must be the path to
    a GFF file from the following repositories: RefSeq, Ensembl, JGI
    (Phytozome and Mycocosm/Phycocosm).

2.  **–genome_fasta (required)**: The input for this argument must be
    the path to a genomic FASTA file from the same repository as the
    GFF. This FASTA file should be from the same annotation version as
    well.

3.  **–repository (required)**: The input for this argument can be any
    of the following options (as a string): Ensembl, RefSeq, JGI-MP
    (Mycocosm/Phycocosm), JGI-Phy (Phytozome).

4.  **–genome_name (required)**: The input for this argument is name of
    the species whose files are being processed. This can be an
    abbreviation/nickname (i.e - PDac for Phoenix dactylifera), but it
    must be the same name that you assign in the snpEff config file if
    snpEff will be used as well.

5.  **–excl**: The input for this argument are the chromosome IDs, as
    strings, for chromosomes that should not be included in the
    resulting GFF if desired. If there are multiple chromosome IDs to be
    excluded, they should be passed as an unquoted, white-space
    separated string after the argument flag

6.  **–vcf**: The input for this argument must be the path to a VCF
    file. It needs to have the same chromosome naming conventions as
    your other input files.

7.  **–config**: The input for this argument must be the path to a
    snpEff appropriate config file, if the intent is to use this tool in
    conjunction with snpEff. Instructions can be found here
    <https://pcingola.github.io/SnpEff/snpeff/build_db/#step-1-configure-a-new-genome>.

8.  **-a**: This argument has no input, as it is meant to act as a flag.
    If the argument is present, all features found in the original GFF
    in the reconstituted GFF. If this flag is not present, then the
    features in the output file will be the ones most commonly
    investigated and compatible with snpEff (gene, mRNA/transcript,
    exon, CDS, UTR, etc.)

## Output

The primary output of Broom consists of five files:

1)  a reconstituted GFF with minimum information that has been resorted
    in the same order as the original file \>
    `{repository}_{species}_interm.gff` ,

2)  a CDS FASTA file where the entry identifiers map back to their mRNA
    parent \> `{repository}_{species}_cds.fa` ,

3)  a protein FASTA file where the entry identifiers map back to their
    mRNA parent \> `{repository}_{species}_protein.fa` ,

4)  a text file containing misformatted features from the GFF \>
    `{repository}_{species}_misformatted_gff_features.txt` ,

5)  a text file containing CDS features that generated translation
    warnings during the creation of the protein FASTA file \>
    `{repository}_{species}_fasta_warnings.txt`.

If Broom is being used in conjunction with snpEff, there are four files
total passed on to the variant annotator: the first three listed above
and the original genomic FASTA file provided. Examples of these outputs
are shown below using RefSeq’s Phoenix dactlyifera data.

<br>

<figure>
<img src="images/gffpic" width="698"
alt="Minimized P. dactlyifera GFF" />
<figcaption aria-hidden="true">Minimized P. dactlyifera GFF</figcaption>
</figure>

<br>

<figure>
<img src="images/protpic" width="559" alt="Minimized protein FASTA" />
<figcaption aria-hidden="true">Minimized protein FASTA</figcaption>
</figure>

<br>

<figure>
<img src="images/cdspic" width="546" alt="Minimized CDS FASTA" />
<figcaption aria-hidden="true">Minimized CDS FASTA</figcaption>
</figure>

<br>

<figure>
<img src="images/gfferrorpic" alt="Misformatted GFF features" />
<figcaption aria-hidden="true">Misformatted GFF features</figcaption>
</figure>

<br>

<figure>
<img src="images/cdserrorpic" alt="Misformatted CDS features" />
<figcaption aria-hidden="true">Misformatted CDS features</figcaption>
</figure>

<br>

In the process of creating these files, there a number of other
intermediate files created including: GFFs for each individual
annotation source (column 2 of GFF), sqlite databases corresponding to
each GFF file, and an unsorted minimized GFF.

## Notes

- The warnings generated during the translation of the CDS sequences are
  almost always due to poor annotations in the original GFF file. This
  issue would occur even if your GFF was unprocessed, and no features
  were deleted from the databases. There is no way to fix this without
  changing the fundamental annotation, which is why features that suffer
  from discordant numbers of bases are written out with warnings, but
  still included in the GFF and FASTA files. It is up to the user to
  check accuracy of variant annotations in these regions.
- There are no substantial differences between v1 and v2. v2 ensures
  that the output files are placed in the same directory that the file
  passed to the GFF argument is located, but v1 could not be replaced at
  the time of publishing due to ongoing scale testing.
