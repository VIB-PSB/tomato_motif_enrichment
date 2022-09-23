# Motif enrichment in tomato

The present GitHub repository contains the code and data necessary to reproduce a computational method covering two protocols: one for gene ID mapping between two tomato genome assemblies, and a second protocol for motif enrichment and GRNs inference in tomato. This computational method is covered in the chapter titled "Prediction of transcription factor regulators and gene regulatory networks in tomato using binding site information" from the 2nd edition of "Plant Gene Regulatory Networks" - Methods in Molecular Biology Series.

## **Folder structure**

The folder structure of the repository consists of three folders with the necessary scripts and files to perform the two main protocols (the gene ID conversion protocol and the motif enrichment protocol) of this chapter and the functional network analysis. Additionally, there is a compressed file named “precomputed_outputs.tar.gz” (with the corresponding instructions to uncompress it) that contains the final result files derived from executing the protocols presented in this chapter. Finally, there is an empty output folder, where the results from reproducing this protocol will be stored.

* **gene_id_conversion**: data and scripts to perform gene ID conversion protocol.
  - data: the gff3 file of query genome, gff3 filtered for gene features of target genome, and test set of ITAG 2.5 gene IDs for conversion to ITAG 4.0 gene IDs. Additionally, to make it more visible, in this folder we have added an output file of the protocol: the conversion table file between ITAG 2.5 and ITAG 4.0
  - scripts: script to parse liftoff output and obtain the gene ID mapping table (getGeneIDFile.py), and script to convert set of genes in ITAG 2.5 to ITAG 4.0 (convert_to_ITAG4.py).
* **enrichment_pipeline**: data and scripts to perform motif enrichment.
  - data: motif mapping files, motif information file and set files.
  - scripts: Nextflow pipeline script, configuration file and scripts used by the Nextflow pipeline
* **functional_network_analysis**: scripts and data to perform functional analysis of predicted networks.
  - data: description of tomato genes and arabidopsis orthologs, and jasmonate-related GO terms.
  - scripts: script to get TGs for specific TF of network (get_TGs.py), and to obtain Cytoscape input visualization files per GO category and TF (getVisualizationFiles_perGO.py and getVisualizationFiles_perTF.py).

Requirements:

* [Nextflow version 21.10.6](https://www.nextflow.io/).
* [Python version 3.6.5](https://www.python.org/downloads/release/python-365/).
* [Python package “pandas” version 1.1.5](https://pypi.org/project/pandas/1.1.5/).
* [Liftoff version 1.6.1](https://github.com/agshumate/Liftoff).
* [Bedtools version 2.2.28](https://bedtools.readthedocs.io/en/latest/content/history.html#version-2-28-0-23-mar-2019).
* [Minimap version 2.18](https://github.com/lh3/minimap2/releases/tag/v2.18).
* [Cytoscape version 3.9.1](https://cytoscape.org/download.html).
* [Git version 2.13.1](https://github.com/git-guides/install-git).

## **Contact**

For questions or remarks about the code or data in this repository please contact Nicolás Manosalva Pérez (nicolas.manosalvaperez@psb.vib-ugent.be) or Klaas Vandepoele (klaas.vandepoele@psb.vib-ugent.be)
