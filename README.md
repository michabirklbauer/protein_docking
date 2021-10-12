# Automatic identification of important interactions and interaction-frequency-based scoring in protein-ligand complexes

## Description

This repository contains code and results for my master's thesis in Data Science and Engineering for the University of Applied Sciences Upper Austria: Hagenberg. Foremost it features the **Protein Interaction Analyzer (PIA)** - a collection of python scripts and workflows to extract interactions and their frequencies from a set of protein-ligand complexes and subsequently score them. This information can be used to predict if a complex is active or not.

PIA largely builds upon [BioPandas](https://github.com/rasbt/biopandas), [RDKit](https://github.com/rdkit/rdkit) and [PLIP](https://github.com/pharmai/plip) and many other python packages. PIA can either be set up with [Anaconda](https://anaconda.org/) using the `environment.yml` file and downloading the PIA scripts or by using the provided Docker image (see [below](#docker-usage)). The PIA scripts are located in `scripts/python/PIA.py` and `scripts/python/PIAScore.py`. However it is recommended to downloaded the newest versions from the [**PIA repository**](https://github.com/michabirklbauer/PIA).

## Repository Structure

- **archive:** Archived files.
- **data:** Data for all 11 analyzed targets including workflows and results.
- **docs:** Thesis, presentations and files that were included in those.
- **scripts:** PIA implementation and some other sample scripts that were used for analyzing the data.
- **tests:** Automated tests for GitHub actions.
- **workflows:** Sample workflows for [scoring](https://github.com/michabirklbauer/protein_docking/tree/master/workflows/scr), [comparing active and inactive frequencies](https://github.com/michabirklbauer/protein_docking/tree/master/workflows/vs) and [extracting interactions and frequencies from a sdf file](https://github.com/michabirklbauer/protein_docking/tree/master/workflows/sdf).

## Protonated Structures

[PLIP](https://github.com/pharmai/plip) will create protonated structures of the target PDB file in a temporary directory when analyzing a complex. The protonated structures may not be deleted automatically after the analysis is complete. In Windows 10 the structures are created in `/Users/*/AppData/Local/Temp` and can be deleted manually.

## Docker Usage

Either manually build the image using the provided Dockerfile:  To manually build the image run the following command:

```bash
docker image build -f Dockerfile -t protein_docking:latest .
```

OR pull it from dockerhub via [michabirklbauer/protein_docking:latest](https://hub.docker.com/r/michabirklbauer/protein_docking):

```bash
docker pull michabirklbauer/protein_docking
```

Once the image is built or downloaded you should create a directory that can be mounted to the container to share files e.g. in Windows 10 create a directory on drive `C:` called `docker_share`. You will have to add this directory in Docker *Settings > Resources > File Sharing* and restart Docker. Then the container can be run with:

```bash
docker run -v C:\docker_share:/exchange -p 8888:8888 -it protein_docking:latest
```

This will mount the directory `C:/docker_share` to the directory `/exchange` in the container and files can be shared freely by copy-pasting into those directories. Python can be run normally in the container and PIA can be imported without needing to install any requirements. To launch [JupyterLab](https://jupyter.org/) in the container you need to run the command:

```bash
jupyter lab --ip=0.0.0.0 --port=8888 --allow-root
```

JupyterLab can then be accessed from any browser via the given link.

## Further Development

PIA will have its standalone repository with more polishing and documentation in [github.com/michabirklbauer/pia](https://github.com/michabirklbauer/PIA). For help or issues with PIA that are not directly connected to the work in this thesis please refer to that repository instead!

## Thesis - Abstract

Molecular docking is an important tool in virtual screening for the discovery and design of new active agents for drug usage. The docking process is influenced by how well molecules fit in the binding site and which interactions occur between the protein and the ligand. Detection of these interactions can be automated with tools like the Protein-Ligand Interaction Profiler (PLIP) by PharmAI. However, identification and assessment of the importance of the different interactions in a protein-ligand complex is still a manual task that requires additional experimental data or domain knowledge about the target. The goals of this thesis are twofold: Firstly, to automatically identify those interactions that have a significant influence on ligand binding, and secondly, to develop a novel scoring function which is able to discriminate active molecules from inactive ones if possible. The underlying data basis were selected targets of the Directory of Useful Decoys: Enhanced (DUD-E) and available structures from the Protein Data Bank (PDB). Specifically 11 targets were analysed: 11-Beta-Hydroxysteroid Dehydrogenase 1 (HSD11B1), Acetylcholinesterase (ACHE), Coagulation Factor XA (FXA), Cyclooxygenase 1 and 2 (COX1/COX2), Dipeptidyl Peptidase IV (DPP4), Monoamine Oxidase B (MAOB), P38 Mitogen-Activated Protein Kinase 14 (MAPK14),  Phosphodiesterase 5 (PDE5A), Protein-Tyrosine Phosphatase 1B (PTP1B) and Soluble Epoxide Hydrolase (SEH). PLIP is used to extract interactions present in a protein-ligand complex and the respective interaction’s frequency is measured across all target structures. Cofactors were excluded from the analysis and hydrophobic interactions were only counted once per residue. Additionally, when analysing docking poses only the pose that had the most interactions contributed to the calculation. Furthermore, four different scoring functions that are based on the differences in frequencies between active and inactive compounds were established and their performance was assessed on an independent test partition containing unseen ligands. The results show that interactions which are known from literature to be important for ligand binding are found for all targets except ACHE, in many cases among the top ranked interactions in terms of frequency. This behaviour implies a relationship between interaction frequency and the interaction’s significance in ligand binding. Interaction-frequency-based scoring was tested in five targets and performed above baseline accuracy in four of the five targets. In all targets scoring led to an enrichment of active compounds and false positive rates fluctuated between 0 and 33%. Interaction frequency analysis and interaction-frequency-based scoring could therefore be used as supporting tools in virtual screening to further enhance results.

[Full Thesis](https://raw.githubusercontent.com/michabirklbauer/protein_docking/master/docs/Thesis_MichaBirklbauer.pdf)

## Contact

- Mail: [micha.birklbauer@gmail.com](mailto:micha.birklbauer@gmail.com)
- Telegram: [https://telegram.me/micha_birklbauer](https://telegram.me/micha_birklbauer)
