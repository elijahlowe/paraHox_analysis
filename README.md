paraHox_analysis
==========

This project is to examine the effects of two of the ParaHox genes--xlox and cdx--on the developing gut in the sea urchin _Strongylocentrotus purpuratus_(_Sp_) and _Patiria miniata_(_Pm_). Knowndowns of xlox and cdx were conducted on both _Sp_ and _Pm_ and sequenced using illumina and compared against the wildtype embryos. 

## Read mapping and Differential expression analysis

Reads were first trimmed using Trimmomatic (v0.33) with the scripts **trim_pm.qsub**, **trim_spcdx.qsub** and **trim_splox.qsub**.
_Sp_ reads were mapped to [Genome sequence (V3.1)] and _Pm_ reads were mapped to the [genome sequence (V1.0) Scaffolds] using Tophat (2.0.8b). After mapping reads were sorted using SamTools (v) and counts wre extracted using HTSeq (v0.6.1). The following scripts were used **sp_cdx.qsub**, **sp_lox48.qsub** and **sp_lox72.qsub** for _Sp_ and the the following scripts **pm_cdx.qsub** and **pm_lox.qsub** were used for _Pm_. 

Differentially expressed genes were identified using DESeq2 and after the differentially expressed genes were annotated using extracted information from <echinobase.org> which are in the **data/** directory and using annot_sp.py and annot_pm.py scripts. 
The scripts are used as such, executing from the script folder:
```
python annot_sp.py _input file_ > _output_
```
[Genome sequence (V3.1)]: http://www.echinobase.org/Echinobase/SpDownloads
[genome sequence (V1.0) Scaffolds]: http://www.echinobase.org/Echinobase/PmDownload

## Gene onthology and ortholgy

Gene onthology was conducted with both [Panther] and [BLAST2GO]. Gene orthology was conducted using [OrthoMCL] and the following speices:
* [Human (Homo sapiens)]
* [Mouse (Mus musculus)]
* [Elephant shark (Callorhinchus milii)]
* [Ciona intestinalis]
* [Platypus (Ornithorhynchus anatinus)]
* [Spotted gar (Lepisosteus oculatus)]
* [Xenopus tropicalis]

[Panther]: http://www.pantherdb.org/
[BLAST2GO]: https://www.blast2go.com/
[OrthoMCL]: http://www.orthomcl.org/orthomcl/
[Human (Homo sapiens)]: ftp://ftp.ensembl.org/pub/release-83/fasta/homo_sapiens/pep/
[Mouse (Mus musculus)]: ftp://ftp.ensembl.org/pub/release-83/fasta/mus_musculus/pep/
[Elephant shark (Callorhinchus milii)]: http://esharkgenome.imcb.a-star.edu.sg/download/
[Ciona intestinalis]: ftp://ftp.ensembl.org/pub/release-83/fasta/ciona_intestinalis/pep/
[Platypus (Ornithorhynchus anatinus)]: ftp://ftp.ensembl.org/pub/release-83/fasta/ornithorhynchus_anatinus/pep/
[Spotted gar (Lepisosteus oculatus)]: ftp://ftp.ensembl.org/pub/release-83/fasta/lepisosteus_oculatus/pep/
[Xenopus tropicalis]: ftp://ftp.ensembl.org/pub/release-83/fasta/xenopus_tropicalis/pep/

#Transcription Factors
Transcription Factors ID's were downloaded from the [DBD], and [Metazoa Ensembl biomart] was used to identify the 

[DBD]: http://www.transcriptionfactor.org/index.cgi?Download
[Metazoa Ensembl biomart]: http://metazoa.ensembl.org/biomart/martview/edc6f0a11d97a43a9e640246579050d6
