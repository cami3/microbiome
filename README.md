# Microbiome App - App microbioma

############# Head over to https://microbiomephylo.com for a Comprehensive Microbiome Amplicon Sequencing Downstream Data Analysis, based on phyloseq. MicrobiomePhylo is designed with an intuitive interface to facilitate the analyses and publication-ready plots. ##################


16S rRNA meta-barcoding NGS data analysis app - bacterial microbiome

Web app for user-friendly pre-processing of 16S rRNA NGS (Next Generation Sequencing - Illumina) raw fastq.gz files, running interactively a pipeline based on the QIIME2 (Quantitative Insights Into Microbial Ecology) analysis package. Developed using python programming and the streamlit library. 
The web app performs OTUs/ASVs (Operational Taxonomic Unit/Amplicon Sequence Variant) clustering/denoising and taxonomy assignment and feature table normalization through rarefaction. It provides the raw along with rarefied ASVs counts feature table and the taxonomy table. It also provides the phylogenetic tree-based clustering of samples. Moreover, the web app renders results of alpha and beta diversity analysis as well as taxa relative abundances data, viewable in app and on the qiime2 platform https://view.qiime2.org.
The web app is also provided with functionalities for interactive tertiary data analysis if OTU, taxonomy and samples' associated metadata tables are provided. It produces tables and plots about the taxonomy associated to features, as well as the number of features, as well as the taxa barplots, as well as the alpha and beta diversity metrics comparisons with test statistics, for each user-defined group of samples at each user-defined taxonomic level.

All results tables and visualizations publication-ready are available for download.
Written in italian and english.
Available as a docker image at https://hub.docker.com/r/bac3/microbiome
Upon executing the Docker container, the application becomes accessible on localhost at port 80.


App per l'analisi di dati NGS da esperimenti di meta-barcoding del gene codificante per rRNA 16S - microbioma batterico

Web application per il pre-processamento intuitivo dei dati grezzi da sequenziamento NGS (Next Generation Sequencing - Illumina) del gene codificante per rRNA 16S, nel formato di files fastq.gz. La web app permette lo svolgimento interattivo dell'analisi dei dati grezzi con una pipeline basata sulla piattaforma di analisi QIIME2 (Quantitative Insights Into Microbial Ecology). E' sviluppata in linguaggio di programmazione python con l'utilizzo della libreria "streamlit". 
La web app esegue il clustering/denoising delle sequenze di DNA in OTUs/ASVs (Operational Taxonomic Units/Amplicon Sequence Variants), l'assegnazione della tassonomia alle features e la normalizzazione della tabella delle features attraverso la rarefazione. Fornisce anche il raggruppamento dei campioni basato sull'albero filogenetico. La web app, inoltre, produce i risultati dell'analisi di alfa e beta diversità ed anche dei dati di abbondanza relativa delle features alle varie categorie tassonomiche, visualizzabili in app e sulla piattaforma qiime2 https://view.qiime2.org.
La web app è dotata anche dell'opzione di analisi terziaria interattiva dei dati, se l'utenza abbia già a disposizione le tabelle delle OTU, della tassonomia ed eventualmente dei metadati associati ai campioni. Produce tabelle e grafici della tassonomia associata ad ogni feature, così come del numero delle features ed anche grafici a barre dei taxa e confronti con test statistici tra le metriche di alfa e beta diversità, per ogni raggruppamento dei campioni definito dall'utente, per ogni categoria tassonomica definita dall'utente.

Tutti i risultati in formato di tabelle .csv e di immagini sono scaricabili e pronti per la pubblicazione.
Scritta in italiano e inglese.
Disponibile su Docker Hub all'indirizzo: https://hub.docker.com/r/bac3/microbiome
Dopo l'esecuzione del container Docker, l'applicazione è accessibile su localhost alla porta 80.
