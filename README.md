# MicrobiomePhylo App

Explore the depths of microbiome research with MicrobiomePhylo, your premier destination for comprehensive microbiome amplicon sequencing downstream data analysis. Powered by the robust capabilities of the phyloseq and vegan packages, our web application stands as a beacon for researchers seeking to navigate the complexities of microbiome data with ease and precision.

### MicrobiomePhylo: A Gateway to Advanced Microbiome Analysis

Our platform, accessible at https://microbiomephylo.com, is meticulously designed to offer an unparalleled user experience. With an intuitive interface at its core, MicrobiomePhylo simplifies the journey from raw data to publication-ready visualizations, ensuring that your research not only meets but exceeds the standards of scientific excellence.

### Seamless Integration with QIIME2 Artifacts

Understanding the importance of seamless data integration, MicrobiomePhylo offers full compatibility with QIIME2 artifacts. This feature empowers researchers to effortlessly analyze their OTU and taxonomy tables, rooted trees, and metadata, all in the .qza format. By bridging the gap between data collection and analysis, our platform ensures that your focus remains on uncovering the insights hidden within your microbiome data.

# Key Features:

### Comprehensive Analysis Tools: 
Leverage the power of phyloseq for in-depth analysis of microbiome amplicon sequencing data.
### Intuitive User Interface: 
Navigate through the analysis process with ease, thanks to our user-friendly platform design.
### Seamless QIIME2 Integration: 
Import your .qza files directly into MicrobiomePhylo for a streamlined analysis experience.
### Publication-Ready Visualizations: 
Generate visually compelling plots that are ready to enhance the impact of your research publications.

Join the Community of Microbiome Researchers

Whether you're delving into the microbiome for the first time or you're an experienced researcher looking for a more efficient analysis solution, MicrobiomePhylo is here to support your scientific endeavors. Visit us at https://microbiomephylo.com and take the first step towards unlocking the full potential of your microbiome data.

Embark on your journey with MicrobiomePhylo and transform the way you analyze microbiome amplicon sequencing data. Together, let's advance the frontiers of microbiome research.




# Microbiome App - App microbiome
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
