# Helper functions for the streamlit-made web application Metabarcoding_16S_V3V4_app.py

import streamlit as st
import pandas as pd
import qiime2
from qiime2 import Artifact, Metadata

from qiime2.plugins import demux, quality_filter, deblur, dada2, metadata, feature_table, vsearch, diversity
from qiime2.plugins.dada2.methods import denoise_single

import os
import tempfile

@st.experimental_singleton(show_spinner=True)
def import_paired_end_fastq_gz_files(_filepath):
    '''
    Import R1 and R2 fastq.gz files for all samples in the project.
    fastq files must be in the Casava 1.8 format.
    '''
    paired_end_sequences = qiime2.Artifact.import_data(
        'SampleData[PairedEndSequencesWithQuality]', 
        _filepath,
        'CasavaOneEightSingleLanePerSampleDirFmt')
    return paired_end_sequences


@st.experimental_memo(show_spinner=True)
def import_single_end_fastq_gz_files(_filepath):
    '''
    Import R1 fastq.gz files for all samples in the project.
    fastq files must be in the Casava 1.8 format.
    '''
    single_end_sequences = qiime2.Artifact.import_data(
        'SampleData[SequencesWithQuality]', 
        _filepath,
        'CasavaOneEightSingleLanePerSampleDirFmt')
    return single_end_sequences

@st.cache(show_spinner=True)
def vsearch_join_pairs(paired_end_sequences, minovlen, minmergelen, maxmergelen, maxee, maxdiffs):#minovlen=20, minmergelen=460, maxmergelen=480, maxee=1, maxdiffs=5):
    '''
    Performs joining of paired end reads based on given paramenters.
    Default parameters are set for reads' length of 250 bases (v2 Illumina MiSeq kit)
    '''
    demux_joined = vsearch.methods.merge_pairs(
        demultiplexed_seqs = paired_end_sequences,
        minovlen=minovlen,
        minmergelen=minmergelen,
        maxmergelen=maxmergelen,
        maxee=maxee,
        maxdiffs=maxdiffs)
    return demux_joined



@st.cache(show_spinner=True)
def quality_filter_paired_end(demux_joined, min_quality, quality_window):
    '''
    Performs quality filtering of paired-end fastq reads based on phred64 Illumina quality score.
    '''
    demux_filter = quality_filter.methods.q_score(demux=demux_joined, min_quality=min_quality, quality_window=quality_window)
    secure_temp_dir_q_filter_summary = tempfile.mkdtemp(prefix="temp_", suffix="_q_filter_summary")
    filter_stats = metadata.visualizers.tabulate(demux_filter.filter_stats.view(qiime2.Metadata))
    filter_stats.visualization.export_data(secure_temp_dir_q_filter_summary)

    df_q_filter = pd.read_csv(secure_temp_dir_q_filter_summary+'/metadata.tsv', sep='\t' , skiprows=[1], index_col='sample-id')
    df_q_filter.columns = ["total_input_reads", 
        "total_retained_reads", 
        "reads_truncated", 
        "reads_too_short_after_truncation",
        "reads_exceeding_maximum_ambiguous_bases"]
    df_q_filter = df_q_filter.eval("""
        perc_retained_reads = ((total_retained_reads * 100)/total_input_reads)
        """)
    cols_renameing = {
        "total_input_reads": 'totale_sequenze_iniziali', 
        "total_retained_reads": 'totale_sequenze_accettabili', 
        "reads_truncated": 'sequenze_troncate', 
        "reads_too_short_after_truncation": 'sequenze_troncate_troppo_corte_scartate',
        "reads_exceeding_maximum_ambiguous_bases": ('sequenze_con_oltre_%s_basi_ambigue_scartate' %(quality_window)),
        'perc_retained_reads': 'percentuale_sequenze_accettabili'
        }
    
    df_q_filter = df_q_filter.rename(cols_renameing, axis=1)
    return demux_filter, df_q_filter, secure_temp_dir_q_filter_summary

@st.cache_resource
def app_alpha_rare_curves(_table, max_depth, metrics):
    alpha_rare_curves = alpha_rarefaction(table = _table, max_depth=max_depth, metrics=metrics)
    return alpha_rare_curves


@st.cache_resource
def app_align_to_tree_mafft_fasttree(_sequences, _table, sampling_depth, _metadata):
    result_alignment = align_to_tree_mafft_fasttree(sequences = _sequences)
    seqs_alignment = result_alignment.alignment
    masked_alignment = result_alignment.masked_alignment
    tree = result_alignment.tree
    rooted_tree = result_alignment.rooted_tree

    secure_temp_dir_phylogenetic_tree = tempfile.mkdtemp(prefix="temp_", suffix="_phylogenetic_tree")
    rooted_tree.export_data(secure_temp_dir_phylogenetic_tree)
    
    result_core_metrics = core_metrics_phylogenetic(table=_table, phylogeny=rooted_tree, sampling_depth=sampling_depth, metadata=_metadata)
    core_metr_phylo = result_core_metrics
    return seqs_alignment, masked_alignment, tree, rooted_tree, core_metr_phylo, secure_temp_dir_phylogenetic_tree


@st.cache(show_spinner=True)
def dada2_denoise_single_joined(_demux_filter, N, trim_TF):
    trunc_len = N
    pooling_method = "independent"
    chimera_method = "consensus"
    n_threads = 0 # all available cores
    if trim_TF:
        dada2_sequences = dada2.methods.denoise_single(
            demultiplexed_seqs = _demux_filter,
            trunc_len=trunc_len,
            pooling_method=pooling_method, 
            chimera_method=chimera_method,
            n_threads=n_threads)
    else:
        dada2_sequences = dada2.methods.denoise_single(
            demultiplexed_seqs = _demux_filter,
            trunc_len=0, # disable trimming
            pooling_method=pooling_method,
            chimera_method=chimera_method,
            n_threads=n_threads)
    return dada2_sequences





@st.cache(show_spinner=True)
def deblur_denoise_trim_paired_end(_demux_filter, N, trim_TF):
    '''
    Performs denoising of data and trimming at position N, which is defined by the user based on
    the filter stats --> position in the reads where the median quality score drops too low.
    Only forward reads are supported at this time, so perform paired-end reads joining first.
    If trim_TF is True, trimming is enabled at length N.
    Otherwise, if trim_TF is False, disable trimming.
    "Deblur operates only on same length sequences. " from the web https://forum.qiime2.org/t/deblur-plugin-error/2172
    '''
    if trim_TF:
        deblur_sequences = deblur.methods.denoise_16S(
            _demux_filter,
            trim_length=N, # disable trimming
            sample_stats=True, 
            jobs_to_start=57)
    else:
        deblur_sequences = deblur.methods.denoise_16S(
            _demux_filter,
            trim_length=-1, # disable trimming
            sample_stats=True,
            jobs_to_start=57)
    return deblur_sequences


@st.experimental_memo(show_spinner=True)
def import_SequencesWithQuality(_filepath):
    '''
    Import R1 fastq.gz files for all samples in the project.
    fastq files must be in the Casava 1.8 format.
    '''
    single_end_sequences = qiime2.Artifact.import_data(
        'SampleData[SequencesWithQuality]', 
        _filepath,
        'SingleLanePerSampleSingleEndFastqDirFmt')
    return single_end_sequences


@st.cache
def import_gg_13_8_pre_trained_classifier(filepath):
    '''
    Import reference Green Genes 13/8 otus sequences
    '''
    ref_gg_13_8_pre_trained_classifier = qiime2.Artifact.import_data(
        'TaxonomicClassifier', 
        filepath)
    return ref_gg_13_8_pre_trained_classifier

@st.cache
def import_ref_gg_13_8_otus_seqs(filepath):
    '''
    Import reference Green Genes 13/8 otus sequences
    '''
    ref_gg_13_8_seqs = qiime2.Artifact.import_data(
        'FeatureData[Sequence]', 
        filepath)
    return ref_gg_13_8_seqs

@st.cache
def import_ref_gg_13_8_otus_taxonomy(filepath):
    '''
    Import reference Green Genes 13/8 taxonomy
    '''
    ref_gg_13_8_taxonomy = qiime2.Artifact.import_data(
        'FeatureData[Taxonomy]', 
        filepath)
    return ref_gg_13_8_taxonomy

