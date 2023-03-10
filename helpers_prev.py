# Helper functions for the streamlit-made web application Metabarcoding_16S_V3V4_app.py

import streamlit as st
import qiime2
from qiime2 import Artifact, Metadata

from qiime2.plugins import demux, quality_filter, deblur, metadata, feature_table, vsearch, diversity
from qiime2.plugins.dada2.methods import denoise_single

import os

#os.environ['TMPDIR'] = '/Volumes/HD_ext_CAMI/streamlit_apps/tmp_dir/'
#TMPDIR = os.getenv('TMPDIR')

#@st.cache(show_spinner=True)
def import_paired_end_fastq_gz_files(filepath):
    '''
    Import R1 and R2 fastq.gz files for all samples in the project.
    fastq files must be in the Casava 1.8 format.
    '''
    paired_end_sequences = qiime2.Artifact.import_data(
        'SampleData[PairedEndSequencesWithQuality]', 
        filepath,
        'CasavaOneEightSingleLanePerSampleDirFmt')
    return paired_end_sequences


@st.cache(show_spinner=True)
def vsearch_join_pairs(paired_end_sequences, minovlen, minmergelen, maxmergelen, maxee, maxdiffs):#minovlen=20, minmergelen=460, maxmergelen=480, maxee=1, maxdiffs=5):
    '''
    Performs joining of paired end reads based on given paramenters.
    Default parameters are set for reads' length of 250 bases (v2 Illumina MiSeq kit)
    '''
    demux_joined = vsearch.methods.join_pairs(
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
    return demux_filter


@st.cache(show_spinner=True)
def dada2_denoise_single_joined(demux_filter, N, trim_TF):
    trunc_len = N
    pooling_method = "independent"
    chimera_method = "consensus"
    n_threads = 0 # all available cores
    if trim_TF:
        dada2_sequences = dada2.methods.denoise_single(
            demultiplexed_seqs = demux_filter,
            trunc_len=trunc_len,
            pooling_method=pooling_method, 
            chimera_method=chimera_method,
            n_threads=n_threads)
    else:
        dada2_sequences = dada2.methods.denoise_single(
            demultiplexed_seqs = demux_filter,
            trunc_len=0, # disable trimming
            pooling_method=pooling_method,
            chimera_method=chimera_method,
            n_threads=n_threads)
    return dada2_sequences





@st.cache(show_spinner=True)
def deblur_denoise_trim_paired_end(demux_filter, N, trim_TF):
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
            demux_filter,
            trim_length=N, # disable trimming
            sample_stats=True, 
            jobs_to_start=57)
    else:
        deblur_sequences = deblur.methods.denoise_16S(
            demux_filter,
            trim_length=-1, # disable trimming
            sample_stats=True,
            jobs_to_start=57)
    return deblur_sequences



@st.cache
def import_ref_gg_13_8_otus_seqs(filepath):
    '''
    Import R1 and R2 fastq.gz files for all samples in the project.
    fastq files must be in the Casava 1.8 format.
    '''
    ref_gg_13_8_seqs = qiime2.Artifact.import_data(
        'FeatureData[Sequence]', 
        filepath)
    return ref_gg_13_8_seqs

@st.cache
def import_ref_gg_13_8_otus_taxonomy(filepath):
    '''
    Import R1 and R2 fastq.gz files for all samples in the project.
    fastq files must be in the Casava 1.8 format.
    '''
    ref_gg_13_8_taxonomy = qiime2.Artifact.import_data(
        'FeatureData[Taxonomy]', 
        filepath)
    return ref_gg_13_8_taxonomy

