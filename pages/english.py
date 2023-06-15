import streamlit as st
import streamlit.components.v1 as components
import streamlit_ext as ste # ste.download_button() permette di non eseguire il reload della pagina ad ogni interazione (st.download_button() innesca invece il reload)
from streamlit.runtime.scriptrunner import get_script_run_ctx

import pandas as pd
from matplotlib import pyplot as plt

import plotly.express as px
import numpy as np
import functools
import altair as alt
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.colors import n_colors
from itertools import cycle

import os, shutil
from os.path import basename
import zipfile
from zipfile import ZipFile

from contextlib import redirect_stderr
import sys
import io
from io import StringIO

from bs4 import BeautifulSoup
import pathlib

import base64

import tempfile
from tempfile import NamedTemporaryFile

from flask import Flask
from flask_sslify import SSLify

import qiime2
from qiime2 import Artifact, Metadata

from qiime2.plugins import demux, quality_filter, deblur, dada2, metadata, feature_table, vsearch, diversity

from qiime2.plugins.dada2.methods import denoise_single

from qiime2.plugins.phylogeny.pipelines import align_to_tree_mafft_fasttree
from qiime2.plugins.diversity.pipelines import core_metrics_phylogenetic
from qiime2.plugins.diversity.visualizers import alpha_rarefaction
from qiime2.plugins.diversity.pipelines import beta
from qiime2.plugins.diversity.visualizers import beta_group_significance
from qiime2.plugins.emperor.visualizers import plot
from qiime2.plugins.diversity.methods import pcoa


from qiime2.plugins.feature_classifier.pipelines import classify_hybrid_vsearch_sklearn

from qiime2.plugins.taxa.visualizers import barplot


# from helpers import filter_dataframe, import_paired_end_fastq_gz_files, import_single_end_fastq_gz_files, \
# quality_filter_paired_end, \
# vsearch_join_pairs, import_SequencesWithQuality, \
# deblur_denoise_trim_paired_end, app_alpha_rare_curves, import_ref_gg_13_8_otus_seqs, \
# import_ref_gg_13_8_otus_taxonomy, dada2_denoise_single_joined, \
# 	app_align_to_tree_mafft_fasttree

import locale
#todo:
# ANALISI SECONDARIA:
# FATTO, MANCANO DA INSERIRE I VALORI NUMERICI: aggiorni i suggerimenti dei parametri di fusione delle reads del frammento 
# full length (in base al tipo di libreria e alla lunghezza delle reads)
# aggiorni l'expander sulla fusione delle reads in real-time: catch standard error del merging di vsearch
# sistemi BLAST per le sequenze ASVs rappresentative
# ANALISI TERZIARIA:
# alfa rarefazione: read html statistiche kruskal wallis tutti i gruppi 
# beta rarefazione

# LE FUNZIONI CUSTOM SONO NEL FILE helpers.py


# GIRA NELL AMBIENTE CONDA: qiime2-2023.2, streamlit version 1.11.1



#locale.setlocale(locale.LC_ALL, 'it_IT')



st.set_page_config(
	page_title="App microbiome metabarcoding data analysis",
	page_icon="🧊",
	layout="wide",
	initial_sidebar_state="expanded",
	menu_items={
		'Get help': 'mailto:cami3taf@gmail.com',
		'Report a bug': 'mailto:cami3taf@gmail.com',
		'About': "# *App microbiome metabarcoding data analysis* Info/contact: cami3taf@gmail.com"
	}
)

# set plotly images config, for high resolution png download
config = {
  'toImageButtonOptions': {
	'format': 'png', # one of png, svg, jpeg, webp
	'filename': 'custom_image',
	'height': 500,
	'width': 700,
	'scale':6 # Multiply title/legend/axis/canvas sizes by this factor
  }
}

# # Hide hamburger menu
# st.markdown(""" <style>
# #MainMenu {visibility: hidden;}
# footer {visibility: hidden;}
# </style> """, unsafe_allow_html=True)


# Condense layout removing the padding between components, to reduce the amount of scrolling
# the users need to do
padding = 2
st.markdown(f""" <style>
	.reportview-container .main .block-container{{
		padding-top: {padding}rem;
		padding-right: {padding}rem;
		padding-left: {padding}rem;
		padding-bottom: {padding}rem;
	}} </style> """, unsafe_allow_html=True)

st.set_option('deprecation.showPyplotGlobalUse', False)


app = Flask(__name__)
sslify = SSLify(app)

# Titolo
st.title('Microbiome App')


st.markdown('''The Microbiome App is a powerful tool for interactive analysis of NGS \
	(Next-Generation Sequencing) data from 16S rRNA gene meta-barcoding experiments. \
		It offers an intuitive interface and utilizes the Qiime2 analysis package to \
			process and interpret microbiome data. With integrated result visualization, \
				the app enables researchers to gain insights into bacterial community composition \
					and explore the complexities of microbiome research.''')

# Descrizione
descr_menu_plchldr = st.empty()
descr_menu_plchldr.info('The sidebar menu displays available options and analysis steps.')


def inject_gad():
	"""Add this in your streamlit app.py
	see https://github.com/streamlit/streamlit/issues/969
	"""
	# new tag method
	GA_ID = "google_adsense"
	# NOTE: you should add id="google_analytics" value in the GA script
	# https://developers.google.com/analytics/devguides/collection/analyticsjs
	GA_JS = """
	<script async src="https://pagead2.googlesyndication.com/pagead/js/adsbygoogle.js?client=ca-pub-################"
	 crossorigin="anonymous"></script>
	"""
	# Insert the script in the head tag of the static template inside your virtual
	index_path = pathlib.Path(st.__file__).parent / "static" / "index.html"
	#logging.info(f'editing {index_path}')
	soup = BeautifulSoup(index_path.read_text(), features="lxml")
	if not soup.find(id=GA_ID):  # if cannot find tag
		bck_index = index_path.with_suffix('.bck')
		if bck_index.exists():
			shutil.copy(bck_index, index_path)  # recover from backup
		else:
			shutil.copy(index_path, bck_index)  # keep a backup
		html = str(soup)
		new_html = html.replace('<head>', '<head>\n' + GA_JS)
		index_path.write_text(new_html)
	
def inject_gad1():
	"""Add this in your streamlit app.py
	see https://github.com/streamlit/streamlit/issues/969
	"""
	# new tag method
	GA_ID = "google_adsense"
	# NOTE: you should add id="google_analytics" value in the GA script
	# https://developers.google.com/analytics/devguides/collection/analyticsjs
	GA_JS = """
	<script async src="https://pagead2.googlesyndication.com/pagead/js/adsbygoogle.js?client=ca-pub-#################"
		crossorigin="anonymous"></script>
	<!-- pubblicità prova -->
	<ins class="adsbygoogle"
		style="display:block"
		data-ad-client="ca-pub-################"
		data-ad-slot="3333170307"
		data-ad-format="auto"
		data-full-width-responsive="true"></ins>
	<script>
		(adsbygoogle = window.adsbygoogle || []).push({});
	</script>
	"""
	# Insert the script in the head tag of the static template inside your virtual
	index_path = pathlib.Path(st.__file__).parent / "static" / "index.html"
	#logging.info(f'editing {index_path}')
	soup = BeautifulSoup(index_path.read_text(), features="lxml")
	if not soup.find(id=GA_ID):  # if cannot find tag
		bck_index = index_path.with_suffix('.bck')
		if bck_index.exists():
			shutil.copy(bck_index, index_path)  # recover from backup
		else:
			shutil.copy(index_path, bck_index)  # keep a backup
		html = str(soup)
		new_html = html.replace('<head>', '<head>\n' + GA_JS)
		index_path.write_text(new_html)
	

def inject_ga():
	GA_ID = "google_analytics"

	# Note: Please replace the id from G-XXXXXXXXXX to whatever your
	# web application's id is. You will find this in your Google Analytics account
	
	GA_JS = """
	<!-- Google tag (gtag.js) -->
	<script async src="https://www.googletagmanager.com/gtag/js?id=G-##########"></script>
	<script>
		window.dataLayer = window.dataLayer || [];
		function gtag(){dataLayer.push(arguments);}
		gtag('js', new Date());

		gtag('config', 'G-##########', { 'anonymize_ip': true });
	</script>
	<!-- Google tag (gtag.js) -->
	<script async src="https://www.googletagmanager.com/gtag/js?id=G-##########"></script>
	<script>
		window.dataLayer = window.dataLayer || [];
		function gtag(){dataLayer.push(arguments);}
		gtag('js', new Date());

		gtag('config', 'G-##########');
	</script>
	"""

	# Insert the script in the head tag of the static template inside your virtual
	index_path = pathlib.Path(st.__file__).parent / "static" / "index.html"
	# logging.info(f'editing {index_path}')
	soup = BeautifulSoup(index_path.read_text(), features="html.parser")
	if not soup.find(id=GA_ID):  # if cannot find tag
		bck_index = index_path.with_suffix('.bck')
		if bck_index.exists():
			shutil.copy(bck_index, index_path)  # recover from backup
		else:
			shutil.copy(index_path, bck_index)  # keep a backup
		html = str(soup)
		new_html = html.replace('<head>', '<head>\n' + GA_JS)
		index_path.write_text(new_html)


# inject_gad()
# inject_gad1()
# inject_ga()


h_plchldr = st.empty()

# Amazon affiliates
HtmlFile = open("test_amzn.html", 'r', encoding='utf-8')
source_code = HtmlFile.read()
# print(source_code)
components.html(source_code, height=200)

# Menù laterale
sidemenus = st.sidebar
sidemenus.header('Menu')

dashboard_name = sidemenus.text_input('Name of the project:', 
	help='At the end of the tertiary data analysis two \
		.html dashboards of results could be downloaded. Created with the library `piesparrow` and named: \
			\"Dashboard-microbiome-[Bar/Pie]Chart-<name of the project>-<samples grouping variable>.html\"',
	key='dashboard_name')

if (('dashboard_name' in st.session_state.keys()) and (st.session_state.dashboard_name != '')):
	pass
else:
	st.warning('Start by inserting a name for the project in the sidebar menu')
	st.stop()


side_subheader_placeholder = sidemenus.empty()
side_placeholder0 = sidemenus.empty()
side_placeholder0a = sidemenus.empty()
sidemenus.markdown('---')
sidemenus.markdown('IMPORTANT (more info ?):')
skip_placeholder = sidemenus.empty()
sidemenus.markdown('---')

library_PE_SE_placeholder = sidemenus.empty()
step_n = 0

skip = skip_placeholder.checkbox('Skip the steps of \"Secondary analysis of raw data\" and \
	go directly to the \"Tertiary analysis of pre-processed data\".',
	help='To be selected to do the tertiary data analysis \
		of OTUs/ASVs and taxonomy tables and associated samples\' metadata. \
			To be always selected at the end of all the steps of the secondary analysis of raw data to speed up the App, \
				after downloading the two .csv files of results (\
					OTUs/ASVs rarefied table inside tab _Normalization_ and \
						taxonomy table inside section _Taxonomic classification_) which \
							would need to be uploaded inside the tertiary data analysis of pre-processed data form. App page re-run occurs \
								by default at each user-interaction. To be selected to reduce wait times.')

library_PE_SE = library_PE_SE_placeholder.radio('Library type:', options=['Paired-end', 'Single-end'], 
	help='Illumina library prep method to be considered for the correct raw data upload. Contact raw data provider.',
	key='library_radio')

hypervar_regions_plchldr = sidemenus.empty()
######################## FUNZIONI UTILI ########################
# Funzioni copiate da helper.py file

#@st.cache_resource(show_spinner=True)
def import_paired_end_fastq_gz_files(_filepath):
	'''
	Import R1 and R2 fastq.gz files for all samples in the project.
	fastq files must be in the Casava 1.8 format.
	'''
	with st.spinner('Analisi in corso'):
		paired_end_sequences = qiime2.Artifact.import_data(
			'SampleData[PairedEndSequencesWithQuality]', 
			_filepath,
			'CasavaOneEightSingleLanePerSampleDirFmt')
	return paired_end_sequences


#@st.cache_resource(show_spinner=True)
def import_single_end_fastq_gz_files(_filepath):
	'''
	Import R1 fastq.gz files for all samples in the project.
	fastq files must be in the Casava 1.8 format.
	'''
	with st.spinner('Analisi in corso'):
		single_end_sequences = qiime2.Artifact.import_data(
			'SampleData[SequencesWithQuality]', 
			_filepath,
			'CasavaOneEightSingleLanePerSampleDirFmt')
	return single_end_sequences

#@st.cache_resource(show_spinner=True)
def vsearch_join_pairs(_paired_end_sequences, minovlen, minmergelen, maxmergelen, maxee, maxdiffs):#minovlen=20, minmergelen=460, maxmergelen=480, maxee=1, maxdiffs=5):
	'''
	Performs joining of paired end reads based on given paramenters.
	Default parameters are set for reads' length of 250 bases (v2 Illumina MiSeq kit)
	'''
	with st.spinner('Analisi in corso'):
		demux_joined = vsearch.methods.merge_pairs(
			demultiplexed_seqs = _paired_end_sequences,
			minovlen=minovlen,
			minmergelen=minmergelen,
			maxmergelen=maxmergelen,
			maxee=maxee,
			maxdiffs=maxdiffs)
	return demux_joined



#@st.cache_resource(show_spinner=True)
def quality_filter_paired_end(_demux_joined, min_quality, quality_window):
	'''
	Performs quality filtering of paired-end fastq reads based on phred64 Illumina quality score.
	'''
	with st.spinner('Analisi in corso'):
		demux_filter = quality_filter.methods.q_score(demux=_demux_joined, min_quality=min_quality, quality_window=quality_window)
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
		# cols_renameing = {
		#     "total_input_reads": 'totale_sequenze_iniziali', 
		#     "total_retained_reads": 'totale_sequenze_accettabili', 
		#     "reads_truncated": 'sequenze_troncate', 
		#     "reads_too_short_after_truncation": 'sequenze_troncate_troppo_corte_scartate',
		#     "reads_exceeding_maximum_ambiguous_bases": ('sequenze_con_oltre_%s_basi_ambigue_scartate' %(quality_window)),
		#     'perc_retained_reads': 'percentuale_sequenze_accettabili'
		#     }
		
		# df_q_filter = df_q_filter.rename(cols_renameing, axis=1)
	return demux_filter, df_q_filter, secure_temp_dir_q_filter_summary

#@st.cache_resource
def app_alpha_rare_curves(_table, max_depth, metrics):
	with st.spinner('Analisi in corso'):
		alpha_rare_curves = alpha_rarefaction(table = _table, max_depth=max_depth, metrics=metrics)
	return alpha_rare_curves


#@st.cache_resource
def app_align_to_tree_mafft_fasttree(_sequences, _table, sampling_depth, _metadata):
	with st.spinner('Analisi in corso'):
		result_alignment = align_to_tree_mafft_fasttree(sequences = _sequences)
		seqs_alignment = result_alignment.alignment
		masked_alignment = result_alignment.masked_alignment
		tree = result_alignment.tree
		rooted_tree = result_alignment.rooted_tree

		result_core_metrics = core_metrics_phylogenetic(table=_table, phylogeny=rooted_tree, sampling_depth=sampling_depth, metadata=_metadata)
		core_metr_phylo = result_core_metrics
	return seqs_alignment, masked_alignment, tree, rooted_tree, core_metr_phylo


#@st.cache_resource(show_spinner=True)
def dada2_denoise_single_joined(_demux_filter, N, trim_TF):
	with st.spinner('Analisi in corso'):
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





#@st.cache_resource(show_spinner=True)
def deblur_denoise_trim_paired_end(_demux_filter, N, trim_TF):
	'''
	Performs denoising of data and trimming at position N, which is defined by the user based on
	the filter stats --> position in the reads where the median quality score drops too low.
	Only forward reads are supported at this time, so perform paired-end reads joining first.
	If trim_TF is True, trimming is enabled at length N.
	Otherwise, if trim_TF is False, disable trimming.
	"Deblur operates only on same length sequences. " from the web https://forum.qiime2.org/t/deblur-plugin-error/2172
	'''
	with st.spinner('Analisi in corso'):
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


#@st.cache_resource(show_spinner=True)
def import_SequencesWithQuality(_filepath):
	'''
	Import R1 fastq.gz files for all samples in the project.
	fastq files must be in the Casava 1.8 format.
	'''
	with st.spinner('Analisi in corso'):
		single_end_sequences = qiime2.Artifact.import_data(
			'SampleData[SequencesWithQuality]', 
			_filepath,
			'SingleLanePerSampleSingleEndFastqDirFmt')
	return single_end_sequences


@st.cache_resource
def import_gg_13_8_pre_trained_classifier(filepath):
	'''
	Import reference Green Genes 13/8 otus sequences
	'''
	ref_gg_13_8_pre_trained_classifier = qiime2.Artifact.import_data(
		'TaxonomicClassifier', 
		filepath)
	return ref_gg_13_8_pre_trained_classifier

@st.cache_resource
def import_ref_gg_13_8_otus_seqs(filepath):
	'''
	Import reference Green Genes 13/8 otus sequences
	'''
	ref_gg_13_8_seqs = qiime2.Artifact.import_data(
		'FeatureData[Sequence]', 
		filepath)
	return ref_gg_13_8_seqs

@st.cache_resource
def import_ref_gg_13_8_otus_taxonomy(filepath):
	'''
	Import reference Green Genes 13/8 taxonomy
	'''
	ref_gg_13_8_taxonomy = qiime2.Artifact.import_data(
		'FeatureData[Taxonomy]', 
		filepath)
	return ref_gg_13_8_taxonomy


#########################
# Funzione per aggiungere una casella da spuntare per filtrare un df interattivamente
from pandas.api.types import (
	is_categorical_dtype,
	is_datetime64_any_dtype,
	is_numeric_dtype,
	is_object_dtype,
)

def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
	"""
	Adds a UI on top of a dataframe to let viewers filter columns

	Args:
		df (pd.DataFrame): Original dataframe

	Returns:
		pd.DataFrame: Filtered dataframe
	"""
	modify = st.checkbox("Add filters to table")

	if not modify:
		return df

	df = df.copy()

	# Try to convert datetimes into a standard format (datetime, no timezone)
	for col in df.columns:
		if is_object_dtype(df[col]):
			try:
				df[col] = pd.to_datetime(df[col])
			except Exception:
				pass

		if is_datetime64_any_dtype(df[col]):
			df[col] = df[col].dt.tz_localize(None)

	modification_container = st.container()

	with modification_container:
		to_filter_columns = st.multiselect("Filter the dataframe based on", df.columns)
		for column in to_filter_columns:
			left, right = st.columns((1, 20))
			# Treat columns with < 10 unique values as categorical
			if is_categorical_dtype(df[column]) or df[column].nunique() < 10:
				user_cat_input = right.multiselect(
					f"Values for {column}",
					df[column].unique(),
					default=list(df[column].unique()),
				)
				df = df[df[column].isin(user_cat_input)]
			elif is_numeric_dtype(df[column]):
				_min = float(df[column].min())
				_max = float(df[column].max())
				step = (_max - _min) / 100
				user_num_input = right.slider(
					f"Values for {column}",
					min_value=_min,
					max_value=_max,
					value=(_min, _max),
					step=step,
				)
				df = df[df[column].between(*user_num_input)]
			elif is_datetime64_any_dtype(df[column]):
				user_date_input = right.date_input(
					f"Values for {column}",
					value=(
						df[column].min(),
						df[column].max(),
					),
				)
				if len(user_date_input) == 2:
					user_date_input = tuple(map(pd.to_datetime, user_date_input))
					start_date, end_date = user_date_input
					df = df.loc[df[column].between(start_date, end_date)]
			else:
				user_text_input = right.text_input(
					f"Substring or regex in {column}",
				)
				if user_text_input:
					df = df[df[column].astype(str).str.contains(user_text_input)]

	return df

# Funzione per la creazione di un archivio .zip di una cartella mantenendo la struttura delle sotto cartelle
def zipfolder(foldername, target_dir):            
	zipobj = zipfile.ZipFile(foldername, 'w', zipfile.ZIP_DEFLATED)
	rootlen = len(target_dir) + 1
	for base, dirs, files in os.walk(target_dir):
		for file in files:
			fn = os.path.join(base, file)
			zipobj.write(fn, fn[rootlen:])

# funzione per la visualizzazione di file PDF in streamlit
def displayPDF(file):
	# Opening file from file path
	with open(file, "rb") as f:
		base64_pdf = base64.b64encode(f.read()).decode('utf-8')
	# Embedding PDF in iframe
	pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="700" height="600" type="application/pdf"></iframe>'
	return pdf_display

# funzione per convertire un pandas dataframe a file csv con encoding utf-8, per il download
@st.cache_data
def convert_df(df): # sembra che non funzioni! Usata libreria streamlit_ext al posto della funzione nativa: ste.download_button
	# IMPORTANT: Cache the conversion to prevent computation on every rerun
	return df.to_csv().encode('utf-8')

# funzione per arrotondare un numero alla decina piu vicina
def myround(x, base=10):
	return base * round(x/base)


# funzione cached per calcolare le statistiche riassuntive dei dati grezzi caricati dall'utente
#@st.cache_resource(show_spinner=True)
def app_demux_visualizers_summarize(_sequences):
	'''
	Funzione cached per calcolare le statistiche riassuntive 
	dei dati grezzi caricati dall'utente
	'''
	with st.spinner('Analisi in corso'):
		demux_summary = demux.visualizers.summarize(_sequences)
		secure_temp_dir_demux_summary = tempfile.mkdtemp(prefix="temp_", suffix="_demux_summary")
		demux_summary.visualization.export_data(secure_temp_dir_demux_summary)
		demux_summary = pd.read_csv(secure_temp_dir_demux_summary+'/per-sample-fastq-counts.tsv', sep='\t' , index_col='sample ID')
		
		with open(secure_temp_dir_demux_summary+"/quality-plot.html", 'r', encoding='utf-8') as HtmlFile:
			source_code_demux_summary = HtmlFile.read() 

	return demux_summary, source_code_demux_summary, secure_temp_dir_demux_summary

# funzione del bottone per cancellare i file caricati in precedenza e la cache associata e 
# i file temporanei associati
def clear_files_and_cache_button_callback():
	'''
	Funzione del bottone per cancellare i file caricati in precedenza e 
	la cache associata e i file temporanei associati, da cliccare per 
	modificare i dati da analizzare. Cancella la cache di:
	import_function (import_paired_end_fastq_gz_files o import_single_end_fastq_gz_files);
	app_demux_visualizers_summarize();
	cancella la st.session_state_key demux_fastq_input;
	cancella le cartelle temporanee imported_sequences_temp_dir e
	secure_temp_dir_demux_summary;
	'''
	try:

		import_paired_end_fastq_gz_files.clear()
		import_function.clear()
		app_demux_visualizers_summarize.clear()
		for i in st.session_state.keys():
			del i
		# del st.session_state['demux_fastq_input']
		# del st.session_state['library_radio']

	except Exception as e:

		st.exception(e)
		st.stop()
	
	try:
		shutil.rmtree(st.session_state.imported_sequences_temp_dir)
	except FileNotFoundError as e:
		st.exception(e)
		
	except NameError as e:
		st.exception(e)
		

# Funzione per effettuare la fusione delle reads paired-end con vsearch e visualizzare i risultati
#@st.cache_resource(show_spinner=True)
def form_callback_vsearch_join(_paired_end_sequences, imported_sequences_temp_dir, secure_temp_dir_demux_summary):
	with st.spinner('Analisi in corso'):
		minovlen = st.session_state.minovlen_slider
		minmergelen = st.session_state.minmergelen_slider
		maxmergelen = st.session_state.maxmergelen_slider
		maxee = st.session_state.maxee_slider
		maxdiffs = st.session_state.maxdiffs_slider
		
		# todo

		demux_joined = vsearch_join_pairs(_paired_end_sequences, minovlen, minmergelen, maxmergelen, maxee, maxdiffs)
		err = 'test'

		secure_temp_dir_joined_summary = tempfile.mkdtemp(prefix="temp_", suffix="_demux_joined_summary")
		demux_joined_summary = demux.visualizers.summarize(demux_joined.merged_sequences)
		demux_joined_summary.visualization.export_data(secure_temp_dir_joined_summary)
		df_demux_joined_summary = pd.read_csv(secure_temp_dir_joined_summary+'/per-sample-fastq-counts.tsv', sep='\t', index_col='sample ID')
		joined_sequences = demux_joined.merged_sequences
		pdf_joined_sequences = displayPDF(secure_temp_dir_joined_summary+'/demultiplex-summary-forward.pdf')
		
	return df_demux_joined_summary, joined_sequences, pdf_joined_sequences, secure_temp_dir_joined_summary, err

# Funzione per creare il dataframe finale dai files OTU table e taxonomy table .tsv
@st.cache_data
def create_final_df(x,y):

	'''
	Funzione per importare OTU .tsv file e taxonomy .tsv file come pandas
	dataframes e per creare il dataframe complessivo di entrambi 
	fusi insieme. Rinomina le colonne dei livelli tassonomici in italiano
	e ritorna 3 pd.DataFrame (final_df, data_tax, data_OTU) ed
	1 lista (taxa_levels).
	indici consentiti taxonomy file: OTU, OTU ID, Feature ID, #OTU ID, features, #NAME, featureID;
	indici consentiti OTU file: #TAXONOMY;
	categorie tassonomiche consentite: {
		'Domain': 'Regno',
		'Kingdom': 'Regno',
		'Class': 'Classe',
		'Order': 'Ordine',
		'Family': 'Famiglia',
		'Genus': 'Genere',
		'Species': 'Specie',
		'feature': 'Variante'
		}
	'''
	
	# Importazione dati come pandas dataframes
	if y.name.endswith('.tsv'):		
		data_tax = pd.read_csv(y, sep='\t', index_col=0) # indici consentiti: OTU, OTU ID, Feature ID, #OTU ID, features, #NAME, featureID
		data_tax.index.name = '#TAXONOMY'
	elif y.name.endswith('.txt'):		
		data_tax = pd.read_csv(y, sep='\t', index_col=0) # indici consentiti: OTU, OTU ID, Feature ID, #OTU ID, features, #NAME, featureID
		data_tax.index.name = '#TAXONOMY'
	elif y.name.endswith('.csv'):
		data_tax = pd.read_csv(y, sep=',', index_col=0) # indici consentiti: OTU, OTU ID, Feature ID, #OTU ID, features, #NAME
		data_tax.index.name = '#TAXONOMY'
	elif y.name.endswith('.qza'): # File in cui sono presenti 3 colonne (moving pictures tutorial): Feature ID, Taxon, Confidence
		with NamedTemporaryFile(dir='.', suffix='.qza') as f:
			f.write(y.getbuffer())
			y_import = Artifact.load(f.name)
			data_tax = y_import.view(pd.DataFrame)
			# Split dela colonna Taxon in base al punto e virgola per creare i livelli tassonomici
			data_tax_tmp = data_tax['Taxon'].str.split(';', expand=True)
			if len(data_tax_tmp.columns) == 6:
				data_tax[['Kingdom',
					'Phylum',
					'Class',
					'Order',
					'Family',
					'Genus']] = data_tax_tmp
				data_tax = data_tax[['Kingdom',
					'Phylum',
					'Class',
					'Order',
					'Family',
					'Genus']]
			elif len(data_tax_tmp.columns) == 7:
				data_tax[['Kingdom',
					'Phylum',
					'Class',
					'Order',
					'Family',
					'Genus',
					'Species']] = data_tax_tmp
				# Mantenimento delle sole colonne nuove create sopra
				data_tax = data_tax[['Kingdom',
					'Phylum',
					'Class',
					'Order',
					'Family',
					'Genus',
					'Species']]
			data_tax.index.name = '#TAXONOMY'


	data_tax = data_tax.fillna('Unclassified')

	
	# Se presente una colonna ulteriore oltre la Specie, come OTU o feature, diventa l indice del dataframe
	try:

		data_tax['OTU'] = data_tax.index
	except Exception as e:

		st.exception(e)
		try:

			data_tax['Feature ID'] = data_tax.index
		except Exception as e:

			st.exception(e)
			try:

				data_tax['OTU ID'] = data_tax.index
			except Exception as e:

				st.exception(e)
				try:

					data_tax['#OTU ID'] = data_tax.index
				except Exception as e:

					st.exception(e)
					try:

						data_tax['features'] = data_tax.index
					except Exception as e:

						st.exception(e)
						try:

							data_tax['featuresID'] = data_tax.index
						except Exception as e:

							st.exception(e)
						
	data_tax.index.name = '#TAXONOMY'			

	# Rinominazione delle colonne dei livelli tassonomici in italiano:
	# final_df_taxonomy_level_ITA_col_names = {
	# 	'Domain': 'Regno',
	# 	'Kingdom': 'Regno',
	# 	'Class': 'Classe',
	# 	'Order': 'Ordine',
	# 	'Family': 'Famiglia',
	# 	'Genus': 'Genere',
	# 	'Species': 'Specie',
	# 	'feature': 'Variante'
	# 	}
	# data_tax = data_tax.rename(final_df_taxonomy_level_ITA_col_names, axis=1)
	
	if x.name.endswith('.tsv'):		
		data_OTU = pd.read_csv(x, sep='\t', index_col=0) # indici consentiti: #NAME
	elif x.name.endswith('.txt'):		
		data_OTU = pd.read_csv(x, sep='\t', index_col=0) # indici consentiti: #NAME
	elif x.name.endswith('.csv'):		
		data_OTU = pd.read_csv(x, sep=',', index_col=0) # indici consentiti: #NAME
	elif x.name.endswith('.qza'):		
		with NamedTemporaryFile(dir='.', suffix='.qza') as f:
			f.write(x.getbuffer())
			x_import = Artifact.load(f.name)
			data_OTU = x_import.view(pd.DataFrame).T
			
	
	# Creazione multiindex df comprensivo dei dati obbligatori 
	# (classificaizone tassonomica ed conte delle reads per OTU per campione):
	# Definizione dei livelli della tassonomia così come dai dati importati
	taxa_levels = data_tax.columns.to_list()
	# Creazione di una copia del df originale della classificazione tassonomica delle OTU,
	# per modificare il formato in un multiindex pandas dataframe
	data_tax_multi = data_tax.copy().reset_index()
	data_tax_multi = pd.DataFrame(data_tax.index)
	# Creazione di un multiindex index gerarchico dei livelli tassonomici
	data_tax_multi.index = pd.MultiIndex.from_frame(data_tax.loc[:], names=data_tax.columns)
	# Dataframe definitivo: fusione del multiindex df della classificazione tassonomica 
	# delle OTU con i dati originali dell'OTU file (reads counts per ogni OTU per ogni campione)
	try:
		final_df = data_tax_multi.merge(data_OTU, left_on='#TAXONOMY', right_index=True, how='inner')
		key_error_warning = ''
	except KeyError as e:
		key_error_warning = 'Sequences of OTU file and taxonomy file do not match %s' %e
		final_df = pd.concat([data_tax_multi,data_OTU])
		final_df[final_df.duplicated('#TAXONOMY', keep=False)]

	otus_col = (final_df.iloc[:,0]).values
	final_df = final_df.iloc[:,1:]
	final_df = final_df.reset_index()
	final_df['OTU'] = otus_col
	
	return final_df, data_tax, data_OTU, taxa_levels, key_error_warning

# estrazione di dataframe pandas da file PDF 
# NON SERVE
def extract_sequence_len_stats_from_pdf(source_code):
	'''
	Funzione per estrarre in formato pandas dataframe
	le informazioni dal file pdf di sommario 
	riguardo alla distribuzione delle lunghezze 
	delle sequenze in numero di nucleotidi.
	'''
	fc = 28.28
	box = [1,1,20,1]
	
	for i in range(0, len(box)):
		box[i] *= fc

	tl = tb.read_pdf(source_code,
		pages = 1,
		area = [box],
		output_format = "dataframe",
		pandas_options = {'header': None},
		stream=True)

	df = tl[0]

	return df

# Annotazione delle OTU presenti in ogni gruppo di campioni e frequenza delle annotazioni.
@st.cache(show_spinner=True)
def OTUs_annots_freqs(_idx_sample_grouping):
	'''
	Annotazione delle OTU presenti in ogni gruppo di campioni 
	e frequenza delle annotazioni.
	Funzione che ritorna un dizionario di liste di dizionari ed un dizionario di pandas dataframe e 
	stampa a schermo: il nome di ogni gruppo di campioni,
	un menu expander denominato "1. nome del primo taxon del livello tassonomico selezionato,
	2. nome del secondo taxon del livello tassonomico selezionato 3. ..." e cosi via
	con al suo interno un grafico a barre delle OTUs associate a tale taxon ed
	una tabella di tutte le OTUs associate a tale taxon, le stesse rappresentate nel grafico a barre.
	La funziona ritorna un dizionario di liste di dizionari e un pandas dataframe:
	taxa_counts_d: {chiave = nome di ciascun gruppo di campioni:  valore = lista di dizionari con [{
		chiave = ogni taxon della lista taxa - contenente tutti quelli osservati nel gruppo di campioni 
		con conte maggiori di zero per il livello tassonomico selezionato - : valore = ogni numero di OTUs
		della lista counts - contenente tutti i numeri di OTUs associati ad ogni taxon per il 
		livello tassonomico selezionato -}
	df1_d: dizionario di dataframes, uno per ciascun gruppo di campioni (chiavi), con 
	colonne = campioni presenti nel gruppo selezionato e livello tassonomico selezionato, 
	righe = OTUs con conte maggiori di zero per il livello tassonomico selezionato ed il livello Phylum
	'''
	if st.session_state.tax_level_radio != 'OTU':
	
		warning_tab_num_of_OTUs_per_taxon = None
		taxa_counts_d = {}
		df1_d = {}
		tabella_df_l_l = []
		for grouping_key, grouped_samples in _idx_sample_grouping.items():
			
			# le righe commentate di stampa a schermo nella applicazione sono spostate fuori dalla funzione, dopo la sua chiamata, per poter sfruttare i vantaggi del caching delle variabili
			# st.markdown('***%s***' %grouping_key)
			
			if type(grouped_samples) == 'object': # se Raggruppamento dei campioni = Tutti i campioni, allora grouped_samples e' gia' una lista. Negli altri casi e' un object type di classe pd.Index
				grouped_samples = list(grouped_samples)
			

			# Rimozione delle OTU che non sono presenti (hanno reads counts = 0) nel gruppo di campioni
			otus_idxs = final_df.loc[:, final_df.columns.isin(grouped_samples)][final_df.loc[:, final_df.columns.isin(grouped_samples)].any(axis=1)].index
			df = pd.DataFrame(final_df.iloc[otus_idxs,:].loc[:,st.session_state.tax_level_radio])
			
			# per la tab grafico delle abbondanze relative
		
			if ((st.session_state.tax_level_radio != 'Phylum') and (st.session_state.tax_level_radio != 'Kingdom')):
				df1 = pd.DataFrame(final_df.iloc[otus_idxs,:][(pd.Index(['Phylum'])).append((pd.Index(grouped_samples))).append((pd.Index([st.session_state.tax_level_radio])))])
				df1_d[grouping_key] = df1
			elif st.session_state.tax_level_radio == 'Phylum':
				df1 = pd.DataFrame(final_df.iloc[otus_idxs,:][(pd.Index(['Kingdom'])).append((pd.Index(grouped_samples))).append((pd.Index([st.session_state.tax_level_radio])))])
				df1_d[grouping_key] = df1
			elif st.session_state.tax_level_radio == 'Kingdom':
				df1 = pd.DataFrame(final_df.iloc[otus_idxs,:][(pd.Index(grouped_samples)).append((pd.Index([st.session_state.tax_level_radio])))])
				df1_d[grouping_key] = df1
			
			# Frequenza delle annotazioni tassonomiche
			taxa, counts = np.unique(df[st.session_state.tax_level_radio].values, return_counts=True)
			taxa_counts_d[grouping_key] = [{i: j} for i, j in zip(taxa, counts)]
			
			tabella_df_l = []
			for number, taxon in enumerate(taxa):
				
				expander_string = ((' '.join([str(number+1), '. ', taxon])))
				
				with st.expander(expander_string):

					otus = (data_tax.iloc[otus_idxs,:][data_tax.iloc[otus_idxs,:][st.session_state.tax_level_radio] == taxon].index)
					tabella_df = pd.DataFrame(data = otus.values, columns = [taxon], index = range(1,len(otus)+1))
					tabella_df.columns.name = taxon
	
					tabella_df_l.append(tabella_df)
				
					# le righe commentate di stampa a schermo nella applicazione sono spostate fuori dalla funzione, dopo la sua chiamata, per poter sfruttare i vantaggi del caching delle variabili
					# st.bar_chart(tabella_df)

					# st.dataframe(tabella_df)
					
					# csv = convert_df(tabella_df)

					# ste.download_button(
					# 	label="Scarica tabella in formato CSV",
					# 	data=csv,
					# 	file_name=('OTUs_per_%s_%s.csv' %(st.session_state.tax_level_radio, taxon)),
					# 	mime='text/csv')
			
			tabella_df_l_l.append(tabella_df_l)

	elif st.session_state.tax_level_radio == 'OTU':
		
		tabella_df_l_l = None
		taxa_counts_d = {}
		df1_d = {}
	
		for grouping_key, grouped_samples in _idx_sample_grouping.items():
			
			# le righe commentate di stampa a schermo nella applicazione sono spostate fuori dalla funzione, dopo la sua chiamata, per poter sfruttare i vantaggi del caching delle variabili
			# tab_num_of_OTUs_per_taxon.markdown('***%s***' %grouping_key)

			grouped_samples = list(grouped_samples)

			# Rimozione delle OTU che non sono presenti (hanno reads counts = 0) nel gruppo di campioni
			otus_idxs = final_df.loc[:, final_df.columns.isin(grouped_samples)][final_df.loc[:, final_df.columns.isin(grouped_samples)].any(axis=1)].index
			df = pd.DataFrame(final_df.iloc[otus_idxs,:].loc[:,st.session_state.tax_level_radio])

			# per la tab grafico delle abbondanze relative
			df1 = pd.DataFrame(final_df.iloc[otus_idxs,:][['Genus'] + grouped_samples + [st.session_state.tax_level_radio]])
			
			df1_d[grouping_key] = df1
			warning_tab_num_of_OTUs_per_taxon = 'It\'s not possible to visualize charts at the OTU level.'

			# Frequenza delle annotazioni tassonomiche
			taxa, counts = np.unique(df[st.session_state.tax_level_radio].values, return_counts=True)
			taxa_counts_d[grouping_key] = [{i: j} for i, j in zip(taxa, counts)]
			
	return taxa_counts_d, df1_d, tabella_df_l_l, warning_tab_num_of_OTUs_per_taxon

	
#@st.cache_resource
def app_classify_hybrid_vsearch_sklearn(_query, _reference_reads, _reference_taxonomy, _classifier, reads_per_batch, randseed):
	with st.spinner('Analisi in corso'):
		seqs_classification = classify_hybrid_vsearch_sklearn(query=_query,
			reference_reads=_reference_reads,
			reference_taxonomy=_reference_taxonomy,
			classifier=_classifier,
			prefilter=False,
			reads_per_batch=reads_per_batch,
			randseed=randseed)
	return seqs_classification

@st.cache
def app_alpha_divs(_table):
	alpha_divs = []
	alpha_result = diversity.pipelines.alpha(table=_table, metric='observed_features')
	alpha_diversity_obs_feat = alpha_result.alpha_diversity
	alpha_result = diversity.pipelines.alpha(table=_table, metric='shannon')
	alpha_diversity_shannon = alpha_result.alpha_diversity
	alpha_result = diversity.pipelines.alpha(table=_table, metric='simpson')
	alpha_diversity_simpson = alpha_result.alpha_diversity
	alpha_result = diversity.pipelines.alpha(table=_table, metric='pielou_e')
	alpha_diversity_pielou_e = alpha_result.alpha_diversity
	
	alpha_divs.append(alpha_diversity_obs_feat)
	alpha_divs.append(alpha_diversity_shannon)
	alpha_divs.append(alpha_diversity_pielou_e)
	alpha_divs.append(alpha_diversity_simpson)

	return alpha_divs
################################################################

if skip is False:

	step_n += 1
	side_subheader_placeholder.subheader('***Secondary analysis raw data***')
	side_placeholder0.info('***%s.*** Data upload' %(step_n))
	side_placeholder0a.markdown(f"<a href='#linkto_{step_n}'>{step_n}. Data upload</a>", unsafe_allow_html=True)

	h_plchldr.header('***Secondary analysis of raw data***')

	# sample_data_secondary_analysis = '/app/microbiome/sample_data/paire_end_sequences'
	# st.info('Download sample raw data for secondary analysis.')
	
	# with open(sample_data_secondary_analysis+"/zip_sample_data.zip", 'rb') as f:
		
	# 	ste.download_button(
	# 		label="Download sample raw data .zip",
	# 		data=f,
	# 		file_name="sample_raw_data_secondary_analysis.zip",
	# 		mime="application/zip")

	# Form di caricamento dati
	st.markdown(f"<div id='linkto_{step_n}'></div>", unsafe_allow_html=True)
	ctx = get_script_run_ctx()
	session_id = ctx.session_id
	with st.form(key='form_demux_fastq_upload', clear_on_submit=False):
		if st.session_state.library_radio == 'Single-end':
			
			library_radio_help_string = 'Upload a R1 file for each sample. \
				\n> Filename format: [SampleID]\_[index]\_[L001]\_R1\_[001].fastq.gz (e.g. Sample1\_S1\_L001\_R1\_001.fastq.gz)'
			import_function = import_single_end_fastq_gz_files
			df_cols_to_rename_for = {'forward sequence count': 'forward sequence count'}
			df_cols_to_rename_rev = {'reverse sequence count': 'forward sequence count'}
			sample_data_dir = '/app/microbiome/sample_data/single_end_sequences'
		elif st.session_state.library_radio == 'Paired-end':
			
			library_radio_help_string = 'Upload a R1 and a R2 file for each sample.\
				\n> Filename format: [SampleID]\_[index]\_[L001]\_R[1-2]\_[001].fastq.gz (e.g. Sample1\_S1\_L001\_R1\_001.fastq.gz)'
			import_function = import_paired_end_fastq_gz_files
			df_cols_to_rename = {'forward sequence count': 'forward sequence count', 'reverse sequence count': 'reverse sequence count'}
			sample_data_dir = '/app/microbiome/sample_data/paire_end_sequences'


		lbl_plchldr = st.empty()

		data_demux_fastq = st.file_uploader(
			label = 'Files fastq.gz:',
			key='demux_fastq_input',
			accept_multiple_files = True,
			help='Upload raw data files to be analyzed. Accepted file formats: fastq.gz, Phred + 33. \
				\n> Illumina library prep method selected in the sidebar menu is __%s__: %s.' %(st.session_state.library_radio, library_radio_help_string))

		
		submit_button = st.form_submit_button(
			label='Upload',  type="primary"
			)
		
		clear_files_and_cache_bttn_plchldr = st.empty()
		# bottone per cancellare i file caricati in precedenza e la cache associata, da cliccare per modificare i files fastq.gz
		clear_files_and_cache_button = clear_files_and_cache_bttn_plchldr.form_submit_button('Modify data files')
	
		sample_data_bttn_plchldr = st.empty()
		# bottone per caricare dati di esempio paired-end
		sample_data_bttn = sample_data_bttn_plchldr.form_submit_button('Load sample data')
	
	if submit_button or (st.session_state['demux_fastq_input'] != []):
		
		try:
			
			del st.session_state.imported_sequences
		except Exception as e:
			pass
		try:

			imported_sequences_temp_dir = tempfile.mkdtemp(prefix="temp_",suffix="_fastq_gz")
			st.session_state.imported_sequences_temp_dir = imported_sequences_temp_dir
			for i in st.session_state['demux_fastq_input']:
				with tempfile.NamedTemporaryFile(dir=st.session_state.imported_sequences_temp_dir, prefix=i.name, suffix='.fastq.gz', delete=False) as f:
					f.write(i.getbuffer())
				os.rename(f.name, '%s/%s' %(st.session_state.imported_sequences_temp_dir, i.name))
				
			imported_sequences = import_function(st.session_state.imported_sequences_temp_dir)
			st.session_state.imported_sequences = imported_sequences

			st.success('Successfully uploaded raw data.')
			lbl_plchldr.info('***ATTENTION!*** To modify previously uploaded files and re-run an analysis from scratch: \
				\n> :one: manually delete all files by clicking on each :x: ; \
				\n> :two: __Select or drag and drop__ new fastq.gz files; \
				\n> :three: click __Modify data files__: the page will be re-ran and new files kept in memory; \
				\n> :four: click Upload')
			
			
				
			
			side_placeholder0.success('***%s.*** Data upload \
				\n> Tab %s. Raw data summary statistics' %(step_n, step_n))
			side_placeholder0a.markdown(f"<a href='#linkto_{step_n}_tab'>Tab {step_n}. Raw data summary statistics</a>", unsafe_allow_html=True)

			


			
			if clear_files_and_cache_button:
				clear_files_and_cache_button_callback()
				st.experimental_rerun()
				
			try:
				
				import_function.clear()
				app_demux_visualizers_summarize.clear()
				quality_filter_paired_end.clear()
				import_SequencesWithQuality.clear()
				app_classify_hybrid_vsearch_sklearn.clear()
				#del st.session_state['demux_fastq_input']
			except Exception as e:
				st.exception(e)
			
			descr_plchldr = st.empty()
			descr_plchldr.markdown('The main page shows results and is to be scrolled. The tabs are to be navigated in order \
			to set the parameters of the raw data analysis steps.')
			
			st.markdown('***You should navigate the tabs below in order to execute step-by-step \
				raw data pre-processing.*** ')
		
		except Exception as e:

			if e.__class__ is qiime2.core.exceptions.ValidationError:
				st.warning('Sequencing files must match Illumina library prep method \
					selected in the sidebar menu (%s). \
					For more infos see the ? in the upper right in the above form.' %(library_PE_SE))
			
			st.stop()

		
		st.balloons()
	elif sample_data_bttn:

		imported_sequences = import_function(sample_data_dir)
		st.session_state.imported_sequences = imported_sequences
		st.session_state.imported_sequences_temp_dir = sample_data_dir
		
		st.success('Successfully uploaded raw data.')

		side_placeholder0.success('***%s.*** Data upload \
				\n> Tab %s. Raw data summary statistics' %(step_n, step_n))
		side_placeholder0a.markdown(f"<a href='#linkto_{step_n}_tab'>Tab {step_n}. Raw data summary statistics</a>", unsafe_allow_html=True)

		st.balloons()

	elif (("imported_sequences" in st.session_state.keys()) and (st.session_state.imported_sequences is not None)):

		pass
	else:
		st.stop()

	
		
	if st.session_state.library_radio == 'Single-end':
		# Definizione di 4 tabs per la visualizzazione degli step di pre-processamento dei dati grezzi
		tab_summary_stats, tab_QC, tab_denoising, tab_taxonomy, tab_rarefaction, tab_diversity, tab_phylogenetic_tree = st.tabs([':one: Raw data summary statistics', 
		':two: Quality control', 
		':three: Denoising', 
		':four: Taxonomic classification',
		':five: Normalization',
		':six: Alpha and beta diversity\'',
		':seven: Phylogenetic tree'])
	elif st.session_state.library_radio == 'Paired-end':
		# Definizione di 5 tabs per la visualizzazione degli step di pre-processamento dei dati grezzi
		tab_summary_stats, tab_join, tab_QC, tab_denoising, tab_taxonomy, tab_rarefaction, tab_diversity, tab_phylogenetic_tree = st.tabs([':one: Raw data summary statistics', 
		':two: R1 and R2 merging', 
		':three: Quality control', 
		':four: Denoising',
		':five: Taxonomic classification', 
		':six: Normalization',
		':seven: Alpha and beta diversity\'',
		':eight: Phylogenetic tree'
		]
		)
	

	# Calcolo delle statistiche riassuntive dei dati grezzi con metodo qiime2 demux.visualizers.summarize()
	try:
		demux_summary, source_code_demux_summary, secure_temp_dir_demux_summary = app_demux_visualizers_summarize(st.session_state.imported_sequences)
		st.session_state.secure_temp_dir_demux_summary = secure_temp_dir_demux_summary
	except Exception as e:
		st.exception(e)
		st.error('Error during computing of summary statistics of raw data')
		st.stop()
	st.balloons()

	
	with tab_summary_stats:
		st.markdown(f"<div id='linkto_{step_n}_tab'></div>", unsafe_allow_html=True)
		step_n += 1
		try: 
			try: # R1
				df = demux_summary.rename(df_cols_to_rename_for, axis=1, errors = "raise")
			except: # R2
				df = demux_summary.rename(df_cols_to_rename_rev, axis=1)
		except: # paired-end
			df = demux_summary.rename(df_cols_to_rename, axis=1, errors = "raise")
			
		st.session_state.demux_summary_df = df
		st.subheader('Sequence count per sample')
		try:
			st.write('Total samples: %s R1, %s R2' %(
			len(~df['forward sequence count'].isna()),
			len(~df['reverse sequence count'].isna())))
		except:
			st.write('Total samples: %s R1' %(
			len(~df['forward sequence count'].isna())))
		
		tabella_df = df
		st.table(tabella_df.style.format(formatter='{:,.0f}'))

		csv = convert_df(tabella_df)
		ste.download_button(
			label="Download CSV table",
			data=csv,
			file_name='sequence_counts_per_sample.csv',
			mime='text/csv')

		st.table(pd.DataFrame(df.sum().rename('Total')).T.style.format(formatter='{:,.0f}'))
		
		cols = st.columns((1,1))
		forw_subheader = cols[0].empty()
		rev_subheader = cols[1].empty()
		# Displaying pdf Files of demux summary visualization
		# definizione del numero di bins per gli istogrammi basato sulla quantita di
		# campioni caricati
		if ((len(df.index) > 4) and (len(df.index) < 10)):
			bins = 4
		elif ((len(df.index) > 1) and (len(df.index) <= 4)):
			bins = 2
		elif (len(df.index) == 1):
			bins = 1
		elif (len(df.index) > 10):
			bins = 10
		try:
			forw_subheader.subheader('Forward reads frequency histogram')
			fig, ax = plt.subplots()
			ax = df['forward sequence count'].plot.hist(bins=bins, color='grey', linewidth=1, edgecolor='black')
			ax.set_xlabel('Number of sequences')
			ax.set_ylabel('Number of samples')

			cols[0].pyplot(fig)

			#cols[0].markdown(displayPDF(secure_temp_dir+'/demultiplex-summary-forward.pdf'), unsafe_allow_html=True)
		except Exception as e:

			st.exception(e)
			pass

		try:
			rev_subheader.subheader('Reverse reads frequency histogram')
			fig, ax = plt.subplots()
			ax = df['reverse sequence count'].plot.hist(bins=bins, color='grey', linewidth=1, edgecolor='black')
			ax.set_xlabel('Number of sequences')
			ax.set_ylabel('Number of samples')

			cols[1].pyplot(fig)

			#cols[1].markdown(displayPDF(secure_temp_dir+'/demultiplex-summary-reverse.pdf'), unsafe_allow_html=True)
		except:
			rev_subheader.empty()
			pass

		
		st.subheader('Sequence counts summary')
		tabella_df = df.describe()
		# .rename(
		# 	{'count': 'num. campioni', 
		# 	'mean': 'media',
		# 	'std': 'dev. standard',
		# 	'min': 'minimo',
		# 	'max': 'massimo'}, axis=0).drop('num. campioni', axis=0)
		st.table(tabella_df.style.format(formatter='{:,.2f}'))

		csv = convert_df(tabella_df)
		ste.download_button(
			label="Download CSV table",
			data=csv,
			file_name='sequenze_counts_summary.csv',
			mime='text/csv')

		st.subheader('Reads length distribution summary')
		source_code = source_code_demux_summary
		
		dff_R1= pd.DataFrame(pd.read_html(source_code,thousands='.')[0])
		# .replace({
		# 	'Total Sequences Sampled': 'Totale sequenze campionate', 
		# 	'50% (Median)': '50% (Mediana)'}
		# 	)
		dff_R1 = dff_R1.rename({'Unnamed: 0': 'index', 'Unnamed: 1': 'Forward reads'}, axis=1).set_index('index')
		try:
		  
			dff_R2 = pd.DataFrame(pd.read_html(source_code,thousands='.')[1])
			# .replace({
			# 	'Total Sequences Sampled': 'Totale sequenze campionate', 
			# 	'50% (Median)': '50% (Mediana)'}
			# 	)
			dff_R2 = dff_R2.rename({'Unnamed: 0': 'index', 'Unnamed: 1': 'Reverse reads'}, axis=1).set_index('index')
		
			tabella_dff = pd.merge(dff_R1, dff_R2, left_index=True, right_index=True)
		except Exception as e:
			
			st.warning('Single-end sequences')
			tabella_dff = dff_R1

		st.table(tabella_dff)
		
		approx_reads_lenght_R1 = myround(int(tabella_dff.iloc[5,0].split(' ')[0]))
		
		try:
			approx_reads_lenght_R2 = myround(int(tabella_dff.iloc[5,1].split(' ')[0]))
		except:
			approx_reads_lenght_R2 = 0

		#st.dataframe(extract_sequence_len_stats_from_pdf()) # Non serve, risolto col codice qui sopra!!
	st.balloons()
		
	try:

		with tab_join:
			
			side_placeh1 = sidemenus.empty()
			side_placeh1a = sidemenus.empty()
			side_placeh1.info('***%s.*** Set R1 and R2 reads merging parameters \
				\n > Tab %s. R1 and R2 merging' %(step_n, step_n))
			side_placeh1a.markdown(f"<a href='#linkto_{step_n}_tab'>Tab {step_n}. R1 and R2 merging</a>", unsafe_allow_html=True)
			st.markdown(f"<div id='linkto_{step_n}_tab'></div>", unsafe_allow_html=True)
			
			hypervar_regions_form = hypervar_regions_plchldr.form('hypervar_regiorns_form')
			with hypervar_regions_form:
				hypervar_regions = st.radio(label='hypervariable 16S rRNA gene regions target of sequencing:',
				options=['V3V4', 'V4'],
				index=0,
				help='Regions target of sequencing to be considered to determine expected fragment length.',
				key='hypervar_regions_radio')
				hypervar_btn = st.form_submit_button('Confirm')
			
			if ((hypervar_btn) or ('hypervar_regions_radio' in st.session_state.keys())):
				pass
			else:
				st.warning('The page is awaiting for you to select a sequencing target region in the sidebar menu.')
				st.stop()

			st.markdown('Suggested parameters for R1 and R2 reads merging:')
			cols = st.columns((1,1))

			# Blocco di codice per impostare i suggerimenti dei parametri di fusione delle reads R1 ed R2 in base alla lunghezza delle letture (cartucce
			# e tecnologia di sequenziamento)
			if ((approx_reads_lenght_R1 == 250) or (approx_reads_lenght_R1 == 300)):

				reads_lenght_R1_input = cols[0].number_input('R1 reads length:', help='Needed to suggest the correct parameters \
					for R1 and R2 reads merging to obtain the expected fragment. \
						\n> Default approximate value based on previous analysis. \
						This value varies based on the sequencing technology. \
						Contact raw data provider.',
						min_value=250, max_value=300, step=50, value=approx_reads_lenght_R1, key='reads_lenght_R1_n_input')
			else:

				st.session_state.reads_lenght_R1_n_input = approx_reads_lenght_R1
				cols[0].markdown('Approximate R1 reads length: %s' %(st.session_state.reads_lenght_R1_n_input))
			
			if ((approx_reads_lenght_R2 == 250) or (approx_reads_lenght_R2 == 300)):

				reads_lenght_R2_input = cols[1].number_input('R2 reads length:', help='R2 reads length should be \
					equal to that of R1 reads. This value is only showed for precision check.',
						min_value=250, max_value=300, step=50, value=approx_reads_lenght_R2, key='reads_lenght_R2_n_input', disabled=True)
			else:

				st.session_state.reads_lenght_R2_n_input = approx_reads_lenght_R2
				cols[1].markdown('Approximate R2 reads length: %s' %(st.session_state.reads_lenght_R2_n_input))

			

			if st.session_state.reads_lenght_R1_n_input == 250:
				if st.session_state.hypervar_regions_radio == 'V3V4':
					minovlen_n = 20
					minmergelen_n = 440
					maxmergelen_n = 480 # average 464
					maxee_n = 1
					maxdiffs_n = 5
				elif st.session_state.hypervar_regions_radio == 'V4':
					minovlen_n = 220
					minmergelen_n = 230
					maxmergelen_n = 270 # average 254
					maxee_n = 1
					maxdiffs_n = 20
				
				suggested_joining_params_label = ' \
					\n> __MiSeq V2 flowcells (250x2 bps):__ \
					\n * Minimum overlapping length between R1 and R2 reads = %s; \
					\n * Minimum expected fragment length = %s;\
					\n * Maximum expected fragment length = %s;\
					\n * Maximum expected errors number = %s;\
					\n * Maximum number of differences = %s' %(minovlen_n, minmergelen_n, maxmergelen_n, maxee_n, maxdiffs_n)
			elif st.session_state.reads_lenght_R1_n_input == 300:
				if st.session_state.hypervar_regions_radio == 'V3V4':
					minovlen_n = 60
					minmergelen_n = 440
					maxmergelen_n = 480
					maxee_n = 1
					maxdiffs_n = 10
				elif st.session_state.hypervar_regions_radio == 'V4':
					minovlen_n = 220
					minmergelen_n = 230
					maxmergelen_n = 270 # average 254
					maxee_n = 1
					maxdiffs_n = 20
				
				suggested_joining_params_label = '\
					\n> __MiSeq V3 flowcells (300x2 bps):__ \
					\n * Minimum overlapping length between R1 and R2 reads  = %s; \
					\n * Minimum expected fragment length  = %s;\
					\n * Maximum expected fragment length = %s;\
					\n * Maximum expected errors number = %s;\
					\n * Maximum number of differences = %s' %(minovlen_n, minmergelen_n, maxmergelen_n, maxee_n, maxdiffs_n)
			else:
				if st.session_state.hypervar_regions_radio == 'V3V4':
					suggested_joining_params_label = 'Expected fragment length = 464 nucleotides.'
				elif st.session_state.hypervar_regions_radio == 'V4':
					suggested_joining_params_label = 'Expected fragment length = 254 nucleotides.'
			
			st.markdown(suggested_joining_params_label)
			

			if (('minovlen_n' in globals()) and ('minmergelen_n' in globals()) and ('maxmergelen_n' in globals()) and ('maxee_n' in globals()) and ('maxdiffs_n' in globals())):

				with st.form(key='reads_merging_parameters_form'):
					st.markdown('Set the parameters:')
					# cols = st.columns(5)
					minovlen = st.slider(label='Minimum overlapping length between R1 and R2 reads:',
						max_value=int(250),
						value=minovlen_n,
						step=1,
						help='%s'%(minovlen_n),
						key='minovlen_slider')
					minmergelen = st.slider(label='Minimum expected fragment length:',
						max_value=int(500),
						value=minmergelen_n,
						step=1,
						key='minmergelen_slider')
					maxmergelen = st.slider(label='Maximum expected fragment length:',
						max_value=int(500),
						value=maxmergelen_n,
						step=1,
						key='maxmergelen_slider')
					maxee = st.slider(label='Maximum expected errors number:',
						max_value=int(500),
						value=maxee_n,
						step=1,
						key='maxee_slider')
					maxdiffs = st.slider(label='Maximum number of differences:',
						max_value=int(500),
						value=maxdiffs_n,
						step=1,
						key='maxdiffs_slider')
					submit_button = st.form_submit_button(label='Merge paired-end reads')
			else:

				with st.form(key='reads_merging_parameters_form'):
					st.markdown('Set the parameters:')
					# cols = st.columns(5)
					minovlen = st.slider(label='Minimum overlapping length between R1 and R2 reads:',
						max_value=int(250),
						step=1,
						help='Minimum length of the overlapping region between R1 and R2 reads during merging.',
						key='minovlen_slider')
					minmergelen = st.slider(label='Minimum expected fragment length:',
						max_value=int(500),
						step=1,
						help='Minimum length of the merged fragment to be retained and not discarded.',
						key='minmergelen_slider')
					maxmergelen = st.slider(label='Maximum expected fragment length:',
						max_value=int(500),
						step=1,
						help='Maximum length of the merged fragment to be retained and not discarded.',
						key='maxmergelen_slider')
					maxee = st.slider(label='Maximum expected errors number:',
						max_value=int(500),
						step=1,
						help='Maximum number of expected errors of the merged fragment to be retained and not discarded.',
						key='maxee_slider')
					maxdiffs = st.slider(label='Maximum number of differences:',
						max_value=int(500),
						step=1,
						help='Maximum number of mismatches in the area of overlap during R1 and R2 reads merging.',
						key='maxdiffs_slider')
					submit_button = st.form_submit_button(label='Merge paired-end reads')

			if ((submit_button) or (st.session_state.minovlen_slider != 0)):

				df_demux_joined_summary, joined_sequences, pdf_joined_sequences, secure_temp_dir_joined_summary, stderr_joining = form_callback_vsearch_join(st.session_state.imported_sequences, st.session_state.imported_sequences_temp_dir, st.session_state.secure_temp_dir_demux_summary)
				st.session_state.joined_sequences = joined_sequences
				# with st.expander('Fusione delle letture R1 ed R2 in corso ...'):
				# 	st.write('prova std err stderr_joining') #stderr_joining)
				# 	st.write('%s'%(stderr_joining))

				side_placeh1.success('***%s.*** Set R1 and R2 reads merging parameters \
					\n > Tab %s. R1 and R2 merging' %(step_n, step_n))
				
				step_n += 1

				df = df_demux_joined_summary.rename({'forward sequence count': 'merged R1R2 sequence count'}, axis=1)

				st.subheader('Merged sequence counts per sample')
				st.write('Total Samples: %s (R1R2 merged)' %(
					len(~df['merged R1R2 sequence count'].isna())))
				df = df.merge(st.session_state.demux_summary_df.loc[:,'forward sequence count'], left_index=True, right_index=True)
				df.eval("""
				perc_starting_seqs_merged = ((`merged R1R2 sequence count`) * 100)/ (`forward sequence count`)
				""",
				inplace = True
				)
				df = df.drop(axis = 1, labels = ['forward sequence count'])
				st.table(df.style.format(formatter='{:,.0f}'))
				df = df.drop(axis = 1, labels= ['perc_starting_seqs_merged'])
				
				csv = convert_df(df)
				ste.download_button(label='Download CSV table',
					data=csv,
					file_name='merged_sequence_counts_per_sample.csv',
					mime='text/csv')

				st.subheader('Merged sequence frequency histogram')
				#st.markdown(pdf_joined_sequences, unsafe_allow_html=True)
				fig, ax = plt.subplots()
				ax = df['merged R1R2 sequence count'].plot.hist(bins=bins, color='grey', linewidth=1, edgecolor='black')
				ax.set_xlabel('Number of sequences')
				ax.set_ylabel('Number of samples')
				cols = st.columns((1,1)) # Serve per mantenere una dimensione ridotta del grafico sotto
				# che altrimenti occupa tutta la pagina. L'argomento figsize in plt.subplots() non funziona
				cols[0].pyplot(fig)
				
				st.subheader('Merged sequence counts summary')
				tabella_df = pd.concat([pd.DataFrame(df.min().rename('Minimum')).T,
					pd.DataFrame(df.median().rename('Median')).T,
					pd.DataFrame(df.mean().rename('Mean')).T,
					pd.DataFrame(df.max().rename('Maximum')).T,
					pd.DataFrame(df.sum().rename('Total')).T])

				st.table(tabella_df.style.format(formatter='{:,.2f}'))
				
				csv = convert_df(tabella_df)
				ste.download_button(label='Download CSV table',
					data=csv,
					file_name='merged_sequence_counts_summary.csv',
					mime='text/csv')

			else:
				st.info('Use the form above to set R1 and R2 reads merging parameters.')
				st.warning("The page is awaiting for you to set the desired parameters for R1 and R2 reads merging \
					using the form above and to click __Merge paired-end reads__.")
				st.stop()
				step_n += 1
			
			

	except Exception as e:

		#st.exception(e)
		
		with tab_summary_stats:
			st.warning('Single-end sequences, skip reads merging step.')
			# st.exception(e)
		
		# step_n += 1
		pass

	with tab_QC:
		
		side_plchldr2 = sidemenus.empty()
		side_plchldr2a = sidemenus.empty()
		side_plchldr2.info('***%s.*** Set quality control parameters \
			\n > Tab %s. Quality control' %(step_n, step_n))
		side_plchldr2a.markdown(f"<a href='#linkto_{step_n}_tab'>Tab {step_n}. Quality control</a>", unsafe_allow_html=True)
		st.markdown(f"<div id='linkto_{step_n}_tab'></div>", unsafe_allow_html=True)
		
		with st.form(key='quality_filter_form'):
			st.markdown('Set the parameters (pre-set defaults)')
			min_quality = st.number_input(
				label='Minimum acceptable Phred quality value, below which a base is considered \"low quality\":', 
				min_value=0, 
				max_value=40, 
				value=20, 
				step=1,
				help='Default value: 20.',
				key='min_quality_num') # minimum acceptable PHRED score
			quality_window = st.number_input(
				label='Maximum number of consecutive \"low quality\" bases observable, above which a sequence is truncated:', 
				min_value=0, 
				max_value=250, 
				value=4, 
				step=1,
				help='Default value: 4.',
				key='quality_window_num') # maximum number of low PHRED scores that can be observed in direct succession before truncating a sequence read.
			min_quality = int(min_quality)
			quality_window = int(quality_window)
			submit_button = st.form_submit_button('Quality control')
		
		if library_PE_SE == 'Paired-end':
			if ((submit_button) or ((st.session_state.min_quality_num != 0) and (st.session_state.quality_window_num != 0)) and ('joined_sequences' in st.session_state.keys())):
			
				try:
					
					demux_filter, df_q_filter, secure_temp_dir_q_filter_summary = quality_filter_paired_end(st.session_state.joined_sequences, min_quality, quality_window)
					st.session_state.demux_filter = demux_filter
					st.session_state.df_q_filter = df_q_filter

				except Exception as e:

					st.warning('Sequenze Single-end')
					# st.exception(e)
					demux_filter, df_q_filter, secure_temp_dir_q_filter_summary = quality_filter_paired_end(st.session_state.imported_sequences, min_quality, quality_window)
					st.session_state.demux_filter = demux_filter
					st.session_state.df_q_filter = df_q_filter
					
				
				side_plchldr2.success('***%s.*** Set quality control parameters \
				\n > Tab %s. Quality control' %(step_n, step_n))
				step_n += 1

				st.subheader('Quality control stats')
				st.write(df_q_filter.style.format(formatter='{:,.0f}'))
				
				csv = convert_df(df_q_filter)
				ste.download_button(label='Download CSV table',
					data=csv,
					file_name='QC_stats.csv',
					mime='text/csv')
				
				st.balloons()


			else:
				st.info('Use the form above to set the quality control parameters.')
				st.warning("The page is awaiting for you to set the desired quality control parameters \
					using the form above.")
			
				step_n += 1
		elif library_PE_SE == 'Single-end':

			if ((submit_button) or ((st.session_state.min_quality_num != 0) and (st.session_state.quality_window_num != 0))):
			
				try:
					
					demux_filter, df_q_filter, secure_temp_dir_q_filter_summary = quality_filter_paired_end(st.session_state.joined_sequences, min_quality, quality_window)
					st.session_state.demux_filter = demux_filter
					st.session_state.df_q_filter = df_q_filter

				except Exception as e:

					st.warning('Sequenze Single-end')
					# st.exception(e)
					demux_filter, df_q_filter, secure_temp_dir_q_filter_summary = quality_filter_paired_end(st.session_state.imported_sequences, min_quality, quality_window)
					st.session_state.demux_filter = demux_filter
					st.session_state.df_q_filter = df_q_filter
					
				
				side_plchldr2.success('***%s.*** Set quality control parameters \
				\n > Tab %s. Quality control' %(step_n, step_n))
				step_n += 1

				st.subheader('Quality control stats')
				st.write(df_q_filter.style.format(formatter='{:,.0f}'))
				
				csv = convert_df(df_q_filter)
				ste.download_button(label='Download CSV table',
					data=csv,
					file_name='QC_stats.csv',
					mime='text/csv')
				
				st.balloons()


			else:
				st.info('Use the form above to set the quality control parameters.')
				st.warning("The page is awaiting for you to set the desired quality control parameters \
					using the form above.")
			
				step_n += 1

	with tab_denoising:

		side_plchldr3 = sidemenus.empty()
		side_plchldr3a = sidemenus.empty()
		side_plchldr3.info('***%s.*** Denoising sequences \
			\n > Tab %s. Denoising' %(step_n, step_n))
		side_plchldr3a.markdown(f"<a href='#linkto_{step_n}_tab'>Tab {step_n}. Denoising</a>", unsafe_allow_html=True)
		st.markdown(f"<div id='linkto_{step_n}_tab'></div>", unsafe_allow_html=True)
		
		side_form_denoising_pipes_plchldr = sidemenus.empty()
		# Valori da aggiornare per aggiungere metodi di denoising nuovi alla applicazione 
		denoising_opts = ['Dada2', 'Deblur']
		denoising_pipes_funcs = [dada2_denoise_single_joined, deblur_denoise_trim_paired_end]
		denoising_pipes_funcs_d = {
			k:v for k,v in zip(denoising_opts, denoising_pipes_funcs)
			}

		with side_form_denoising_pipes_plchldr.form('side_form_denoising_pipes'):
			denoising_pipes = st.multiselect('Desired pipelines for the denoising of sequences:',
			options = denoising_opts,
			default = 'Deblur',
			key='denoising_pipes_multisel',
			help='Enter the desired denoising method/s to obtain the ASVs.\
				Available pipelines at the moment: Deblur, Dada2. It\'s possible to select all the two methods \
					to compare the results.')

			submit_button = st.form_submit_button('Denoising')
			
		if ((submit_button) or (('denoising_pipes_multisel' in st.session_state.keys()) and (st.session_state.denoising_pipes_multisel != []))):
		
			# trimming globale delle sequenze durante il denoising con deblur 
			trim_filtered_seqs = st.checkbox(label='Global trimming of all sequences during denoising?',
				value=True,
				key= 'trim_checkbox')
			if st.session_state.trim_checkbox:
				
				try:

					if st.session_state.hypervar_regions_radio == 'V3V4':
						value_N = 460 # average length is 464
					elif st.session_state.hypervar_regions_radio == 'V4':
						value_N = 250 # average length is 254
					
				except: # single-end sequences
					
					if 'hypervar_regions_radio' not in st.session_state.keys():
						hypervar_regions_form = hypervar_regions_plchldr.form('hypervar_regiorns_form')
					with hypervar_regions_form:
						hypervar_regions = st.radio(label='Hypervariable 16S rRNA gene regions target of sequencing:',
						options=['V3V4', 'V4'],
						index=0,
						help='Regions target of sequencing to be considered to determine expected fragment length.',
						key='hypervar_regions_radio')
						hypervar_btn = st.form_submit_button('Confirm')
					
					if ((hypervar_btn) or ('hypervar_regions_radio' in st.session_state.keys())):
						pass
					else:
						st.warning('The page is awaiting for you to select the region target of sequencing in the sidebar menu.')
						st.stop()
						
					st.markdown('Consider the length of the reads to set a global trimming value under the mean length:')
					reads_lenght_R1_input = st.number_input('Mean forward reads length:', help='To set the correct length \
						for global reads trimming during denoising. \
							\n> Pre-set mean value based on the previous analyses. \
							This value varies based on the sequencing technology. \
							Contact the raw data provider.',
							min_value=100, max_value=300, step=50, value=approx_reads_lenght_R1, key='reads_lenght_R1_n_input')

					if st.session_state.hypervar_regions_radio == 'V3V4':
						value_N = st.session_state.reads_lenght_R1_n_input # average length is 464
					elif st.session_state.hypervar_regions_radio == 'V4':
						value_N = st.session_state.reads_lenght_R1_n_input # average length is 254
					

				with st.form('global trim length'):

					# lunghezza del taglio basata sull abbassamento della qualita solitamente riscontrato nelle R2  verso la fine della lettura (3' estremita')
					N = st.number_input(
						label='Global trimming length:', 
						min_value=0, 
						max_value=500, 
						value=value_N,
						help='region/s %s: mean expected length %s'%(st.session_state.hypervar_regions_radio, value_N),
						step=1,
						key='N_global_trimming')
					st.session_state.N = int(st.session_state.N_global_trimming)
				
			
					submit_button = st.form_submit_button('Global trim reads')
			else:

				st.session_state.N = 0
			
			if (submit_button or ((st.session_state.trim_checkbox is True) and ('demux_filter' in st.session_state.keys()))):
				
				secure_temp_dir_filtered_seqs = tempfile.mkdtemp(prefix="temp_", suffix="_filtered_seqs")
				uuid = st.session_state.demux_filter.filtered_sequences.uuid
				st.session_state.demux_filter.filtered_sequences.save(secure_temp_dir_filtered_seqs+'/filtered_seqs.qza')
				st.session_state.demux_filter.filtered_sequences.extract(filepath = secure_temp_dir_filtered_seqs+'/filtered_seqs.qza', output_dir='%s'%secure_temp_dir_filtered_seqs)
				demux_filter_1 = import_SequencesWithQuality(_filepath=secure_temp_dir_filtered_seqs+'/%s/data/'%uuid)
				st.session_state.filtered_sequences = demux_filter_1
				
				# if len(st.session_state.denoising_pipes_multisel > 1):
					
				cols = st.columns(len(st.session_state.denoising_pipes_multisel))
				
				denoising_fs = [denoising_pipes_funcs_d[k] for k in st.session_state.denoising_pipes_multisel]
				for i, (denoising_f, denoising_pipe) in enumerate(zip(denoising_fs, st.session_state.denoising_pipes_multisel)):
					
					with cols[i]:
						
						st.header(denoising_pipe)
						try:
							denoised_sequences = denoising_f(
								_demux_filter=st.session_state.filtered_sequences, 
								N=st.session_state.N, 
								trim_TF=st.session_state.trim_checkbox
							)
						except Exception as e:
							print(e)
							st.warning('The %s denoising of merged reads throwed an exception.'%(denoising_pipe))
							st.stop()
						
						st.session_state['denoised_sequences_%s'%(denoising_pipe)] = denoised_sequences
						st.session_state['feature_table_%s'%(denoising_pipe)] = denoised_sequences.table
						st.session_state['rep_seqs_%s'%(denoising_pipe)] = denoised_sequences.representative_sequences
						
						if st.session_state.library_radio == 'Paired-end':

							st.success('Merged R1R2 sequences succesfully denoised into ASVs with %s pipeline.'%(denoising_pipe))
						elif st.session_state.library_radio == 'Single-end':

							st.success('R1 sequences succesfully denoised into ASVs with %s pipeline.'%(denoising_pipe))
						
						# Tabella delle statistiche del processo di Deblurring
						secure_temp_dir_denoising = tempfile.mkdtemp(prefix="temp_", suffix="_denoising")
						try:
							denoising_stats = deblur.visualizers.visualize_stats(st.session_state['denoised_sequences_%s'%(denoising_pipe)].stats)
							denoising_stats.visualization.export_data(secure_temp_dir_denoising)
							with open(secure_temp_dir_denoising+"/index.html", 'r', encoding='utf-8') as HtmlFile:
								source_code = HtmlFile.read()
						
							tabella_df= pd.DataFrame(pd.read_html(source_code,thousands='.')[0]).rename({'Unnamed: 0': 'Index'},
								# ,
								# 'sample-id': 'ID campione', 
								# 'reads-raw': 'num totale letture iniziali',
								# 'fraction-artifact-with-minsize': 'frazione di letture artefatte incluse di lunghezza minima',
								# 'fraction-artifact': 'frazione di sequenze artefatte',
								# 'fraction-missed-reference': 'frazione di letture mancante di sequenze di riferimento positivo',
								# 'unique-reads-derep': 'letture uniche dopo dereplicazione',
								# 'reads-derep': 'num letture dopo dereplicazione (escluso singoletti)',
								# 'unique-reads-deblur': 'letture uniche dopo il Deblur',
								# 'reads-deblur': 'num di letture dopo il Deblur',
								# 'unique-reads-hit-artifact': 'letture uniche che combaciano con sequenze artefatte negative',
								# 'reads-hit-artifact': 'num di letture che combaciano con sequenze artefatte negative',
								# 'unique-reads-chimeric': 'letture uniche chimeriche',
								# 'reads-chimeric': 'num di letture chimeriche',
								# 'unique-reads-hit-reference': 'letture uniche di Deblur che combaciano con sequenze di riferimento positivo',
								# 'reads-hit-reference': 'num di letture di Deblur che combaciano con sequenze di riferimento positivo',
								# 'unique-reads-missed-reference': 'letture uniche di Deblur mancanti di sequenza di riferimento positivo',
								# 'reads-missed-reference': 'letture di Deblur mancanti di sequenza di riferimento positivo'},
								axis=1).drop('Index', axis=1)
							tabella_df = tabella_df.set_index('sample-id')
							st.session_state['stats_%s_per_sample'%(denoising_pipe)] = tabella_df
							st.subheader('%s denoising stats per sample'%(denoising_pipe))
							st.session_state['stats_denoising_%s_per_sample'%(denoising_pipe)] = tabella_df
							st.dataframe(tabella_df.style.format(formatter='{:,.0f}'))
							
							csv = convert_df(tabella_df)
							ste.download_button(label='Download CSV table',
								data=csv,
								file_name='%s_denoising_stats_per_sample.csv'%(denoising_pipe),
								mime='text/csv')
						except: # dada2
							denoising_stats = st.session_state['denoised_sequences_%s'%(denoising_pipe)].denoising_stats
							denoising_stats.export_data(secure_temp_dir_denoising)
							tabella_df = pd.read_csv(secure_temp_dir_denoising+"/stats.tsv", sep='\t', skiprows=[1], index_col=0)
							st.subheader('%s denoising stats per sample'%(denoising_pipe))
							tabella_df = tabella_df
							# .rename({
							# 	'input': 'letture iniziali',
							# 	'filtered': 'filtrate',
							# 	'percentage of input passed filter': 'percentuale letture iniziali post-filtro',
							# 	'denoised': 'post-denoising',
							# 	'non-chimeric': 'non-chimeriche',
							# 	'percentage of input non-chimeric': 'percentuale letture iniziali non-chimeriche'}, axis = 1)
							st.session_state['stats_denoising_%s_per_sample'%(denoising_pipe)] = tabella_df
							
							st.dataframe(tabella_df.style.format(formatter='{:,.0f}'))
							
							csv = convert_df(tabella_df)
							ste.download_button(label='Download CSV table',
								data=csv,
								file_name='%s_denoising_stats_per_sample.csv'%(denoising_pipe),
								mime='text/csv')
						
						# Tabella delle Features ASVs
						st.subheader('ASVs table - Amplicon Sequence Variants (features)')
						try: # Deblur
							df_feature_table = st.session_state['feature_table_%s'%(denoising_pipe)].view(pd.DataFrame).T
							df_feature_table.index.name = '#NAME'
							st.write(df_feature_table.style.format(formatter='{:,.0f}'))
							
							csv = convert_df(df_feature_table)
							ste.download_button(label='Download CSV table',
								data=csv,
								file_name='ASVs_feature_table_%s.csv'%(denoising_pipe),
								mime='text/csv')
						except Exception as e:
							st.exception(e)
							st.warning('Dada2 denoising exception, feature table visualization as table, look for comment # Tabella delle Features ASVs')
							

						# sequenze ASVs rappresentative
						st.subheader('ASVs representative sequences table')
						df_ASVs_rep_seqs = st.session_state['rep_seqs_%s'%(denoising_pipe)].view(pd.Series) 
						df_ASVs_rep_seqs.name = 'Sequence'
						# st.write(df_ASVs_rep_seqs) # Non si visualizza una semplice serie ma un oggetto di tipo DNA con metadata e statistiche!

						csv = convert_df(df_ASVs_rep_seqs)
						ste.download_button(label='Download CSV table',
							data=csv,
							file_name='ASVs_representative_sequences_%s.csv'%(denoising_pipe),
							mime='text/csv')

						# Tabella delle statistiche della Feature Table
						st.header('ASVs summary statistics')
						st.subheader('ASVs table summary')
						try: # Deblur
							output_viz = feature_table.visualizers.summarize(st.session_state['feature_table_%s'%(denoising_pipe)])
							output_viz.visualization.export_data(secure_temp_dir_denoising)
							st.success("Feature table summary statistics computed.")
							with open(secure_temp_dir_denoising+"/index.html", 'r', encoding='utf-8') as HtmlFile:
								source_code = HtmlFile.read()
							tabella_df = pd.DataFrame(pd.read_html(source_code,thousands=',', decimal='.')[0])
							tabella_df = tabella_df.set_index('Metric', drop = True)
							# tabella_df = tabella_df.rename({'Number of samples': 'Numero di campioni',
							# 'Number of features': 'Numero di OTUs o ASVs',
							# 'Total frequency': 'Frequeza totale'}, axis=0).rename({'Sample': 'Campione'}, axis = 1)
							st.session_state['stats_%s_per_sample'%(denoising_pipe)] = tabella_df
							st.subheader('%s stats per sample'%(denoising_pipe))
							st.table(tabella_df.style.format(formatter='{:,.0f}'))

							csv = convert_df(tabella_df)
							ste.download_button(label='Download CSV table',
								data=csv,
								file_name='ASVs_feature_table_%s_summary.csv'%(denoising_pipe),
								mime='text/csv')
						except: # dada2
							output_viz = feature_table.visualizers.summarize(st.session_state['denoised_sequences_%s'%(denoising_pipe)].table)
							output_viz.visualization.export_data(secure_temp_dir_denoising)
							st.success("Feature table summary statistics computed.")
							with open(secure_temp_dir_denoising+"/index.html", 'r', encoding='utf-8') as HtmlFile:
								source_code = HtmlFile.read()
							tabella_df= pd.DataFrame(pd.read_html(source_code,thousands='.')[0])
							# tabella_df = tabella_df.rename({'Number of samples': 'Numero di campioni',
							# 'Number of features': 'Numero di OTUs o ASVs',
							# 'Total frequency': 'Frequeza totale'}, axis=0).rename({'Sample': 'Campione'}, axis = 1)
							st.session_state['stats_%s_per_sample'%(denoising_pipe)] = tabella_df
							st.subheader('%s stats per sample'%(denoising_pipe))
							st.table(tabella_df.style.format(formatter='{:,.0f}'))
							
							csv = convert_df(tabella_df)
							ste.download_button(label='Download CSV table',
								data=csv,
								file_name='ASVs_feature_table_%s_summary.csv'%(denoising_pipe),
								mime='text/csv')
						
						st.subheader('ASVs frequency per sample summary')
						tabella_df = pd.DataFrame(pd.read_html(source_code,thousands=',', decimal='.')[1])
						tabella_df = tabella_df.set_index('Unnamed: 0', drop = True)
						# tabella_df = tabella_df.rename({'Frequency': 'Frequenza'}, axis = 1).rename({'Minimum frequency': 'Frequenza minima',
						# '1st quartile': '1o quartile',
						# 'Median frequency': 'Frequenza mediana',
						# '3rd quartile': '3o quartile',
						# 'Maximum frequency': 'Frequenza massima',
						# 'Mean frequency': 'Frequenza media'}, axis = 0)
						st.table(tabella_df.style.format(formatter='{:,.0f}'))
						
						csv = convert_df(tabella_df)
						ste.download_button(label='Download CSV table',
							data=csv,
							file_name='ASVs_freqs_per_sample_%s_summary.csv'%(denoising_pipe),
							mime='text/csv')

						
						
						st.subheader('Samples frequency per ASV summary')
						tabella_df = pd.DataFrame(pd.read_html(source_code,thousands=',', decimal='.')[2])
						tabella_df = tabella_df.set_index('Unnamed: 0', drop = True)
						# tabella_df = tabella_df.rename({'Frequency': 'Frequenza'}, axis = 1).rename({'Minimum frequency': 'Frequenza minima',
						# '1st quartile': '1o quartile',
						# 'Median frequency': 'Frequenza mediana',
						# '3rd quartile': '3o quartile',
						# 'Maximum frequency': 'Frequenza massima',
						# 'Mean frequency': 'Frequenza media'}, axis = 0)
						st.table(tabella_df.style.format(formatter='{:,.0f}'))
						
						csv = convert_df(tabella_df)
						ste.download_button(label='Download CSV table',
							data=csv,
							file_name='sample_freqs_per_ASV_%s_summary.csv'%(denoising_pipe),
							mime='text/csv')
						
						
						# Generate deblur representative sequences summary visualizations
						st.header('Representative ASVs sequences stats')
						st.subheader('Sequences summary')
						deblur_feature_table_summary = feature_table.visualizers.tabulate_seqs(st.session_state['rep_seqs_%s'%(denoising_pipe)])
						deblur_feature_table_summary.visualization.export_data(secure_temp_dir_denoising)
						with open(secure_temp_dir_denoising+"/index.html", 'r', encoding='utf-8') as HtmlFile:
							source_code = HtmlFile.read() 
						tabella_df = pd.DataFrame(pd.read_html(source_code,thousands=',')[0]).T
						tabella_df = tabella_df.rename({0: 'Value'}, axis=1)
						# .rename({
						# 	'Sequence Count': 'Conta sequenze',
						# 	'Min Length': 'Lunghezza minima', 
						# 	'Max Length': 'Lunghezza massima',
						# 	'Mean Length': 'Lunghezza media',
						# 	'Range': 'Intervallo',
						# 	'Standard Deviation': 'Deviazione standard'}, axis = 0)
						st.table(tabella_df.style.format(formatter='{:,.0f}'))
						
						csv = convert_df(tabella_df)
						ste.download_button(label='Download CSV table',
							data=csv,
							file_name='ASVs_representative_sequences_%s_summary.csv'%(denoising_pipe),
							mime='text/csv')

						st.subheader('Sequence length distribution')
						tabella_df = pd.DataFrame(pd.read_html(source_code,thousands=',', header=0)[1]).T
						tabella_df_hdr = tabella_df.iloc[0,:]
						tabella_df = tabella_df.iloc[1:,:]
						tabella_df.columns = tabella_df_hdr
						tabella_df.index.name = 'Percentile'
						tabella_df = tabella_df#.rename({'Length* (nts):': 'Lunghezza in nucleotidi'}, axis=1)
						st.table(tabella_df)#.style.format(formatter='{:,.0f}'))
						
						csv = convert_df(tabella_df)
						ste.download_button(label='Download CSV table',
							data=csv,
							file_name='ASVs_representative_sequences_%s_length_distribution.csv'%(denoising_pipe),
							mime='text/csv')
						
						st.subheader('Table of sequences lengths and nucleotides with links to ncbi BLAST')
						st.markdown('For ncbi BLAST see sequences [].')
						tabella_df = pd.DataFrame(pd.read_html(source_code,thousands=',', header=0)[2])
						tabella_df = tabella_df#.rename({'Feature ID': 'ID OTUs o ASVs', 'Sequence Length': 'Lunghezza sequenza', 'Sequence': 'Sequenza'}, axis = 1)
						
						sequenza_link_BLAST = ['[%s](http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?ALIGNMENT_VIEW=Pairwise&PROGRAM=blastn&DATABASE=nt&CMD=Put&QUERY=%s)'%(seq, seq) for seq in tabella_df['Sequence'].values]
						
						tabella_df['Sequence'] = sequenza_link_BLAST
						st.dataframe(tabella_df)
						
						csv = convert_df(tabella_df)
						ste.download_button(label='Download CSV table',
							data=csv,
							file_name='ASVs_%s_lengths_and_nt_seqs_table.csv'%(denoising_pipe),
							mime='text/csv')
					

			side_plchldr3.success('***%s.*** Denoising \
			\n > Tab %s. Denoising' %(step_n, step_n))
			step_n += 1

			st.balloons()
		else:

			st.warning('The page is awaiting for you to enter the denosing method in the sidebar menu form.')
			step_n += 1
			st.stop()

	
	with tab_taxonomy:

		side_plchldr4 = sidemenus.empty()
		side_plchldr4a = sidemenus.empty()
		side_plchldr4.info('***%s.*** Taxonomic classification \
			\n > Tab %s. Taxonomic classification' %(step_n, step_n))
		side_plchldr4a.markdown(f"<a href='#linkto_{step_n}_tab'>Tab {step_n}. Taxonomic classification</a>", unsafe_allow_html=True)
		st.markdown(f"<div id='linkto_{step_n}_tab'></div>", unsafe_allow_html=True)
		

		# form_upload_classifier_data_plchldr = st.empty()
		# with form_upload_classifier_data_plchldr.form('upload_classifier_data'):
			
		# 	path_pre_trained_classifier = st.file_uploader('Select a pre-trained classifier file:', 
		# 	key='pre_trained_classifier_path_input',
		# 	accept_multiple_files = False,
		# 	type='qza',
		# 	#value='/app/pre_trained_classifier/gg-13-8-99-nb-weighted-classifier.qza',
		# 	help='Pre-trained classifier file, available at https://docs.qiime2.org/2023.2/data-resources/, file format: .qza (artifact qiime2). \
		# 		Default: pre_trained_classifier/gg-13-8-99-nb-weighted-classifier.qza')

		# 	path_reference_taxonomy = st.file_uploader('Select a reference taxonomy file:',
		# 	key='reference_taxonomy_path_input',
		# 	accept_multiple_files = False,
		# 	type='txt',
		# 	#value='/app/microbiome/reference_seqs/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt',
		# 	help='Reference taxonomy file of OTU sequences of the pre-trained classifier, formato del file: .txt. \
		# 		Default: /app/microbiome/reference_seqs/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt')

		# 	path_reference_otus_seqs = st.file_uploader('Select a reference OTU file:',
		# 	key='reference_otus_seqs_path_input',
		# 	accept_multiple_files = False,
		# 	type='fasta',
		# 	#value='app/microbiome/reference_seqs/gg_13_8_otus/rep_set/99_otus.fasta',
		# 	help='OTU reference sequences file of the pre-trained classifier, file format: .fasta. \
		# 		Default: /app/microbiome/reference_seqs/gg_13_8_otus/rep_set/99_otus.fasta')

		# 	submit_button = st.form_submit_button('Upload')

		# if ((submit_button) or (
		# 	('pre_trained_classifier_path_input' in st.session_state.keys()) and
		# 	('reference_taxonomy_path_input' in st.session_state.keys()) and
		# 	('reference_otus_seqs_path_input' in st.session_state.keys()) and
		# 	(st.session_state.pre_trained_classifier_path_input is not None) and
		# 	(st.session_state.reference_taxonomy_path_input is not None) and
		# 	(st.session_state.reference_otus_seqs_path_input is not None))):
			
		# 	with NamedTemporaryFile(dir='.', suffix='.qza') as f:
		# 		f.write(st.session_state.pre_trained_classifier_path_input.getbuffer())
		# 		classifier = Artifact.load(f.name)
		# 	with NamedTemporaryFile(dir='.', suffix='.txt') as f:
		# 		f.write(st.session_state.reference_taxonomy_path_input.getbuffer())
		# 		ref_taxonomy = import_ref_gg_13_8_otus_taxonomy(f.name)
		# 	with NamedTemporaryFile(dir='.', suffix='.fasta') as f:
		# 		f.write(st.session_state.reference_otus_seqs_path_input.getbuffer())
		# 		ref_otus_seqs = import_ref_gg_13_8_otus_seqs(f.name)
		st.session_state.pre_trained_classifier_path_input = '/app/microbiome/pre_trained_classifier/gg-13-8-99-nb-weighted-classifier.qza'
		st.session_state.reference_taxonomy_path_input = '/app/microbiome/reference_seqs/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt'
		st.session_state.reference_otus_seqs_path_input = '/app/microbiome/reference_seqs/gg_13_8_otus/rep_set/99_otus.fasta'

		classifier = Artifact.load(st.session_state.pre_trained_classifier_path_input)
		ref_taxonomy = import_ref_gg_13_8_otus_taxonomy(st.session_state.reference_taxonomy_path_input)
		ref_otus_seqs = import_ref_gg_13_8_otus_seqs(st.session_state.reference_otus_seqs_path_input)
		
		st.session_state.classifier = classifier
		st.session_state.ref_taxonomy = ref_taxonomy
		st.session_state.ref_otus_seqs = ref_otus_seqs
		
		cols = st.columns(len(st.session_state.denoising_pipes_multisel))

		denoising_fs = [denoising_pipes_funcs_d[k] for k in st.session_state.denoising_pipes_multisel]
		for i, (denoising_f, denoising_pipe) in enumerate(zip(denoising_fs, st.session_state.denoising_pipes_multisel)):
			
			with cols[i]:
				
				st.header(denoising_pipe)
				
				try:

					seqs_classification = app_classify_hybrid_vsearch_sklearn(
						_query=st.session_state['rep_seqs_%s'%(denoising_pipe)],
						_reference_reads=st.session_state.ref_otus_seqs,
						_reference_taxonomy=st.session_state.ref_taxonomy,
						_classifier=st.session_state.classifier,
						reads_per_batch=1000,
						randseed=3)
					
					st.session_state['classification_%s'%(denoising_pipe)] = seqs_classification.classification

				except Exception as e:
					
					st.exception(e)
					st.warning('Sequences classification throwed an exception.')
					st.stop()
				
				st.success("Successfully classified sequences.")
				st.balloons()

				st.header('Taxonomic classification of ASVs')
				tabella_df = st.session_state['classification_%s'%(denoising_pipe)].view(pd.DataFrame)['Taxon'].str.split(';', expand=True).rename({0: 'Kingdom',
				1: 'Phylum', 2: 'Class', 3: 'Order', 4: 'Family', 5: 'Genus', 6: 'Species'}, axis=1)
				tabella_df.index.name = '#TAXONOMY'
				#filter_dataframe(tabella_df)
				st.write(tabella_df)

				csv = convert_df(tabella_df)
				ste.download_button(label='Download CSV table',
					data=csv,
					file_name='taxonomic_classification_ASVs_%s.csv'%(denoising_pipe),
					mime='text/csv')



		side_plchldr4.success('***%s.*** Taxonomic classification \
			\n > Tab %s. Taxonomic classification' %(step_n, step_n))
		step_n += 1
		
		# else:
		# 	st.info('Use the form above to select or drag and drop the file of the pre-trained classifier, available at https://docs.qiime2.org/2023.2/data-resources/.')
		# 	st.warning('The page is awaiting for you to upload the files for the taxonomic classifier.')
		# 	# st.stop()
		# 	step_n += 1


	with tab_rarefaction:

		side_plchldr5 = sidemenus.empty()
		side_plchldr5a = sidemenus.empty()
		side_plchldr5.info('***%s.*** Normalization \
			\n > Tab %s. Normalization' %(step_n, step_n))
		side_plchldr5a.markdown(f"<a href='#linkto_{step_n}_tab'>Tab {step_n}. Normalization</a>", unsafe_allow_html=True)
		st.markdown(f"<div id='linkto_{step_n}_tab'></div>", unsafe_allow_html=True)
		
		metrics = {'observed_features', 'shannon', 'pielou_e', 'simpson'}
		
		cols = st.columns(len(st.session_state.denoising_pipes_multisel))

		for i, (denoising_f, denoising_pipe) in enumerate(zip(denoising_fs, st.session_state.denoising_pipes_multisel)):
			
			try: # Deblur
				max_depth_suggested = myround(np.mean(st.session_state['stats_denoising_%s_per_sample'%(denoising_pipe)]['reads-hit-reference']))
			except: # Dada2
				max_depth_suggested = myround(np.mean(st.session_state['stats_denoising_%s_per_sample'%(denoising_pipe)].loc[:,'non-chimeric']))
				
			st.session_state['max_depth_suggested_%s'%(denoising_pipe)] = max_depth_suggested
			max_depth = st.session_state['max_depth_suggested_%s'%(denoising_pipe)]
			with cols[i]:
				
				st.header(denoising_pipe)
				st.info('Default: maximum depth suggested for alpha rarefaction curves, equal to the approximate mean number of sequences per sample: %s'%(st.session_state['max_depth_suggested_%s'%(denoising_pipe)]))
				try:
					alpha_rare_curves = app_alpha_rare_curves(
						_table = st.session_state['feature_table_%s'%(denoising_pipe)], 
						max_depth = max_depth, 
						metrics = metrics)
					
				except Exception as e:
					st.exception(e)
					st.stop()

				st.session_state['alpha_rare_curves_%s'%(denoising_pipe)] = alpha_rare_curves
		
				secure_temp_dir_rarefaction = tempfile.mkdtemp(prefix="temp_", suffix="_alpha_rarefaction")
				alpha_rare_curves.visualization.export_data(secure_temp_dir_rarefaction)
				alpha_rare_curves.visualization.save(secure_temp_dir_rarefaction+'/alpha_rare_curves.qzv')

				st.subheader('Alpha rarefaction curves')
				st.markdown('You can download a .qzv file to visualize the alpha rarefaction curves \
					 at https://view.qiime2.org, to evaluate the preferred rarefaction depth \
						for the ASVs table, based on the dispersion of the computed alpha diversity metrics value \
							and on the number of excluded samples (you would want to retain as much samples and sequences as possible).')
				
				# aggiungere i file .jsonp .csv e .html e le cartelle q2templateassets e dist al bottone di download (.zip file), altrimenti index.html da solo non si visualizza correttamente nel browser
				# zipfolder(secure_temp_dir_rarefaction+"/zip_alpha_rare_curves.zip", secure_temp_dir_rarefaction)

				with open(secure_temp_dir_rarefaction+"/alpha_rare_curves.qzv", 'rb') as f:
					ste.download_button(
						label="Download alpha rarefaction curves .qzv",
						data=f,
						file_name="alpha_rarefaction_curves_%s.qzv" %(denoising_pipe),
						mime="application/qzv")
				

				with st.form('rarefaction_form_%s'%(denoising_pipe)):

					try: # Deblur
						sampling_depth_suggested = np.min(st.session_state['stats_denoising_%s_per_sample'%(denoising_pipe)]['reads-hit-reference'])
					except: # Dada2
						sampling_depth_suggested = np.min(st.session_state['stats_denoising_%s_per_sample'%(denoising_pipe)]['non-chimeric'])

					st.info('Default: Rarefaction depth pre-set, equals to the minimum number of sequences per sample: %s'%sampling_depth_suggested)
					st.info('ATTENTION: See alpha rarefaction curves to pick a sequencing depth that \
						allows to keep as much samples as possible and in the mean time as much sequences as possible (suggestion to exclude \
							those samples with very low number of sequences compared to the others, if present, and setting \
								a greater sequencing depth with respect to default, to retain more information of sequences).')
					sampling_depth = st.number_input(label='Choose the rarefaction depth:',
						min_value=0, 
						max_value=None, 
						value=sampling_depth_suggested, 
						step=1,
						help='Number of sequences for random subsampling of each sample. Consider alpha rarefaction curves.',
						key = 'sampling_depth_%s_input'%(denoising_pipe))
					st.session_state['sampling_depth_%s'%(denoising_pipe)] = int(sampling_depth)

					submit_button = st.form_submit_button('Rarefaction')
				
				if ((submit_button) or (('sampling_depth_%s_input'%(denoising_pipe) in st.session_state.keys()) and
				(st.session_state['sampling_depth_%s_input'%(denoising_pipe)] != 0))):
					
					try:

						rarefy_result = feature_table.methods.rarefy(table=st.session_state['feature_table_%s'%(denoising_pipe)], sampling_depth=st.session_state['sampling_depth_%s'%(denoising_pipe)])
						rarefied_table = rarefy_result.rarefied_table
						st.session_state['rarefied_table_%s'%(denoising_pipe)] = rarefied_table
					except:

						st.warning('Rarefaction throwed an exception.')
						st.stop()

				else:
					st.info('To choose the sampling depth in the form above, consider the alpha rarefaction curves.')
					st.warning('The page is awaiting for you to select a sampling depth using the form above.')

				st.balloons()

				
				
				st.success('Rarefaction at depth of %s sequences per sample successfully completed.' %(sampling_depth))
				st.subheader('ASVs table normalized with rarefaction method')
				tabella_df = st.session_state['rarefied_table_%s'%(denoising_pipe)].view(pd.DataFrame).T
				tabella_df.index.name = '#NAME'
				st.write(tabella_df.style.format(formatter='{:,.0f}'))

				csv = convert_df(tabella_df)
				ste.download_button(label='Download CSV table',
					data=csv,
					file_name='ASVs_feature_table_normalized_%s.csv'%(denoising_pipe),
					mime='text/csv')

				
		side_plchldr5.success('***%s.*** Normalization \
			\n > Tab %s. Normalization' %(step_n, step_n))
		step_n += 1


	with tab_diversity:

		side_plchldr6 = sidemenus.empty()
		side_plchldr6a = sidemenus.empty()
		side_plchldr6.info('***%s.*** Alpha and beta diversity \
			\n > Tab %s. Alpha and beta diversity' %(step_n, step_n))
		side_plchldr6a.markdown(f"<a href='#linkto_{step_n}_tab'>Tab {step_n}. Alpha and beta diversity</a>", unsafe_allow_html=True)
		st.markdown(f"<div id='linkto_{step_n}_tab'></div>", unsafe_allow_html=True)
		

		
		
		st.header('Alpha Diversity')
		if len(st.session_state.denoising_pipes_multisel) == 1:

			subhdr_plchldr = st.empty()
			cols_alpha_divs = st.columns(len(metrics)) # metrics definite prima in tab_rarefaction
		else:
			
			subhdr_plchldr = st.empty()
			cols = st.columns(len(st.session_state.denoising_pipes_multisel))			

		for j, denoising_pipe in enumerate(st.session_state.denoising_pipes_multisel):
			
			try:

				subhdr_plchldr.subheader(denoising_pipe)
				alpha_divs = app_alpha_divs(
					_table=st.session_state['rarefied_table_%s'%(denoising_pipe)])
				st.session_state.alpha_divs_df = pd.concat([i.view(pd.Series) for i in alpha_divs], axis=1)
				st.session_state.alpha_divs = alpha_divs
				for index, i in enumerate(alpha_divs):
					
					with cols_alpha_divs[index]:

						tabella_df = i.view(pd.Series)
						st.write(tabella_df)
						
						csv = convert_df(tabella_df)
						ste.download_button(label='Download CSV table',
							data=csv,
							file_name='alpha_diversity_%s_%s.csv'%(i.view(pd.Series).name, denoising_pipe),
							mime='text/csv')


			except Exception as e:

				with cols[j]:
					
					st.subheader(denoising_pipe)
					alpha_divs = app_alpha_divs(_table=st.session_state['rarefied_table_%s'%(denoising_pipe)])
					st.session_state.alpha_divs = alpha_divs
					st.session_state.alpha_divs_df = pd.concat([i.view(pd.Series) for i in alpha_divs], axis=1)
					for index, i in enumerate(alpha_divs):
						
						tabella_df = i.view(pd.Series)
						st.write(tabella_df)
						
						csv = convert_df(tabella_df)
						ste.download_button(label='Download CSV table',
							data=csv,
							file_name='alpha_diversity_%s_%s.csv'%(i.view(pd.Series).name, denoising_pipe),
							mime='text/csv')
			
		st.header('Beta Diversity')
		cols = st.columns(len(st.session_state.denoising_pipes_multisel))			

		for j, denoising_pipe in enumerate(st.session_state.denoising_pipes_multisel):
			
			try:

				st.subheader(denoising_pipe)
				st.markdown('***Bray Curtis***')
				st.session_state.beta_meta_df = pd.DataFrame(st.session_state.alpha_divs_df)
				st.session_state.beta_meta_df.index.name = '#SampleID'
				beta_bc_rare_vis = diversity.visualizers.beta_rarefaction(
					table = st.session_state['feature_table_%s'%(denoising_pipe)], 
					sampling_depth = st.session_state['sampling_depth_%s'%(denoising_pipe)], 
					metric = 'braycurtis',
					metadata = Metadata(st.session_state.beta_meta_df),
					clustering_method ='upgma',
					correlation_method='pearson'
					)
				
				secure_temp_dir_rarefaction_beta_bc = tempfile.mkdtemp(prefix="temp_", suffix="_beta_bc_rarefaction")

				beta_bc_rare_vis.visualization.export_data(secure_temp_dir_rarefaction_beta_bc)
				# Non funziona lo zipping della cartella beta rarefaction secure_temp_dir_rarefaction_beta!!
				#todo workaround temporaneo:
				st.markdown('You can download a .zip file with visualizations of results of beta rarefaction \
					at the depth of %s sequences per sample. \
					The file index.html shows: an Emperor jackknifed PCoA plot, \
						a Heatmap of Pearson correlation among rarefactions for beta diversity metric \
							***Bray-Curtis*** and a Newick tree of UPGMA clustering of samples. Go to \
									http://etetoolkit.org/treeview/ to visualize it.' %(st.session_state['sampling_depth_%s'%(denoising_pipe)]))
				zipfolder(secure_temp_dir_rarefaction_beta_bc+"/zip_beta_rare.zip", secure_temp_dir_rarefaction_beta_bc)

				with open(secure_temp_dir_rarefaction_beta_bc+"/zip_beta_rare.zip", 'rb') as f:
					
					ste.download_button(
						label="Download beta rarefaction Bray Curtis .zip",
						data=f,
						file_name="beta_rarefaction_braycurtis_%s.zip" %(denoising_pipe),
						mime="application/zip")
				
				st.markdown('***Jaccard***')
				beta_j_rare_vis = diversity.visualizers.beta_rarefaction(
					table = st.session_state['feature_table_%s'%(denoising_pipe)], 
					sampling_depth = st.session_state['sampling_depth_%s'%(denoising_pipe)], 
					metric = 'jaccard',
					metadata = Metadata(st.session_state.beta_meta_df),
					clustering_method ='upgma',
					correlation_method='pearson'
					)
				
				secure_temp_dir_rarefaction_beta_j = tempfile.mkdtemp(prefix="temp_", suffix="_beta_j_rarefaction")

				beta_j_rare_vis.visualization.export_data(secure_temp_dir_rarefaction_beta_j)
				# Non funziona lo zipping della cartella beta rarefaction secure_temp_dir_rarefaction_beta!!
				#todo workaround temporaneo:
				st.markdown('You can download a .zip file with visualizations of results of beta rarefaction \
					at the depth of %s sequences per sample. \
					The file index.html shows: an Emperor jackknifed PCoA plot, \
						a Heatmap of Pearson correlation among rarefactions for beta diversity metric \
							***Jaccard*** and a Newick tree of UPGMA clustering of samples. Go to \
									http://etetoolkit.org/treeview/ to visualize it.' %(st.session_state['sampling_depth_%s'%(denoising_pipe)]))
				
				zipfolder(secure_temp_dir_rarefaction_beta_j+"/zip_beta_rare.zip", secure_temp_dir_rarefaction_beta_j)

				with open(secure_temp_dir_rarefaction_beta_j+"/zip_beta_rare.zip", 'rb') as f:
					
					ste.download_button(
						label="Download beta rarefaction Jaccard .zip",
						data=f,
						file_name="beta_rarefaction_jaccard_%s.zip" %(denoising_pipe),
						mime="application/zip")
			except Exception as e:

				st.exception(e)	
				st.stop()
			
		st.balloons()

		side_plchldr6.success('***%s.*** Alpha and beta diversity \
				\n > Tab %s. Alpha and beta diversity' %(step_n, step_n))
		step_n += 1


	with tab_phylogenetic_tree:

		side_plchldr7 = sidemenus.empty()
		side_plchldr7a = sidemenus.empty()
		side_plchldr7.info('***%s.*** Phylogenetic tree construction \
			\n > Tab %s. Phylogenetic tree' %(step_n, step_n))
		side_plchldr7a.markdown(f"<a href='#linkto_{step_n}_tab'>Tab {step_n}. Phylogenetic tree</a>", unsafe_allow_html=True)
		st.markdown(f"<div id='linkto_{step_n}_tab'></div>", unsafe_allow_html=True)
		
		st.subheader('Core alpha and beta diversity metrics, based on the phylogenetic tree')
		st.markdown('_De novo_ construction of a phylogenetic tree from sequences, through multiple sequence alignment using MAFFT and fasttree.')
		cols = st.columns(len(st.session_state.denoising_pipes_multisel))			
		
		metadata_uploader = st.file_uploader('Upload a metadata file .tsv', accept_multiple_files = False, key = 'metadata_input')
		if st.session_state.metadata_input is not None:
			df_meta = pd.read_csv(st.session_state.metadata_input, sep='\t', index_col=0)
			df_meta.index.name = 'SampleID'
			st.session_state.metadata = Metadata(df_meta)
		elif st.button('Load sample metadata'):
			df_meta = pd.read_csv('/app/microbiome/sample_data/metadata_ASVs.txt', sep='\t', index_col=0)
			df_meta.index.name = 'SampleID'
			st.session_state.metadata = Metadata(df_meta)
		else:
			st.warning('The page is awaiting for you to upload a metadata file to be used for emperor plots of core beta diversity metrics.')
			st.stop()
		for j, denoising_pipe in enumerate(st.session_state.denoising_pipes_multisel):
			
			try:

				st.session_state.seqs_alignment, st.session_state.masked_alignment, st.session_state.tree, st.session_state.rooted_tree, st.session_state.core_metr_phylo = app_align_to_tree_mafft_fasttree(
					_sequences = st.session_state['rep_seqs_%s'%(denoising_pipe)],
					_table = st.session_state['feature_table_%s'%(denoising_pipe)],
					sampling_depth = st.session_state['sampling_depth_%s'%(denoising_pipe)],
					_metadata = st.session_state.metadata)
			
				secure_temp_dir_phylogenetic_tree = tempfile.mkdtemp(prefix="temp_", suffix="_phylogenetic_tree")
				st.session_state.rooted_tree.export_data(secure_temp_dir_phylogenetic_tree)
	
				
				st.info('You can visualize the phylogenetic tree at: \
					\n> * http://etetoolkit.org/treeview/ \
					\n> * https://icytree.org/\
					\n> * https://www.iroki.net/viewer\
							')
				zipfolder(secure_temp_dir_phylogenetic_tree+"/zip_phylogenetic_mafft_alignment_tree.zip", secure_temp_dir_phylogenetic_tree)

				with open(secure_temp_dir_phylogenetic_tree+"/zip_phylogenetic_mafft_alignment_tree.zip", 'rb') as f:
						
					ste.download_button(
						label="Download alignment and phylogenetic tree .zip",
						data=f,
						file_name="alignment_seqs_and_phylogenetic_tree_%s.zip" %(denoising_pipe),
						mime="application/zip")

				st.subheader('Core alpha and beta diversity metrics phylogenetic')
				st.subheader('Faith pd')
				st.table(st.session_state.core_metr_phylo.faith_pd_vector.view(pd.Series))
				
				secure_temp_dir_core_metr_phylo = tempfile.mkdtemp(prefix="temp_", suffix="_core_metr_phylo")
				secure_temp_dir_core_metr_phylo_qzv = tempfile.mkdtemp(prefix="temp_", suffix="_core_metr_phylo_qzv")

				st.session_state.core_metr_phylo.rarefied_table.export_data(secure_temp_dir_core_metr_phylo+'/rarefied_table')
				st.session_state.core_metr_phylo.faith_pd_vector.export_data(secure_temp_dir_core_metr_phylo+'/faith_pd_vector')
				st.session_state.core_metr_phylo.observed_features_vector.export_data(secure_temp_dir_core_metr_phylo+'/observed_features_vector')
				st.session_state.core_metr_phylo.shannon_vector.export_data(secure_temp_dir_core_metr_phylo+'/shannon_vector')
				st.session_state.core_metr_phylo.evenness_vector.export_data(secure_temp_dir_core_metr_phylo+'/evenness_vector')
				st.session_state.core_metr_phylo.unweighted_unifrac_distance_matrix.export_data(secure_temp_dir_core_metr_phylo+'/unweighted_unifrac_distance_matrix')
				st.session_state.core_metr_phylo.weighted_unifrac_distance_matrix.export_data(secure_temp_dir_core_metr_phylo+'/weighted_unifrac_distance_matrix')
				st.session_state.core_metr_phylo.jaccard_distance_matrix.export_data(secure_temp_dir_core_metr_phylo+'/jaccard_distance_matrix')
				st.session_state.core_metr_phylo.bray_curtis_distance_matrix.export_data(secure_temp_dir_core_metr_phylo+'/bray_curtis_distance_matrix')
				st.session_state.core_metr_phylo.unweighted_unifrac_pcoa_results.export_data(secure_temp_dir_core_metr_phylo+'/unweighted_unifrac_pcoa_results')
				st.session_state.core_metr_phylo.weighted_unifrac_pcoa_results.export_data(secure_temp_dir_core_metr_phylo+'/weighted_unifrac_pcoa_results')
				st.session_state.core_metr_phylo.jaccard_pcoa_results.export_data(secure_temp_dir_core_metr_phylo+'/jaccard_pcoa_results')
				st.session_state.core_metr_phylo.bray_curtis_pcoa_results.export_data(secure_temp_dir_core_metr_phylo+'/bray_curtis_pcoa_results')
				st.session_state.core_metr_phylo.weighted_unifrac_emperor.save(secure_temp_dir_core_metr_phylo_qzv+'/weighted_unifrac_emperor.qzv')
				st.session_state.core_metr_phylo.unweighted_unifrac_emperor.save(secure_temp_dir_core_metr_phylo_qzv+'/unweighted_unifrac_emperor.qzv')
				st.session_state.core_metr_phylo.jaccard_emperor.save(secure_temp_dir_core_metr_phylo_qzv+'/jaccard_emperor.qzv')
				st.session_state.core_metr_phylo.bray_curtis_emperor.save(secure_temp_dir_core_metr_phylo_qzv+'/bray_curtis_emperor.qzv')
				
				#zipfolder(secure_temp_dir_core_metr_phylo_qzv+"/zip_core_metr_phylo.zip", secure_temp_dir_core_metr_phylo_qzv)

				with open(secure_temp_dir_core_metr_phylo_qzv+"/bray_curtis_emperor.qzv", 'rb') as f:
						
					ste.download_button(
						label="Download visualization emperor bray curtis .qzv",
						data=f,
						file_name="bray_curtis_emperor_%s.qzv" %(denoising_pipe),
						mime="application/qzv")
				with open(secure_temp_dir_core_metr_phylo_qzv+"/jaccard_emperor.qzv", 'rb') as f:
						
					ste.download_button(
						label="Download visualization emperor jaccard .qzv",
						data=f,
						file_name="jaccard_emperor_%s.qzv" %(denoising_pipe),
						mime="application/qzv")
				with open(secure_temp_dir_core_metr_phylo_qzv+"/weighted_unifrac_emperor.qzv", 'rb') as f:
						
					ste.download_button(
						label="Download visualization emperor weightedunifrac .qzv",
						data=f,
						file_name="weighted_unifrac_emperor_%s.qzv" %(denoising_pipe),
						mime="application/qzv")
				with open(secure_temp_dir_core_metr_phylo_qzv+"/unweighted_unifrac_emperor.qzv", 'rb') as f:
						
					ste.download_button(
						label="Download visualization emperor unweighted unifrac .qzv",
						data=f,
						file_name="unweighted_unifrac_emperor_%s.qzv" %(denoising_pipe),
						mime="application/qzv")
						
				st.subheader('Taxa barplots')
				secure_temp_dir_taxa_barplots = tempfile.mkdtemp(prefix="temp_", suffix="_taxa_barplots_%s"%(denoising_pipe))

				from qiime2.plugins.taxa.visualizers import barplot
				taxa_barplot = barplot(
					table=st.session_state['rarefied_table_%s'%(denoising_pipe)],
					taxonomy=st.session_state['classification_%s'%(denoising_pipe)],
					metadata=st.session_state.metadata)
				taxa_barplot.visualization.save(secure_temp_dir_taxa_barplots+'/taxa_barplots.qzv')

				with open(secure_temp_dir_taxa_barplots+"/taxa_barplots.qzv", 'rb') as f:
						
					ste.download_button(
						label="Download taxa barplots .qzv",
						data=f,
						file_name="taxa_barplots_%s.qzv" %(denoising_pipe),
						mime="application/qzv")
					
			except Exception as e:
				
				st.exception(e)
				st.warning('Phylogenetic tree construction from representative sequences identified with %s \
					throwed an exception.' %(denoising_pipe))
				st.stop()

			
			st.balloons()

		side_plchldr7.success('***%s.*** Phylogenetic tree construction \
			\n > Tab %s. Phylogenetic tree' %(step_n, step_n))
		step_n += 1

	if not st.session_state.imported_sequences_temp_dir == sample_data_dir:
		try:
			
			shutil.rmtree(st.session_state.imported_sequences_temp_dir)
		except FileNotFoundError as e:
			# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
			# st.exception(e)
			pass
		except NameError as e:
			# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
			# st.exception(e)
			pass

	try:
		shutil.rmtree(st.session_state.secure_temp_dir_demux_summary)
	except FileNotFoundError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	except NameError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass

	try:
		shutil.rmtree(secure_temp_dir_q_filter_summary)
	except FileNotFoundError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	except NameError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass

	try:
		shutil.rmtree(secure_temp_dir_joined_summary)
	except FileNotFoundError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	except NameError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass

	try:
		shutil.rmtree(secure_temp_dir_denoising)
	except FileNotFoundError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	except NameError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	
	try:
		shutil.rmtree(secure_temp_dir_rarefaction)
	except FileNotFoundError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	except NameError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass

	try:
		shutil.rmtree(secure_temp_dir_rarefaction_beta_bc)
	except FileNotFoundError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	except NameError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass

	try:
		shutil.rmtree(secure_temp_dir_rarefaction_beta_j)
	except FileNotFoundError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	except NameError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass

	
	try:
		
		shutil.rmtree(secure_temp_dir_phylogenetic_tree)
	except FileNotFoundError as e:
		st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		st.exception(e)
		pass
	except NameError as e:
		st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		st.exception(e)
		pass


	try:
		shutil.rmtree(secure_temp_dir_core_metr_phylo)
	except FileNotFoundError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	except NameError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass

	try:
		shutil.rmtree(secure_temp_dir_core_metr_phylo_qzv)
	except FileNotFoundError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	except NameError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass

	try:
		shutil.rmtree(secure_temp_dir_taxa_barplots)
	except FileNotFoundError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass
	except NameError as e:
		# st.info('Comportamento corretto, le eccezioni sono mostrate a schermo per motivi di sviluppo dell\'applicazione.')
		# st.exception(e)
		pass

else:

	library_PE_SE_placeholder.empty()
	side_placeholder0.empty()
	step_n += 1

	side_subheader_placeholder.subheader('***Secondary analysis raw data***')
	h_plchldr.header('***Tertiary analysis of pre-processed data***')
	
side_subheader_placeholder1 = sidemenus.empty()
side_subheader_placeholder1.subheader('***Tertiary analysis pre-processed data***')

side_placeholder = sidemenus.empty()
side_placeholdera = sidemenus.empty()
side_placeholder1 = sidemenus.empty()
side_placeholder2 = sidemenus.empty()
side_placeholder3 = sidemenus.empty()
side_placeholder4 = sidemenus.empty()
side_placeholder5 = sidemenus.empty()

side_placeholder.info('***%s.*** Data upload' %(step_n))
side_placeholdera.markdown(f"<a href='#ter_linkto_{step_n}'>{step_n}. Data upload</a>", unsafe_allow_html=True)

if skip is False:
	st.markdown("""<hr style="height:8px;border:none;color:#333;background-color:#333;" /> """, unsafe_allow_html=True)
	h1_plchldr = st.empty()
	h1_plchldr = st.header('***Tertiary analysis of pre-processed data***')
	st.markdown('__ATTENTION!__ In the sidebar menu select: Skip the steps of "Secondary analysis of raw data" and go \
		directly to the "Tertiary analysis of pre-processed data" to speed-up the App.')

def delete_session_state_data_input_keys():
	try:
		del st.session_state.data_OTU_input
		del st.session_state.data_tax_input
		del st.session_state.data_meta_input
	except Exception as e:
		print(e)
		pass

sample_data_tertiary_analysis = '/app/microbiome/sample_data/tertiary_analysis_data'
st.info('Download sample data for tertiary analysis.')

with open(sample_data_tertiary_analysis+"/zip_sample_data.zip", 'rb') as f:
	
	ste.download_button(
		label="Download sample data .zip",
		data=f,
		file_name="sample_data_tertiary_analysis.zip",
		mime="application/zip")


# Form di caricamento dati
st.markdown(f"<div id='ter_linkto_{step_n}'></div>", unsafe_allow_html=True)
with st.form(key='form_data_upload'):


	#todo sistemare il sample dataframe da visualizzare nell helper del caricamento OTU file
	sample_df = pd.DataFrame(data={'Sample1': [15, 0, 0], 'Sample2 ...': [6857, 0, 0]}, index=['Otu1', 'Otu2', 'Otu3 ...'])
	sample_df.index.name = '#NAME'
	sample_df_mkdwn = sample_df.to_markdown()
	
	cols = st.columns((1,1))
	data_OTU = cols[0].file_uploader('OTU file:', 
		key='data_OTU_input', 
		help='Reads counts for each OTU/ASV. Accepted file formats: .tsv, .csv, .qza.\
			\n> Columns are samples, rows are OTUs. Header is: \
			\n> #NAME or featureID. \
			\n> __Attention!__ Upload rarefied or otherwise normalized data for correct results.' )# %(sample_df))
	#todo
	data_tax = cols[1].file_uploader('Taxonomy file:', 
		key='data_tax_input', 
		help='Taxonomic classification of OTUs/ASVs. Accepted file formats: .tsv, .csv, .qza. \
			\n> Header is: #TAXONOMY. \
			\n> Columns are taxonomic levels, rows are OTUs.')
	data_meta = cols[0].file_uploader('Metadata file (optional):', 
		key='data_meta_input', 
		help='Optional file with samples\' metadata. Accepted file formats: .tsv.\
			\n> Header is: #NAME.')

	submit_button = st.form_submit_button(
		label='Upload',
		type = 'primary'
		)
	
	delete_sessionstate_bttn_plchldr = st.empty()
	# bottone per cancellare la session state dei file taxonomy, OTU e metadata in input
	delete_sessionstate_bttn = delete_sessionstate_bttn_plchldr.form_submit_button('Modify data')
		
if delete_sessionstate_bttn:
	delete_session_state_data_input_keys()
	try:
		del st.session_state.alpha_divs
	except:
		pass
	st.experimental_rerun()
	

if submit_button or ((st.session_state.data_OTU_input is not None) and (st.session_state.data_tax_input is not None)):

	if ((st.session_state.data_OTU_input is not None) and 
		(st.session_state.data_tax_input is not None)):

		st.success('Successfully uploaded OTU and taxonomy files.')
		side_placeholder.success('***%s.*** Data upload' %(step_n))
		step_n += 1
		side_placeholder1.info('***%s.*** Set the taxonomic level' %(step_n))

	elif ((st.session_state.data_OTU_input is None) and (st.session_state.data_tax_input is not None)):

		st.warning('Select an OTU file to continue.')
		st.stop()
	elif ((st.session_state.data_OTU_input is not None) and (st.session_state.data_tax_input is None)):

		st.warning('Select a taxonomy file to continue.')
		st.stop()
	elif ((st.session_state.data_OTU_input is None) and (st.session_state.data_tax_input is None)):

		st.warning('Select both an OTU and a taxonomy file to continue.')
		st.stop()
	if (st.session_state.data_meta_input is not None):

		# Import metadata as pandas dataframe
		data_meta = pd.read_csv(data_meta, sep='\t', index_col=0)
		data_meta = data_meta[data_meta.index != '#q2:types'] # filtro la riga contenente l'intestazione del manifest file di qiime2
		data_meta = data_meta.dropna(axis=1, how='all')
		data_meta = data_meta.fillna('NA', axis=0)

		st.success('Successfully uploaded metadata file.')
	else:

		st.warning('No metadata file.')


else:

	st.stop()

# Descrizione
descr_plchldr = st.empty()

dashboard_info_dwnld_plchldr = st.empty()

dashboard_barre_dwnld_plchldr = st.empty()
dashboard_torte_dwnld_plchldr = st.empty()

descr_plchldr.markdown('Scroll the main page and navigate the tabs below to \
	visualize the results of the analyses.')


final_df, data_tax, data_OTU, taxa_levels, key_error_warning = create_final_df(st.session_state.data_OTU_input, st.session_state.data_tax_input)
if key_error_warning != '':
	st.warning(key_error_warning)
st.session_state.final_df = final_df
st.session_state.data_tax_df = data_tax
st.session_state.data_OTU_df = data_OTU
st.session_state.taxa_levels_l = taxa_levels
try:
	data_meta = data_meta.loc[st.session_state.data_OTU_df.columns,:]
	st.session_state.data_meta_df = data_meta
except Exception as e:
	# st.exception(e)
	st.warning('No metadata file.')

# Controllo corrispondenza dei campioni dei metadati coi campioni dell OTU file
#data_OTU_samples_equals_data_meta_samples = np.equal(st.session_state.data_meta_df.index, st.session_state.data_OTU_df.columns)
#data_OTU_samples_equals_data_meta_samples

# Menù laterale
# Creazione di un selettore radio del livello tassonomico preferito
tax_level = sidemenus.radio(label='Taxonomic level', options=[i for i in data_tax.columns], key='tax_level_radio',
	help='Attention! Selecting the OTU level could slow down the App.')
side_placeholder1.success('***%s.*** Select taxonomic level' %(step_n))
step_n += 1

# Creazione di un selettore radio della variabile metadati da usare per il raggruppamento dei campioni
if st.session_state.data_meta_input is not None:
	
	side_placeholder2.info('***%s.*** Select sample type' %(step_n))
	side_placeholder3.info('***%s.*** Select colour for bar charts' %(step_n+1))

	sample_grouping_opts = [i for i in data_meta.columns]
	sample_grouping_opts.append('All samples')
	sample_grouping = sidemenus.radio(label='Samples\' grouping', 
	options=sample_grouping_opts, 
	help='Metadata category to be used to group the samples. If no \
		metadata file, all samples are visualized.',
	key='sample_grouping_radio')
	st.session_state.sequenced_samples = data_meta.index.to_list()
else:

	st.session_state.sequenced_samples = data_OTU.columns.to_list()
	sample_grouping = sidemenus.radio(label='Samples\' grouping', 
	options=["All samples"], 
	help='Metadata category to be used to group the samples. If no \
		metadata file, all samples are visualized.',
	key='sample_grouping_radio')
	side_placeholder3.info('***%s.*** Select colour for bar charts' %(step_n+1))

	pass

side_placeholder2.success('***%s.*** Select sample type' %(step_n))
step_n += 1

# To be used for bar charts - 
# Creazione della lista delle Otu corrispondenti al livello tassonomico selezionato col selettore radio
# e il dizionario dei campioni corrispondenti ai gruppi nella variabile metadati selezionata col selettore radio:
# index of tax_level in list taxa_levels 
idx_tax_level = taxa_levels.index(st.session_state.tax_level_radio)
# index of sample_grouping in list samples
if st.session_state.sample_grouping_radio == 'All samples':
	
	st.session_state.idx_sample_grouping = {'All samples': st.session_state.sequenced_samples}
	
else:
	
	# dizionario {"gruppo": [ID campioni]}
	st.session_state.idx_sample_grouping = st.session_state.data_meta_df.groupby(st.session_state.sample_grouping_radio).groups

# Creazione di cinque tabs
tabs = tab_taxonomy_of_OTUs, tab_num_of_OTUs_per_taxon, tab_rel_ab, tab_alpha_div, tab_beta_div = st.tabs([
	'OTUs taxonomic classification', 
	'Number of OTUs per taxon',
	'Taxa relative abundance charts per samples\' group',
	'Alpha diversity charts',
	'Beta diversity charts'

	])



with tab_taxonomy_of_OTUs:
	
	st.markdown('***%s***' %(st.session_state.tax_level_radio))
	if st.session_state.tax_level_radio == 'OTU':
		
		st.table(data_tax.reset_index().iloc[:,:-1])
		taxa_counts_d, df1_d, tabella_df_l_l, warning_tab_num_of_OTUs_per_taxon = OTUs_annots_freqs(st.session_state.idx_sample_grouping)
		
		# stampa a schermo sull applicazione separata dalla funzione in modo tale da poter sfruttare i vantaggi del caching dei dati
		for grouping_key, grouped_samples in st.session_state.idx_sample_grouping.items():
			
			tab_num_of_OTUs_per_taxon.markdown('***%s***' %grouping_key)
			tab_num_of_OTUs_per_taxon.warning(warning_tab_num_of_OTUs_per_taxon)
		
	else:

		with st.spinner('Wait, data analysis ...'):

			taxa_counts_d, df1_d, tabella_df_l_l, warning_tab_num_of_OTUs_per_taxon = OTUs_annots_freqs(st.session_state.idx_sample_grouping)
			
			# stampa a schermo sull applicazione separata dalla funzione in modo tale da poter sfruttare i vantaggi del caching dei dati
			for i, (grouping_key, grouped_samples) in enumerate(st.session_state.idx_sample_grouping.items()):

				st.markdown('***%s***' %grouping_key)
				for number, taxon_df in enumerate(tabella_df_l_l[i]):
					
					taxon = taxon_df.columns.name
					expander_string = ((' '.join([str(number+1), '. ', taxon])))
					with st.expander(expander_string):
						
						st.bar_chart(taxon_df)

						st.dataframe(taxon_df)
						
						csv = convert_df(taxon_df)
						ste.download_button(
							label="Download CSV table",
							data=csv,
							file_name=('OTUs_per_%s_%s.csv' %(st.session_state.tax_level_radio, taxon)),
							mime='text/csv')

		st.balloons()
		
with tab_num_of_OTUs_per_taxon:

	if tax_level != 'OTU':

		with st.spinner('Wait, data analysis ...'):
			
			st.markdown('***%s***' %(st.session_state.tax_level_radio))
			for grouping_key, grouped_samples in st.session_state.idx_sample_grouping.items():

				st.markdown('***%s***' %(grouping_key))
				#st.dataframe(pd.DataFrame(taxa_counts_d[grouping_key]).T)
				df = pd.DataFrame(np.diag(pd.DataFrame(taxa_counts_d[grouping_key]).T), index = (pd.DataFrame(taxa_counts_d[grouping_key]).T.index)).reset_index(drop=False)
				df = df.rename({'index': ('Taxonomic classification level %s' %(st.session_state.tax_level_radio)), 0:'number of OTUs'}, axis=1)
				
				c = alt.Chart(df).mark_circle().encode(
					x=('Taxonomic classification level %s' %(st.session_state.tax_level_radio)), y='number of OTUs', 
					size='number of OTUs', color='number of OTUs', tooltip=[('Taxonomic classification level %s' %(st.session_state.tax_level_radio)), 'number of OTUs'])

				st.altair_chart(c, use_container_width=True)
				
				# Riformattazione della tabella
				
				tabella_df = df.sort_values(by='number of OTUs', ascending=False).reset_index(drop=True)
				st.table(tabella_df.style.format('{:,.0f}', subset= ['number of OTUs']))
				csv = convert_df(tabella_df)
				ste.download_button(
					label="Download CSV table",
					data=csv,
					file_name=('number_of_OTUs_per_%s.csv' %(st.session_state.tax_level_radio)),
					mime='text/csv')

			st.balloons()



with tab_rel_ab:

	st.header('Bar charts')
	with st.spinner('Wait, data analysis ...'):

		# Creazione dei grafici a barre impilate delle conte delle reads delle Otu nei campioni per ogni 
		# annotazione tassonomica del livello selezionato e con indicazione del Phylum corrispondente
		st.markdown('***%s*** - ***%s***' %(
			st.session_state.tax_level_radio, 
			st.session_state.sample_grouping_radio))

		
		if st.session_state.tax_level_radio == 'Phylum':
			
			
			stacked_grouped_bars_radio_plchldr = sidemenus.empty()
			side_form_color_plchldr = sidemenus.empty()
			color_picker_plchldr = st.empty()
			with side_form_color_plchldr.form(key='color_picker_phylum_form'):

				for grouping_key, grouped_samples in st.session_state.idx_sample_grouping.items():
					
					grouped_samples = list(grouped_samples)
								
					if ((st.session_state.data_meta_input is not None) and (grouping_key != 'All samples')):

						color_picker_col = st.radio('Bars colour for group %s' %(grouping_key), options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
						key='color_picker_col_%s_%s_radio'%(tax_level, grouping_key))#+[data_meta.index.name])
						submit_button_label = ''
					elif (st.session_state.data_meta_input is None):

						color_picker_col = None
						st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)] = None
						submit_button_label = 'No metadata file. One default different bar colour \
							for each sample.'
						
						palette_seq_or_qual = st.radio('Bar charts colour', 
						options=['Sequential', 'Gradual'],
						help='Sequential: one different colour for each bar representing a sample (max. 25 colours repeated cyclically). \
							Gradual: rainbow colours shades gradually applied for all the bars representing the samples.',
						key='palette_%s_%s_radio'%(tax_level, grouping_key))
						color_picker_col = palette_seq_or_qual
						st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)] = st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)]
						
						side_placeholder3.success('***%s.*** Select colour for barcharts' %(step_n))

						#palette = cycle(px.colors.qualitative.Alphabet)
					elif ((st.session_state.data_meta_input is not None) and (grouping_key == 'All samples')):
						
						stacked_grouped_bars_button_plchldr = st.columns((1,1))
						stacked_grouped_bars_radio = stacked_grouped_bars_radio_plchldr.radio('Two \
							different types of visualizations and colouring of barcharts available:', options= ['Stacked barcharts', 'Grouped barcharts'], help='Select the desired visualization: \
								\n> Stacked barcharts: to select the sequential or gradual colouring of stacked barcharts.\
								\n> Grouped barcharts: to select the colouting of grouped barcharts based on a metadata category.',
								key='stacked_grouped_bars_radio',
								index=0)
						
						color_radio_plchldr = st.empty()
						# color_picker_col = color_radio_plchldr.radio('Colore delle barre dei grafici per il gruppo %s' %(grouping_key), 
						# 	options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
						# 	key='color_picker_col_%s_%s_radio' %(tax_level, grouping_key))#+[data_meta.index.name])
						
						submit_button_label = 'NOTE: No samples\' grouping.'

						if st.session_state.stacked_grouped_bars_radio == 'Stacked barcharts':
							
							color_picker_col = None
							try:
								del st.session_state['color_picker_col_%s_%s_radio' %(tax_level, grouping_key)]
							except Exception as e:
								pass
							
							palette_seq_or_qual = color_radio_plchldr.radio('Bar charts colour', 
							options=['Sequential', 'Gradual'],
							help='Sequential: one different colour for each bar representing a sample (max. 25 colours repeated cyclically). \
								Gradual: rainbow colours shades gradually applied for all the bars representing the samples.',
							key='palette_%s_%s_radio'%(tax_level, grouping_key))
							
							submit_button_label = 'NOTE: No samples\' grouping. One default different bar colour \
							for each sample.'
						elif st.session_state.stacked_grouped_bars_radio == 'Grouped barcharts':

							color_picker_col = None

							color_picker_col = color_radio_plchldr.radio('Bars colour for group %s' %(grouping_key), 
								options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
								key='color_picker_col_%s_%s_radio' %(tax_level, grouping_key))#+[data_meta.index.name])
							
							try:
								del st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)]			
							except Exception as e:
								pass
							submit_button_label = 'NOTE: No samples\' grouping.'

				
				try:
					st.markdown(submit_button_label)
				except:
					pass
				submit_button = st.form_submit_button('OK, colour', type='primary')
			if ((color_picker_col is not None) and (st.session_state['color_picker_col_%s_%s_radio' %(tax_level, grouping_key)] == st.session_state.sample_grouping_radio)):

				color_picked = color_picker_plchldr.color_picker('Bar charts colour:', value = '#00FFAA', key = 'color_picked')
			if (submit_button or (color_picker_col is not None)) :
				
				side_placeholder3.success('***%s.*** Select colour for barcharts' %(step_n))
				for grouping_key, grouped_samples in st.session_state.idx_sample_grouping.items():
					grouped_samples = list(grouped_samples)
					

					st.markdown('***%s***' %grouping_key)
					
					df = df1_d[grouping_key]
					
					if ((color_picker_col is not None) and (data_meta is not None)):
						
						barmode = 'group'
						palette = cycle(px.colors.qualitative.Alphabet)
						# palette = cycle(px.colors.sequential.PuBu)
						colors_d = data_meta.loc[grouped_samples, :].groupby(st.session_state['color_picker_col_%s_%s_radio' %(tax_level, grouping_key)]).groups
						colors_d_inv = {i:v for v,c in colors_d.items() for i in c}
						
						colors0 = {}
						vs = []
						for ind, (c, v) in enumerate(colors_d_inv.items()):
							
							vs.append(v)
							if ind > 0:
								if v != vs[ind-1]:
									colors0[c] = (next(palette), v)
								else:
									colors0[c] = (None, v)
							else:
								colors0[c] = (next(palette), v)

						if None in [i[0] for i in colors0.values()]:
							colors1 = pd.DataFrame(pd.DataFrame(colors0, index=[0,1], columns=colors0.keys()).T)
							
							colors = colors1.ffill()
							
							#colors.reset_index(drop=True, inplace=True)
							colors = {v[1]:v[0] for c,v in colors.iterrows()}
							
						else:
							colors = {v[1]:v[0] for c,v in colors0.items()}
						
					else:
						
						barmode = 'stack'
						if st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)] == 'Sequential':
							palette = cycle(px.colors.qualitative.Alphabet)
						
						elif st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)] == 'Gradual':
							# sequenza di colori personalizzata fra blu e rosso basata sulla quantita di campioni
							n_colors = len(data_OTU.columns)+1
							custom_colors_grad = px.colors.sample_colorscale("turbo", [n/(n_colors -1) for n in range(n_colors)])
							palette = cycle(custom_colors_grad)
						side_placeholder3.success('***%s.*** Select colour for barcharts' %(step_n))

						colors = {c: next(palette) for c in data_OTU.columns}

	
					# subplot setup
					fig = make_subplots(rows=1, cols=1)
					
					
					# add bars, stesso colore per tutte le barre dei campioni di ciascun gruppo in grouping_key
					counter = []
					
					for cols in df.loc[:, df.columns.isin(grouped_samples)]:
						
						if st.session_state.data_meta_input is None:
							leg_gr_col = cols
						else:
							try:
								leg_gr_col = data_meta.loc[cols, st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)]]
							except: # Se Raggruppameto dei campioni e' Tutti i campioni
								leg_gr_col = cols
										
							#leg_gr_col = colors_d_inv[cols]
						
						if colors[leg_gr_col] not in counter:
							showlegend = True
						else:
							showlegend = False
						
						if 'color_picked' in st.session_state.keys():
							marker_color = st.session_state.color_picked
						else:
							marker_color = colors[leg_gr_col]

						try:
							
							fig.add_trace(go.Bar(
								y = [df[taxa_levels[idx_tax_level-1]],
								df[st.session_state.tax_level_radio]],
								x = df[cols],
								name = leg_gr_col,
								showlegend = showlegend,
								legendgroup = cols,
								orientation = 'h',
								marker_color = marker_color,
								offsetgroup = leg_gr_col), 
							row = 1, col = 1)
							counter.append(colors[leg_gr_col])
						except ValueError:
							st.warning('Choose a different metadata category for barcharts colour')
							st.stop()
					title_plot = len(grouped_samples)
					fig.update_layout(title = 'N = %s'%(title_plot), barmode=barmode, plot_bgcolor = None, bargroupgap = 0.1)

					st.plotly_chart(figure_or_data=fig, height=900, use_container_width=True, config=config)
				step_n += 1
			else:

				st.warning('The page is awaiting for you to pick a colour for barcharts for each \
					samples\' grouping from the sidebar menu.')
			
				step_n += 1

		elif st.session_state.tax_level_radio == 'Kingdom':

			
			stacked_grouped_bars_radio_plchldr = sidemenus.empty()
			color_picker_plchldr = st.empty()
			with sidemenus.form(key='color_picker_form_%s'%(tax_level)):

				for grouping_key, grouped_samples in st.session_state.idx_sample_grouping.items():
					
					if ((st.session_state.data_meta_input is not None) and (grouping_key != 'All samples')):

						color_picker_col = st.radio('Bars colour for group %s' %(grouping_key), options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
						key='color_picker_col_%s_%s_radio' %(tax_level, grouping_key))#+[data_meta.index.name])
						submit_button_label = ''
					elif (st.session_state.data_meta_input is None):

						color_picker_col = None
						st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)] = None
						submit_button_label = 'No metadata file. One default different bar colour \
							for each sample.'
				
						palette_seq_or_qual = st.radio('Bar charts colour', 
						options=['Sequential', 'Gradual'],
						help='Sequential: one different colour for each bar representing a sample (max. 25 colours repeated cyclically). \
								Gradual: rainbow colours shades gradually applied for all the bars representing the samples.',
						key='palette_%s_%s_radio'%(tax_level, grouping_key))
						color_picker_col = palette_seq_or_qual
						st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)] = st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)]
					elif ((st.session_state.data_meta_input is not None) and (grouping_key == 'All samples')):
						
						stacked_grouped_bars_button_plchldr = st.columns((1,1))
						stacked_grouped_bars_radio = stacked_grouped_bars_radio_plchldr.radio('Two \
							different types of visualizations and colouring of barcharts available:', options= ['Stacked barcharts', 'Grouped barcharts'], help='Select the desired visualization: \
								\n> Stacked barcharts: to select the sequential or gradual colouring of stacked barcharts.\
								\n> Grouped barcharts: to select the colouting of grouped barcharts based on a metadata category.',
								key='stacked_grouped_bars_radio',
								index=0)
						
						color_radio_plchldr = st.empty()
						# color_picker_col = color_radio_plchldr.radio('Colore delle barre dei grafici per il gruppo %s' %(grouping_key), 
						# 	options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
						# 	key='color_picker_col_%s_%s_radio' %(tax_level, grouping_key))#+[data_meta.index.name])
						
						submit_button_label = 'NOTE: No samples\' grouping.'
						if stacked_grouped_bars_radio == 'Stacked barcharts':
							
							color_picker_col = None
							try:
								del st.session_state['color_picker_col_%s_%s_radio' %(tax_level, grouping_key)]
							except Exception as e:
								pass
							
							palette_seq_or_qual = color_radio_plchldr.radio('Bar charts colour', 
							options=['Sequential', 'Gradual'],
							help='Sequential: one different colour for each bar representing a sample (max. 25 colours repeated cyclically). \
								Gradual: rainbow colours shades gradually applied for all the bars representing the samples.',
							key='palette_%s_%s_radio'%(tax_level, grouping_key))
							
							submit_button_label = 'NOTE: No samples\' grouping. One default different bar colour \
							for each sample.'
						elif stacked_grouped_bars_radio == 'Grouped barcharts':

							color_picker_col = None

							color_picker_col = color_radio_plchldr.radio('Bars colour for group %s' %(grouping_key), 
								options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
								key='color_picker_col_%s_%s_radio' %(tax_level, grouping_key))#+[data_meta.index.name])
							
							try:
								del st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)]			
							except Exception as e:
								pass
							submit_button_label = 'NOTE: No samples\' grouping.'

				try:
					st.markdown(submit_button_label)
				except:
					pass
				submit_button = st.form_submit_button('OK, colour', type = 'primary')
			if ((color_picker_col is not None) and (color_picker_col == sample_grouping)):

				color_picked = color_picker_plchldr.color_picker('Bar charts colour:', value = '#00FFAA', key = 'color_picked')
			if (submit_button or (color_picker_col is not None)):
			
				side_placeholder3.success('***%s.*** Select colour for barcharts' %(step_n))
				for grouping_key, grouped_samples in st.session_state.idx_sample_grouping.items():				
					
					if ((color_picker_col is not None) and (data_meta is not None)):
					
						barmode = 'group'
						palette = cycle(px.colors.qualitative.Alphabet)
						# palette = cycle(px.colors.sequential.PuBu)
						colors_d = data_meta.loc[grouped_samples, :].groupby(st.session_state['color_picker_col_%s_%s_radio' %(tax_level, grouping_key)]).groups
						colors_d_inv = {i:v for v,c in colors_d.items() for i in c}
						
						
						colors0 = {}
						vs = []
						for ind, (c, v) in enumerate(colors_d_inv.items()):
							
							vs.append(v)
							if ind > 0:
								if v != vs[ind-1]:
									colors0[c] = (next(palette), v)
								else:
									colors0[c] = (None, v)
							else:
								colors0[c] = (next(palette), v)

						if None in [i[0] for i in colors0.values()]:
							colors1 = pd.DataFrame(pd.DataFrame(colors0, index=[0,1], columns=colors0.keys()).T)
							
							colors = colors1.ffill()
							
							#colors.reset_index(drop=True, inplace=True)
							colors = {v[1]:v[0] for c,v in colors.iterrows()}
							
						else:
							colors = {v[1]:v[0] for c,v in colors0.items()}
					else:
						
						barmode = 'stack'
						if st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)] == 'Sequential':
							palette = cycle(px.colors.qualitative.Alphabet)
						elif st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)] == 'Gradual':
							# sequenza di colori personalizzata fra blu e rosso basata sulla quantita di campioni
							n_colors = len(data_OTU.columns)+1
							custom_colors_grad = px.colors.sample_colorscale("turbo", [n/(n_colors -1) for n in range(n_colors)])
							palette = cycle(custom_colors_grad)
						side_placeholder3.success('***%s.*** Select colour for barcharts' %(step_n))

						colors = {c: next(palette) for c in data_OTU.columns}

					# subplot setup
					fig = make_subplots(rows=1, cols=1)

					st.markdown('***%s***' %grouping_key)

					df = df1_d[grouping_key]
					
					grouped_samples = list(grouped_samples)
					
					# add bars, stesso colore per tutte le barre dei campioni di ciascun gruppo in grouping_key
					counter = []
					
					for cols in df[grouped_samples]:
						if st.session_state.data_meta_input is None:
							leg_gr_col = cols
						else:
							try:
								leg_gr_col = data_meta.loc[cols, st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)]]
							except: # Se Raggruppameto dei campioni e' Tutti i campioni
								leg_gr_col = cols
								
							#leg_gr_col = colors_d_inv[cols]
						
						if colors[leg_gr_col] not in counter:
							showlegend = True
						else:
							showlegend = False
							
						if 'color_picked' in st.session_state.keys():
							marker_color = st.session_state.color_picked
						else:
							marker_color = colors[leg_gr_col]

						try:
							
							fig.add_trace(go.Bar(
								y = df[st.session_state.tax_level_radio],
								x = df[cols],
								name = leg_gr_col,
								showlegend = showlegend,
								legendgroup = cols,
								orientation = 'h',
								marker_color = marker_color,
								offsetgroup = leg_gr_col), 
							row = 1, col = 1)
							counter.append(colors[leg_gr_col])
						except ValueError:
							st.warning('Choose a different metadata category for barcharts colour')
							st.stop()
					title_plot = len(grouped_samples)
					fig.update_layout(title = 'N = %s'%(title_plot), barmode=barmode, plot_bgcolor = None, bargroupgap = 0.1)

					st.plotly_chart(figure_or_data=fig, height=900, use_container_width=True, config=config)
				step_n += 1
			else:
				
				st.warning('The page is awaiting for you to pick a colour for barcharts for each \
					samples\' grouping from the sidebar menu.')
			
				step_n += 1

		elif st.session_state.tax_level_radio == 'OTU':
			
			
			stacked_grouped_bars_radio_plchldr = sidemenus.empty()
			color_picker_plchldr = st.empty()
			with sidemenus.form(key='color_picker_form_%s'%(tax_level)):

				for grouping_key, grouped_samples in st.session_state.idx_sample_grouping.items():
					
					grouped_samples = list(grouped_samples)
					
					if ((st.session_state.data_meta_input is not None) and (grouping_key != 'All samples')):

						color_picker_col = st.radio('Bars colour for group %s' %(grouping_key), options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
						key='color_picker_col_%s_%s_radio'%(tax_level, grouping_key))#+[data_meta.index.name])
						submit_button_label = ''
					elif (st.session_state.data_meta_input is None):

						color_picker_col = None
						st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)] = None
						submit_button_label = 'No metadata file. One default different bar colour \
							for each sample.'
				
						palette_seq_or_qual = st.radio('Bar charts colour', 
						options=['Sequential', 'Gradual'],
						help='Sequential: one different colour for each bar representing a sample (max. 25 colours repeated cyclically). \
								Gradual: rainbow colours shades gradually applied for all the bars representing the samples.',
						key='palette_%s_%s_radio'%(tax_level, grouping_key))
						color_picker_col = palette_seq_or_qual
						st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)] = st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)]
					elif ((st.session_state.data_meta_input is not None) and (grouping_key == 'All samples')):
						
						stacked_grouped_bars_button_plchldr = st.columns((1,1))
						stacked_grouped_bars_radio = stacked_grouped_bars_radio_plchldr.radio('Two \
							different types of visualizations and colouring of barcharts available:', options= ['Stacked barcharts', 'Grouped barcharts'], help='Select the desired visualization: \
								\n> Stacked barcharts: to select the sequential or gradual colouring of stacked barcharts.\
								\n> Grouped barcharts: to select the colouting of grouped barcharts based on a metadata category.',
								key='stacked_grouped_bars_radio',
								index=0)
						
						color_radio_plchldr = st.empty()
						# color_picker_col = color_radio_plchldr.radio('Colore delle barre dei grafici per il gruppo %s' %(grouping_key), 
						# 	options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
						# 	key='color_picker_col_%s_%s_radio' %(tax_level, grouping_key))#+[data_meta.index.name])
						
						submit_button_label = 'NOTE: No samples\' grouping.'
						if stacked_grouped_bars_radio == 'Stacked barcharts':
							
							color_picker_col = None
							try:
								del st.session_state['color_picker_col_%s_%s_radio' %(tax_level, grouping_key)]
							except Exception as e:
								pass
							
							palette_seq_or_qual = color_radio_plchldr.radio('Bar charts colour', 
							options=['Sequential', 'Gradual'],
							help='Sequential: one different colour for each bar representing a sample (max. 25 colours repeated cyclically). \
								Gradual: rainbow colours shades gradually applied for all the bars representing the samples.',
							key='palette_%s_%s_radio'%(tax_level, grouping_key))
							
							submit_button_label = 'NOTE: No samples\' grouping. One default different bar colour \
							for each sample.'
						elif stacked_grouped_bars_radio == 'Grouped barcharts':

							color_picker_col = None

							color_picker_col = color_radio_plchldr.radio('Bars colour for group %s' %(grouping_key), 
								options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
								key='color_picker_col_%s_%s_radio' %(tax_level, grouping_key))#+[data_meta.index.name])
							
							try:
								del st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)]			
							except Exception as e:
								pass
							submit_button_label = 'NOTE: No samples\' grouping.'

				try:
					st.markdown(submit_button_label)
				except:
					pass
				submit_button = st.form_submit_button('OK, colour', type = 'primary')
			if ((color_picker_col is not None) and (color_picker_col == sample_grouping)):

				color_picked = color_picker_plchldr.color_picker('Bar charts colour:', value = '#00FFAA', key = 'color_picked')
			if (submit_button or (color_picker_col is not None)):
				
				side_placeholder3.success('***%s.*** Select colour for barcharts' %(step_n))
				for grouping_key, grouped_samples in st.session_state.idx_sample_grouping.items():
					
					grouped_samples = list(grouped_samples)
					st.markdown('***%s***' %grouping_key)
						
					df = df1_d[grouping_key]
					
					if ((color_picker_col is not None) and (data_meta is not None)):
					
						barmode = 'group'
						palette = cycle(px.colors.qualitative.Alphabet)
						# palette = cycle(px.colors.sequential.PuBu)
						colors_d = data_meta.loc[grouped_samples, :].groupby(st.session_state['color_picker_col_%s_%s_radio' %(tax_level, grouping_key)]).groups
						colors_d_inv = {i:v for v,c in colors_d.items() for i in c}
						
						
						colors0 = {}
						vs = []
						for ind, (c, v) in enumerate(colors_d_inv.items()):
							
							vs.append(v)
							if ind > 0:
								if v != vs[ind-1]:
									colors0[c] = (next(palette), v)
								else:
									colors0[c] = (None, v)
							else:
								colors0[c] = (next(palette), v)

						if None in [i[0] for i in colors0.values()]:
							colors1 = pd.DataFrame(pd.DataFrame(colors0, index=[0,1], columns=colors0.keys()).T)
							
							colors = colors1.ffill()
							
							#colors.reset_index(drop=True, inplace=True)
							colors = {v[1]:v[0] for c,v in colors.iterrows()}
							
						else:
							colors = {v[1]:v[0] for c,v in colors0.items()}
					else:
						
						barmode = 'stack'
						if st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)] == 'Sequential':
							palette = cycle(px.colors.qualitative.Alphabet)
						elif st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)] == 'Gradual':
							# sequenza di colori personalizzata fra blu e rosso basata sulla quantita di campioni
							n_colors = len(data_OTU.columns)+1
							custom_colors_grad = px.colors.sample_colorscale("turbo", [n/(n_colors -1) for n in range(n_colors)])
							palette = cycle(custom_colors_grad)
						side_placeholder3.success('***%s.*** Select colour for barcharts' %(step_n))

						colors = {c: next(palette) for c in data_OTU.columns}

					
					# subplot setup
					fig = make_subplots(rows=1, cols=1)
					
					
					# add bars, stesso colore per tutte le barre dei campioni di ciascun gruppo in grouping_key
					counter = []
					
					# colors_d_inv

					# df
					# idx_tax_level-1
					for cols in df.loc[:, df.columns.isin(grouped_samples)]:
						
						if st.session_state.data_meta_input is None:
							leg_gr_col = cols
						else:
							try:
								leg_gr_col = data_meta.loc[cols, st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)]]
							except: # Se Raggruppameto dei campioni e' Tutti i campioni
								leg_gr_col = cols
								
							#leg_gr_col = colors_d_inv[cols]
						
						if colors[leg_gr_col] not in counter:
							showlegend = True
						else:
							showlegend = False

						if 'color_picked' in st.session_state.keys():
							marker_color = st.session_state.color_picked
						else:
							marker_color = colors[leg_gr_col]

						try:
							
							fig.add_trace(go.Bar(
								y = [df[taxa_levels[idx_tax_level-2]],
								df[st.session_state.tax_level_radio]],
								x = df[cols],
								name = leg_gr_col,
								showlegend = showlegend,
								legendgroup = cols,
								orientation = 'h',
								marker_color = marker_color,
								offsetgroup = leg_gr_col), 
							row = 1, col = 1)
							counter.append(colors[leg_gr_col])
						except ValueError:
							st.warning('Choose a different metadata category for barcharts colour')
							st.stop()
					title_plot = len(grouped_samples)
					fig.update_layout(title = 'N = %s'%(title_plot), barmode=barmode, plot_bgcolor = None, bargroupgap = 0.1)

					st.plotly_chart(figure_or_data=fig, height=900, use_container_width=True, config=config)
				step_n += 1
			
			else:
				
				st.warning('The page is awaiting for you to pick a colour for barcharts for each \
					samples\' grouping from the sidebar menu.')
			
				step_n += 1
		else:
			
			stacked_grouped_bars_radio_plchldr = sidemenus.empty()
			color_picker_plchldr = st.empty()
			with sidemenus.form(key='color_picker_form_anytax'):

				for grouping_key, grouped_samples in st.session_state.idx_sample_grouping.items():
					
					if ((st.session_state.data_meta_input is not None) and (grouping_key != 'All samples')):

						color_picker_col = st.radio('Bars colour for group %s' %(grouping_key), options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
						key='color_picker_col_%s_%s_radio' %(tax_level, grouping_key))#+[data_meta.index.name])
						submit_button_label = ''
					elif (st.session_state.data_meta_input is None):
						
						color_picker_col = None
						st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)] = None
						submit_button_label = 'No metadata file. One default different bar colour \
							for each sample.'
						
						palette_seq_or_qual = st.radio('Bar charts colour', 
						options=['Sequential', 'Gradual'],
						help='Sequential: one different colour for each bar representing a sample (max. 25 colours repeated cyclically). \
								Gradual: rainbow colours shades gradually applied for all the bars representing the samples.',
						key='palette_%s_%s_radio'%(tax_level, grouping_key))
						color_picker_col = palette_seq_or_qual
						st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)] = st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)]
					elif ((st.session_state.data_meta_input is not None) and (grouping_key == 'All samples')):
						
						stacked_grouped_bars_button_plchldr = st.columns((1,1))
						stacked_grouped_bars_radio = stacked_grouped_bars_radio_plchldr.radio('Two \
							different types of visualizations and colouring of barcharts available:', options= ['Stacked barcharts', 'Grouped barcharts'], help='Select the desired visualization: \
								\n> Stacked barcharts: to select the sequential or gradual colouring of stacked barcharts.\
								\n> Grouped barcharts: to select the colouting of grouped barcharts based on a metadata category.',
								key='stacked_grouped_bars_radio',
								index=0)
						
						color_radio_plchldr = st.empty()
						# color_picker_col = color_radio_plchldr.radio('Colore delle barre dei grafici per il gruppo %s' %(grouping_key), 
						# 	options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
						# 	key='color_picker_col_%s_%s_radio' %(tax_level, grouping_key))#+[data_meta.index.name])
						
						submit_button_label = 'NOTE: No samples\' grouping.'
						if stacked_grouped_bars_radio == 'Stacked barcharts':
							
							color_picker_col = None
							try:
								del st.session_state['color_picker_col_%s_%s_radio' %(tax_level, grouping_key)]
							except Exception as e:
								pass
							
							palette_seq_or_qual = color_radio_plchldr.radio('Bar charts colour', 
							options=['Sequential', 'Gradual'],
							help='Sequential: one different colour for each bar representing a sample (max. 25 colours repeated cyclically). \
								Gradual: rainbow colours shades gradually applied for all the bars representing the samples.',
							key='palette_%s_%s_radio'%(tax_level, grouping_key))
							
							submit_button_label = 'NOTE: No samples\' grouping. One default different bar colour \
							for each sample.'
						elif stacked_grouped_bars_radio == 'Grouped barcharts':

							color_picker_col = None
							
							color_picker_col = color_radio_plchldr.radio('Bars colour for group %s' %(grouping_key), 
								options=[i for i in data_meta.loc[grouped_samples,:].columns.values], 
								key='color_picker_col_%s_%s_radio' %(tax_level, grouping_key))#+[data_meta.index.name])
							
							try:
								del st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)]			
							except Exception as e:
								pass
							submit_button_label = 'NOTE: No samples\' grouping.'
	
				try:
					st.markdown(submit_button_label)
				except:
					pass
				submit_button = st.form_submit_button('OK, colour', type = 'primary')
			if ((color_picker_col is not None) and (color_picker_col == sample_grouping)):

				color_picked = color_picker_plchldr.color_picker('Bar charts colour:', value = '#00FFAA', key = 'color_picked')
			if (submit_button or (color_picker_col is not None)):
				
				side_placeholder3.success('***%s.*** Select colour for barcharts' %(step_n))
				for grouping_key, grouped_samples in st.session_state.idx_sample_grouping.items():
					
					df = df1_d[grouping_key]
					
					if ((color_picker_col is not None) and (data_meta is not None)):
					
						barmode = 'group'
						palette = cycle(px.colors.qualitative.Alphabet)
						# palette = cycle(px.colors.sequential.PuBu)
						colors_d = data_meta.loc[grouped_samples, :].groupby(st.session_state['color_picker_col_%s_%s_radio' %(tax_level, grouping_key)]).groups
						colors_d_inv = {i:v for v,c in colors_d.items() for i in c}
						
						
						colors0 = {}
						vs = []
						for ind, (c, v) in enumerate(colors_d_inv.items()):
							
							vs.append(v)
							if ind > 0:
								if v != vs[ind-1]:
									colors0[c] = (next(palette), v)
								else:
									colors0[c] = (None, v)
							else:
								colors0[c] = (next(palette), v)

						if None in [i[0] for i in colors0.values()]:
							colors1 = pd.DataFrame(pd.DataFrame(colors0, index=[0,1], columns=colors0.keys()).T)
							
							colors = colors1.ffill()
							
							#colors.reset_index(drop=True, inplace=True)
							colors = {v[1]:v[0] for c,v in colors.iterrows()}
							
						else:
							colors = {v[1]:v[0] for c,v in colors0.items()}
					else:
						
						barmode = 'stack'
						if st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)] == 'Sequential':
							palette = cycle(px.colors.qualitative.Alphabet)
						elif st.session_state['palette_%s_%s_radio'%(tax_level, grouping_key)] == 'Gradual':
							# sequenza di colori personalizzata fra blu e rosso basata sulla quantita di campioni
							n_colors = len(data_OTU.columns)+1
							custom_colors_grad = px.colors.sample_colorscale("turbo", [n/(n_colors -1) for n in range(n_colors)])
							palette = cycle(custom_colors_grad)
						side_placeholder3.success('***%s.*** Select colour for barcharts' %(step_n))

						colors = {c: next(palette) for c in data_OTU.columns}

					
					# subplot setup
					fig = make_subplots(rows=1, cols=1)

					st.markdown('***%s***' %grouping_key)

					grouped_samples = list(grouped_samples)
					# add bars
					counter = []
					
					for cols in df.loc[:, grouped_samples]:
						
						if st.session_state.data_meta_input is None:
							leg_gr_col = cols
						else:
							try:
								leg_gr_col = data_meta.loc[cols, st.session_state['color_picker_col_%s_%s_radio'%(tax_level, grouping_key)]]
							except: # Se Raggruppameto dei campioni e' Tutti i campioni
								leg_gr_col = cols

							#leg_gr_col = colors_d_inv[cols]
						
						if colors[leg_gr_col] not in counter:
							showlegend = True
						else:
							showlegend = False
						
						if 'color_picked' in st.session_state.keys():
							marker_color = st.session_state.color_picked
						else:
							marker_color = colors[leg_gr_col]
# todo treemap livello anytax
						try:
							
							fig.add_trace(go.Bar(
								y = [df['Phylum'],
								df[st.session_state.tax_level_radio]],
								x = df[cols],
								name = leg_gr_col,
								showlegend = showlegend, 
								legendgroup = cols,
								orientation = 'h',
								marker_color = marker_color,
								offsetgroup = leg_gr_col
								), 
							row = 1, col = 1)
							counter.append(colors[leg_gr_col])
						except ValueError:
							st.warning('Choose a different metadata category for barcharts colour')
							st.stop()
					title_plot = len(grouped_samples)
					fig.update_layout(title = 'N = %s'%(title_plot), barmode=barmode, plot_bgcolor = None, bargroupgap = 0.1)
					
					st.plotly_chart(figure_or_data=fig, height=900, use_container_width=True, config=config)
			
							
					# from IPython.core import display
					# r_fig = roughviz.roughviz.stackedbar(df[st.session_state.tax_level_radio], df[grouped_samples],plot_svg=True)

					# html_obj = display.HTML(r_fig)
					# raw_html = html_obj._repr_html_()
					# components.html(raw_html)

				step_n += 1	
			else:
				
				st.warning('The page is awaiting for you to pick a colour for barcharts for each \
					samples\' grouping from the sidebar menu.')
			
				step_n += 1		
		

		###################################### venn diagrams #########################
		# definisco i set di otu per ogni raggruppamento di campioni basato sui metadati per fare i diagrammi
		# di venn online sul sito https://bioinformatics.psb.ugent.be/webtools/Venn/:
		if 'sample_grouping_radio' in st.session_state.keys():
			
			# deifinisco il df dei campioni raggruppati per la variabile dei metadati scelta con le annotazioni delle OTU tassonomiche
			try: # if metadata are present
				df_sunburst = st.session_state.final_df.iloc[:,8:].T.merge(st.session_state.data_meta_df, left_index=True, right_index=True)
				df_sunburst.columns = [i for i in list(st.session_state.final_df.OTU.to_list() + st.session_state.data_meta_df.columns.to_list())]
			except:
				st.warning('No metadata file')
				df_sunburst = st.session_state.final_df.iloc[:,8:].T
				df_sunburst.columns = [i for i in list(st.session_state.final_df.OTU.to_list())]
				
			try:
				df_sunburst = df_sunburst.groupby(st.session_state.sample_grouping_radio).sum()
			except:
				df_sunburst = df_sunburst

			df_sunburst = df_sunburst.T
			# deifinisco il df dei taxa raggruppati per il livello tassonomico scelto per ogni raggruppamento dei campioni
			df_sunburst = df_sunburst.merge(st.session_state.final_df.loc[:, taxa_levels], left_index=True, right_on='OTU')
			df_venn = df_sunburst.iloc[:, :-8]
			df_venn[st.session_state.tax_level_radio] = df_sunburst[st.session_state.tax_level_radio]
			df_venn = df_venn.groupby(st.session_state.tax_level_radio, axis=0).sum()
			st.info('Downloadable __final table__ of samples grouped per the selected _metadata category_ \
				of taxa grouped per the selected _taxonomic level_. Format: .csv.')
			tabella_df = df_venn
			
			csv = convert_df(tabella_df)
			ste.download_button(
				label="Download table of samples grouped per %s of taxonomic level %s" %(st.session_state.sample_grouping_radio, st.session_state.tax_level_radio),
				data=csv,
				file_name= 'final_table_%s_%s_%s.csv' %(st.session_state.dashboard_name, st.session_state.sample_grouping_radio, st.session_state.tax_level_radio),
				mime="text/csv")
			# definisco le conte percentuali
			
			
			st.session_state.tabella_df_filtered = filter_dataframe(tabella_df.reset_index(level=0))
			st.write(st.session_state.tabella_df_filtered)
			csv = convert_df(st.session_state.tabella_df_filtered)
			ste.download_button(
				label="Download final table filtered of samples grouped per %s of taxonomic level %s" %(st.session_state.sample_grouping_radio, st.session_state.tax_level_radio),
				data=csv,
				file_name= 'final_table_filtered_%s_%s_%s.csv' %(st.session_state.dashboard_name, st.session_state.sample_grouping_radio, st.session_state.tax_level_radio),
				mime="text/csv")

			################### PIECHARTS ##############
			st.session_state.df_percentages = df_venn.apply(lambda x: ( x * 100 / (x.sum())).round(0))

			st.header('Pie charts')
			for index, i in enumerate(st.session_state.df_percentages.columns):

				df_torta = st.session_state.df_percentages[i][st.session_state.df_percentages[i] > 0]
				fig = px.pie(df_torta, values=df_torta.name, names=df_torta.index, title='Relative abundances %s - %s'%(st.session_state.tax_level_radio, i))
				
				st.plotly_chart(figure_or_data = fig)

				tabella_df = df_torta
				csv = convert_df(tabella_df)
				ste.download_button(
					label="Download percentage relative abundances table %s" %(i),
					data=csv,
					file_name= 'relative_abundances_perc_%s_%s_%s.csv' %(st.session_state.dashboard_name, st.session_state.tax_level_radio, i),
					mime="text/csv")
			# ############################################
			# st.header('Diagrammi di Venn')
			# # definisco il df dei taxa ripetuti il numero di volte percentuale in cui sono presenti in ciascun gruppo
			# st.session_state.dfs = []
			# for index, i in enumerate(st.session_state.df_percentages.columns):
			# 	df_repeated_taxa_perc = st.session_state.df_percentages.reindex(st.session_state.df_percentages.index.repeat(st.session_state.df_percentages['%s'%i]))
			# 	columns = [i, st.session_state.tax_level_radio]
			# 	df_to_be_appended = df_repeated_taxa_perc.iloc[:, index].reset_index(drop=False)
			# 	df_to_be_appended.columns = columns
			# 	df_to_be_appended = df_to_be_appended.iloc[:, 0]
			# 	st.session_state.dfs.append(df_to_be_appended)
			
			# st.info('E\' possibile scaricare le tabelle delle ripetizioni percentuali dei taxa per ogni raggruppamento di campioni \
			# 		per il livello tassonomico selezionato e caricarle come set differenti sul sito https://bioinfogp.cnb.csic.es/tools/venny/ \
			# 			o https://bioinformatics.psb.ugent.be/webtools/Venn/ \
			# 			per ottenere il __diagramma di Venn__ dei taxa condivisi e non fra gruppi.')
				
			# for i in st.session_state.dfs:
			# 	tabella_df = i
			# 	csv = convert_df(tabella_df)
			# 	ste.download_button(
			# 		label="Download tabella set taxa per diagrammi di Venn - %s" %(i.name),
			# 		data=csv,
			# 		file_name= 'set_taxa_Venn_%s_%s_%s.csv' %(st.session_state.dashboard_name, st.session_state.tax_level_radio, i.name),
			# 		mime="text/csv")
		else:
			
			st.warning('No metadata file to compute venn diagrams')

		# ###########################################################################
		# # sunburst prova
		# fig = px.sunburst(
		# 	st.session_state.final_df, 
		# 	path=taxa_levels, values=st.session_state.final_df.iloc[:, 8])
					
		# st.plotly_chart(figure_or_data=fig)

		# ############################################################################


	st.balloons()

with tab_alpha_div: # 4 METRICHE: shannon, simpson, pielou evenness, observed features
	
	side_placeholder4.info('***%s.*** Alpha diversity charts' %(step_n))

	st.header('Alpha Diversity')
	st.subheader('All samples')
	st.info('Visualization of groups\' comparisons based on the provided metadata: \
	 \n > At the end of the page you can download files of interactive visualizations of comparisons of alpha diversity metrics among groups. \
	 Go to https://view.qiime2.org/ and load a .qzv file to visualize the comparisons of alpha diversity metrics among samples\' metadata groups.')
	
	
	try:

		# Calcolo delle alfa diversita', importazione della OTU table trasposta come artifact qiime2 per calcolare le 4 metriche per ogni campione
		artifact_table = Artifact.import_data("FeatureTable[Frequency]", st.session_state.data_OTU_df.T)
		st.session_state.artifact_table = artifact_table
		alpha_divs = app_alpha_divs(_table=st.session_state.artifact_table)
		st.session_state.alpha_divs = alpha_divs
		
		secure_temp_dir_alpha_gr_sign = tempfile.mkdtemp(prefix="temp_", suffix="_alpha_gr_sign")
		
		cols_alpha_divs = st.columns(4)

		try:
# todo: ricalcoli le core metrics phylogenetic importando le sequenze rappresentative, la feature table e impostando la sampling depth come dai passaggi di analisi secondaria dei dati grezzi
# valuti se valga la pena importare i dati a questo punto oppure se spostare questa sezione nella analisi secondaria - magari no. Per ora questo passaggio funziona solo dopo aver corso la analisi secondaria
			alpha_group_sign = diversity.visualizers.alpha_group_significance(alpha_diversity=st.session_state.core_metr_phylo.faith_pd_vector, metadata=Metadata(st.session_state.data_meta_df))
			alpha_group_sign.visualization.export_data(secure_temp_dir_alpha_gr_sign)
			alpha_group_sign.visualization.save(secure_temp_dir_alpha_gr_sign+'/faith_pd.qzv')
		except:
			pass

		for index, i in enumerate(st.session_state.alpha_divs):
			
			tabella_df = pd.DataFrame(i.view(pd.Series))
			# Modifica nomi delle colonne
			tabella_df.columns = map(lambda x: ' '.join(x.split('_')).upper(), tabella_df.columns)

			# cols_alpha_divs[index].table(tabella_df)
			if 'data_meta_df' in st.session_state.keys():

				tabella_df = pd.DataFrame(tabella_df)
				tabella_df.columns = map(lambda x: ' '.join(x.split('_')).upper(), tabella_df.columns)
				tabella_df = tabella_df.merge(st.session_state.data_meta_df, left_index=True, right_index=True)
			else:
				
				st.warning('No metadata file.')
			if st.session_state.sample_grouping_radio != 'All samples':
				
				# due tipologie di grafici a box and whiskers: uno generico delle metriche per tutti i campioni insieme (prima riga) che mostra anche la media e la dev st
				# uno a box raggruppati per gruppi di campioni in base al raggruppamento dei campioni selezionato nel menu radio a lato (seconda riga)
				fig_boxplot = go.Figure()
				fig_boxplot.add_trace(go.Box(
					y=tabella_df.iloc[:,0],
					name = tabella_df.columns[0],
					#marker_color='darkblue',
					boxmean='sd'))

				fig_boxplot1 = px.box(
					tabella_df, 
					y=tabella_df.columns[0],
					color=st.session_state.sample_grouping_radio)
				cols_alpha_divs[index].subheader('%s'%(tabella_df.columns[0]))
				cols_alpha_divs[index].plotly_chart(fig_boxplot, use_container_width=True, config=config)
				cols_alpha_divs[index].subheader('%s - %s' %(tabella_df.columns[0], st.session_state.sample_grouping_radio))
				cols_alpha_divs[index].plotly_chart(fig_boxplot1, use_container_width=True, config=config)
				
				data_Metadata = st.session_state.data_meta_df
				data_Metadata.index.name = '#SampleID'
				st.session_state.data_Metadata = Metadata(data_Metadata)
				alpha_group_sign = diversity.visualizers.alpha_group_significance(alpha_diversity=i, metadata=st.session_state.data_Metadata)
				alpha_group_sign.visualization.export_data(secure_temp_dir_alpha_gr_sign)
				if index == 0:
					metric_name = 'observed_features'
				elif index == 1:
					metric_name = 'shannon_entropy'
				elif index == 2:
					metric_name = 'pielou_evenness'
				elif index == 3:
					metric_name = 'simpson'
				alpha_group_sign.visualization.save(secure_temp_dir_alpha_gr_sign+'/%s.qzv'%(metric_name))

				with open(secure_temp_dir_alpha_gr_sign+'/index.html', 'r', encoding='utf-8') as HtmlFile:
					source_code = HtmlFile.read()
				
				#components.html(source_code, scrolling=True)
				cols_alpha_divs[index].write(pd.DataFrame(pd.read_html(source_code, decimal='.', index_col=0)[0]))
				
				cols_alpha_divs[index].subheader('Kruskal-Wallis (pairwise)')
				try:

					tabella_df = pd.read_csv(secure_temp_dir_alpha_gr_sign+'/kruskal-wallis-pairwise-%s.csv'%(st.session_state.sample_grouping_radio))
					cols_alpha_divs[index].table(tabella_df)
				except:
					
					cols_alpha_divs[index].warning('Single pool of samples. It\'s impossible to compare among pair groups of samples.')
			
				try:

					with open(secure_temp_dir_alpha_gr_sign+"/%s.qzv"%(metric_name), 'rb') as f:
						ste.download_button(
							label="Download alpha diversity metric %s groups comparison .qzv" %(metric_name),
							data=f,
							file_name="alfa_div_%s_groups_comparison.qzv" %(metric_name),
							mime="application/qzv")
				except Exception as e:
					
					st.exception(e)
					st.warning('No metadata file. There are no groups to be compared.')
					pass

			else:
				
				# se Tutti i campioni: una sola tipologia di grafici a box and whiskers di tutti i campioni insieme, che mostra anche la media e la dev st (una riga)
				fig_boxplot = go.Figure()
				fig_boxplot.add_trace(go.Box(
					y=tabella_df.iloc[:,0],
					name = tabella_df.columns[0],
					#marker_color='darkblue',
					boxmean='sd'))
				
				cols_alpha_divs[index].subheader(st.session_state.sample_grouping_radio)
				cols_alpha_divs[index].plotly_chart(fig_boxplot, use_container_width=True, config=config)
					
	except Exception as e:
		
		st.warning('All samples selected. It\'s impossible to compare among groups of samples.')
		st.exception(e)
		
	side_placeholder4.success('***%s.*** Alpha diversity charts' %(step_n))
	step_n += 1

	
	st.balloons()

with tab_beta_div:

	side_placeholder5.info('***%s.*** Beta diversity charts' %(step_n))

	st.header('Beta Diversity')
	
	st.info('If a samples\' grouping is selected from the sidebar menu, at the end of the page you can download files of interactive \
	 visualizations of comparisons of beta diversity metrics among groups. \
	 \n> Go to https://view.qiime2.org/ and load a .qzv file to visualize the comparisons of beta diversity metrics among samples\' metadata groups.')
	
	try:

		if st.session_state.data_meta_df is not None:
			
				
			st.session_state['feat_table'] = Artifact.import_data(type="FeatureTable[Frequency]", view=st.session_state.data_OTU_df.T)

			
			beta_result = beta(
				table = st.session_state['feat_table'],
				metric = 'braycurtis',
				pseudocount = 1,
				n_jobs = 'auto')
			st.session_state.bc = beta_result.distance_matrix
			
			st.session_state.bc_gr_sig = beta_group_significance(
				distance_matrix = st.session_state.bc,
				metadata = Metadata(st.session_state.data_meta_df).get_column(st.session_state.sample_grouping_radio),
				method = 'permanova',
				permutations = 999)
			
			st.session_state.bc_pcoa = pcoa(
				distance_matrix = st.session_state.bc,
				number_of_dimensions=3
				)
			
			st.session_state.bc_emperor = plot(
				pcoa = st.session_state.bc_pcoa.pcoa,
				metadata=Metadata(st.session_state.data_meta_df),
				ignore_missing_samples=True
				)
			

			secure_temp_dir_beta_gr_sig = tempfile.mkdtemp(prefix="temp_", suffix="_beta_gr_sig")
		
			st.session_state.bc_gr_sig.visualization.save(secure_temp_dir_beta_gr_sig+'/bray_curtis.qzv')

			st.session_state.bc_emperor.visualization.save(secure_temp_dir_beta_gr_sig+'/bray_curtis_emperor.qzv')
			
			beta_result = beta(
				table = st.session_state['feat_table'],
				metric = 'jaccard',
				pseudocount = 1,
				n_jobs = 'auto')
			st.session_state.j = beta_result.distance_matrix
			
			st.session_state.j_gr_sig = beta_group_significance(
				distance_matrix = st.session_state.j,
				metadata = Metadata(st.session_state.data_meta_df).get_column(st.session_state.sample_grouping_radio),
				method = 'permanova',
				permutations = 999)
			
			st.session_state.j_pcoa = pcoa(
				distance_matrix = st.session_state.j,
				number_of_dimensions=3
				)
			
			st.session_state.j_emperor = plot(
				pcoa = st.session_state.j_pcoa.pcoa,
				metadata=Metadata(st.session_state.data_meta_df),
				ignore_missing_samples=True
				)
			
			st.session_state.j_gr_sig.visualization.save(secure_temp_dir_beta_gr_sig+'/jaccard.qzv')
			
			st.session_state.j_emperor.visualization.save(secure_temp_dir_beta_gr_sig+'/jaccard_emperor.qzv')

			with open(secure_temp_dir_beta_gr_sig+"/bray_curtis.qzv", 'rb') as f:
				ste.download_button(
					label="Download comparison beta diversity metric Bray Curtis among groups .qzv",
					data=f,
					file_name="beta_div_bray_curtis_groups_comparison.qzv",
					mime="application/qzv")
			with open(secure_temp_dir_beta_gr_sig+"/bray_curtis_emperor.qzv", 'rb') as f:
				ste.download_button(
					label="Download PCoA Bray Curtis .qzv",
					data=f,
					file_name="beta_div_bray_curtis_PCoA.qzv",
					mime="application/qzv")
			with open(secure_temp_dir_beta_gr_sig+"/jaccard.qzv", 'rb') as f:
				ste.download_button(
					label="Download comparison beta diversity metric Jaccard among groups .qzv",
					data=f,
					file_name="beta_div_jaccard_groups_comparison.qzv",
					mime="application/qzv")
			with open(secure_temp_dir_beta_gr_sig+"/jaccard_emperor.qzv", 'rb') as f:
				ste.download_button(
					label="Download PCoA Jaccard .qzv",
					data=f,
					file_name="beta_div_jaccard_PCoA.qzv",
					mime="application/qzv")
		else:
			st.warning('No metadata file. It\'s impossible to compare among groups of samples.')

	except:

		st.warning('All samples selected or single pool grouping of samples. It\'s impossible to compare among groups of samples.')
		pass

	side_placeholder5.success('***%s.*** Beta diversity charts' %(step_n))
	step_n += 1

	
	st.balloons()

try:
	shutil.rmtree(secure_temp_dir_alpha_gr_sign)
except FileNotFoundError as e:
	pass
	
except NameError as e:
	pass



try:
	shutil.rmtree(secure_temp_dir_beta_gr_sig)
except FileNotFoundError as e:
	pass
	
except NameError as e:
	pass
################################################################### PIESPARROW ##############################################################################################################

# PIE CHARTS
import piesparrow as ps

st.session_state.sequenced_samples = st.session_state.data_OTU_df.columns.tolist()

secure_temp_dir_dashboard_bars = tempfile.mkdtemp(prefix="temp_", suffix="_dashboard_barre")
secure_temp_dir_dashboard_donuts = tempfile.mkdtemp(prefix="temp_", suffix="_dashboard_torte")

ps.init(
	filename=secure_temp_dir_dashboard_donuts+'/Dashboard-microbiome-PieCharts-%s-%s'%(st.session_state.dashboard_name, st.session_state.sample_grouping_radio), 
	title='Dashboard pies - pieSparrow', 
	basetheme=ps.light, 
	charttheme=ps.sparrow_light
)

ps.row(
	ps.colxl(align='left', type='card', content = 
		ps.h1(ps.bold('Bacterial Microbiome Dashboard'))
	
	+	ps.p('This dashboard shows data analyses results \
		from the Microbiome App through interactive visualizations \
			- developed in python with the library piesparrow.')
	)
)
ps.row(
	ps.colsm(content=ps.p(''))
+	ps.colmd(type='card', align='center',content=
		ps.h2('%s' %(st.session_state.dashboard_name))
	+	ps.h3('%s' %(st.session_state.sample_grouping_radio))
	)
)

df = st.session_state.final_df[[st.session_state.tax_level_radio]+ st.session_state.sequenced_samples].groupby(st.session_state.tax_level_radio).sum()
df = df.T.reset_index(level=0)
df_table = df

	
try:

	df = pd.merge(df, st.session_state.data_meta_df, how='left', left_on='index', right_index=True)
	df = df.groupby(st.session_state.sample_grouping_radio).sum().reset_index(level=0)
	# numero di colonne dei metadati numeriche
	n_num_cols = len(st.session_state.data_meta_df.select_dtypes(include=np.number).columns.tolist())
	# seleziono solo le colonne relative ai taxa escludendo le colonne di metadati numeriche mergiate alla fine del dataframe
	df_table = df.iloc[:, :-n_num_cols]
		

		
except Exception as e:
	df_table = df
	# st.exception(e)
	pass



ps.row(

	ps.colxl(type='card', align='left', content = 
		ps.h1('%s' %(st.session_state.tax_level_radio))
	+	ps.table(df=df_table)
	)
)



for j, s_group in enumerate(df.iterrows()):

	try:
		gr = df.loc[:, st.session_state.sample_grouping_radio][j]
		df_pie= pd.DataFrame(df.iloc[j,:]).T
		df_pie= df_pie.set_index(st.session_state.sample_grouping_radio).dropna(how='all', axis=1)
		
	except Exception as e: # 'Tutti i campioni'
		# st.exception(e)
		gr = s_group[1][0]
		df_pie= pd.DataFrame(df.iloc[j,:]).T.set_index('index').dropna(how='all', axis=1)
		
	

	
	ps.row(

		ps.colxl(type='card', align='center', content =
			ps.h2(ps.bold('%s' %(gr)))
		)
	+	ps.colxl(align='center', content =
			ps.donut(
				title = '%s, %s'%(st.session_state.tax_level_radio, gr),
				df = df_pie,
				columns = df.columns[1:], 
			)
		)
	)
	
		
#################################################################################################
# PIESPARROW
# BARCHARTS


ps.init(
	filename=secure_temp_dir_dashboard_bars+'/Dashboard-microbiome-BarCharts-%s-%s'%(st.session_state.dashboard_name, st.session_state.sample_grouping_radio), 
	title='Dashboard bars - pieSparrow', 
	basetheme=ps.light, 
	charttheme=ps.sparrow_light
)

ps.row(
	ps.colxl(align='left', type='card', content = 
		ps.h1(ps.bold('Bacterial Microbiome Dashboard'))
	
	+	ps.p('This dashboard shows data analyses results \
		from the Microbiome App through interactive visualizations \
			- developed in python with the library piesparrow.')
	)
)

ps.row(
	ps.colsm(content=ps.p(''))
+	ps.colmd(type='card', align='center',content=
		ps.h2('%s' %(st.session_state.dashboard_name))
	+	ps.h3('%s' %(st.session_state.sample_grouping_radio))
	)
)



df = st.session_state.final_df[[st.session_state.tax_level_radio]+ st.session_state.sequenced_samples].groupby(st.session_state.tax_level_radio).sum()
df = df.T.reset_index(level=0)
df_table = df

try:

	df = pd.merge(df, st.session_state.data_meta_df, how='left', left_on='index', right_index=True)
	df = df.groupby(st.session_state.sample_grouping_radio).sum().reset_index(level=0)
	df_table = df

	df_barchart = df.T[1:]
	df_barchart_hdr = df.T.iloc[0,:]
	df_barchart.columns = df_barchart_hdr
	df_barchart.index.name = st.session_state.tax_level_radio
	df_barchart = df_barchart.reset_index(level=0)
	df_barchart.index = range(0, len(df_barchart.index))
	df_table = df_barchart
	
except Exception as e: # if sample_grouping_radio = 'Tutti i campioni'

	# st.exception(e)
	df_barchart = df
	df_table = df_barchart
	pass


ps.row(

	ps.colxl(type='card', align='center', content = 

		ps.h1(ps.bold('%s' %(st.session_state.tax_level_radio)))
	+
	
		ps.table(df=df_table)
	)
)
ps.row(
	ps.chart(
			type='bar',
			title='barchart', 
			df=df_barchart, 
			columns=df_barchart.columns, 
			xcolumn=[st.session_state.tax_level_radio], 
			height=700,
		
	)
)

with dashboard_info_dwnld_plchldr:
	st.info('You can download two HTML dashboards of taxa relative abundances.\
		 Interactive visualizations as bar charts and pie charts.')
with dashboard_barre_dwnld_plchldr:
	
	with open(secure_temp_dir_dashboard_bars+'/Dashboard-microbiome-BarCharts-%s-%s.html'%(st.session_state.dashboard_name, st.session_state.sample_grouping_radio), 'rb') as f:
		dwnld_bttn_bars = ste.download_button(
			label="Download Dashboard HTML Bar charts",
			data=f,
			file_name="Dashboard_Bar_charts__%s_%s.html" %(st.session_state.dashboard_name, st.session_state.sample_grouping_radio),
			mime="text/html")

with dashboard_torte_dwnld_plchldr:

	with open(secure_temp_dir_dashboard_donuts+'/Dashboard-microbiome-PieCharts-%s-%s.html'%(st.session_state.dashboard_name, st.session_state.sample_grouping_radio), 'rb') as f:
		dwnld_bttn_donuts = ste.download_button(
			label="Download Dashboard HTML Pie charts",
			data=f,
			file_name="Dashboard_Pie_charts__%s_%s.html" %(st.session_state.dashboard_name, st.session_state.sample_grouping_radio),
			mime="text/html")
		
# if dwnld_bttn_donuts:
try:
	shutil.rmtree(secure_temp_dir_dashboard_donuts)
except FileNotFoundError as e:
	st.exception(e)
	
except NameError as e:
	st.exception(e)

# if dwnld_bttn_bars:
try:
	shutil.rmtree(secure_temp_dir_dashboard_bars)
except FileNotFoundError as e:
	st.exception(e)
	
except NameError as e:
	st.exception(e)

