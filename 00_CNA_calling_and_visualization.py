################################################################################
# analysis.py
#
# Python analysis pipeline supporting the single-cell aneuploidy study.
#
# This script covers:
#   1.  Chromosome-arm CNA calling via ASCETS
#   2.  Reordering ASCETS output files into canonical arm order
#   3.  Building per-sample enrichment q-value matrices from DGE results
#   4.  Visualizing enrichment matrices as clustered heatmaps
#   5.  Visualizing copy number profiles as cell × arm / cell × gene heatmaps,
#       with cells grouped by cell type, AS group, or CN subclone
#
# NOTE: Input data (CNA matrices, cell metadata, enrichment CSVs) are produced
# by the R scripts (01_main_analysis.R, 02_recurrent_aneuploidies.R).
# The gene expression data are publicly available; see the paper for accessions.
#
# All file paths should be updated to your local directory structure.
################################################################################


# ==============================================================================
# 0.  Imports & shared constants
# ==============================================================================

from multiprocessing import Pool
import os
import glob
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import scipy  # noqa: F401 — available for downstream use

sys.setrecursionlimit(10000)

# Canonical chromosome-arm sort order (autosomes + sex chromosomes, p before q)
ARM_ORDER = [
    '1p',  '1q',  '2p',  '2q',  '3p',  '3q',  '4p',  '4q',
    '5p',  '5q',  '6p',  '6q',  '7p',  '7q',  '8p',  '8q',
    '9p',  '9q',  '10p', '10q', '11p', '11q', '12p', '12q',
    '13p', '13q', '14p', '14q', '15p', '15q', '16p', '16q',
    '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q',
    '21p', '21q', '22p', '22q', 'Xp',  'Xq',  'Yp',  'Yq',
]

# CNA value → integer mapping used in arm-level heatmaps
CNA_VALUE_MAP = {"AMP": 1, "DEL": -1, "NEUTRAL": 0, "NC": 0, "LOWCOV": 0}


# ==============================================================================
# 1.  Chromosome-Arm CNA Calling via ASCETS
# ==============================================================================
# Runs the ASCETS R script on every .seg file found under the study tree.
# Each .seg file represents one sample; ASCETS outputs arm-level call files
# and weighted-average segment-mean files.
#
# Input:  .seg files under <study_base>/<study>/<sample>_CNA/
# Output: per-sample ASCETS output directory next to each .seg file
#
# ascets_script  — absolute path to run_ascets.R
# arm_coords     — absolute path to genomic_arm_coordinates_hg38.txt
# n_workers      — number of parallel processes (default 30)

def run_ascets(study_base, ascets_script, arm_coords, n_workers=30):
    """Discover all .seg files and run ASCETS on each in parallel."""
    seg_files = glob.glob(
        os.path.join(study_base, '*', '*_CNA', '*.seg')
    )
    print(f"Running ASCETS on {len(seg_files)} .seg files")
    with Pool(n_workers) as pool:
        pool.starmap(
            _run_ascets_one_file,
            [(f, ascets_script, arm_coords) for f in seg_files]
        )


def _run_ascets_one_file(filename, ascets_script, arm_coords):
    """
    Run ASCETS on a single .seg file.

    Input:
        filename     — full path to the .seg file
        ascets_script — path to run_ascets.R
        arm_coords   — path to genomic_arm_coordinates_hg38.txt

    Output directory: <seg_stem>/
    """
    sample_id  = os.path.basename(filename).split('.')[0]
    seg_stem   = filename.replace('.seg', '')
    output_dir = seg_stem
    output_pfx = os.path.join(output_dir, sample_id)

    os.makedirs(output_dir, exist_ok=True)

    cmd = (
        f'Rscript "{ascets_script}"'
        f' -i {filename}'
        f' -c "{arm_coords}"'
        f' -m 0.0'
        f' -o {output_pfx}'
    )
    os.system(cmd)
    print(f'{sample_id} completed!')


# ==============================================================================
# 2.  Reorder ASCETS Output into Canonical Arm Order
# ==============================================================================
# After ASCETS runs, both the weighted-average segment-means file and the
# arm-level calls file are reordered to match ARM_ORDER and saved with an
# _ordered suffix.  This ensures consistent column ordering across all samples
# for downstream merging and visualisation.
#
# Input:  ASCETS output directories matching <study_base>/*/*_CNA/*_defaultParamsBoc0.1/
# Output: *_weighted_average_segmeans_ordered.txt and *_arm_level_calls_ordered.txt

def reorder_ascets_outputs(study_base, n_workers=1):
    """Reorder all ASCETS output files in parallel."""
    dirs = glob.glob(
        os.path.join(study_base, '*', '*_CNA', '*', '')
    )
    print(f"Reordering outputs in {len(dirs)} directories")
    with Pool(n_workers) as pool:
        pool.map(_reorder_one_dir, dirs)


def _reorder_one_dir(dirpath):
    """
    Reorder arm columns in weighted-average and arm-level call files for one sample.

    Input:  dirpath — path to a single ASCETS output directory
    Output: _ordered.txt files written alongside originals
    """
    for pattern in [
        '*_weighted_average_segmeans.txt',
        '*_arm_level_calls.txt',
    ]:
        for filename in glob.glob(os.path.join(dirpath, pattern)):
            parts       = filename.split('.')
            ordered_path = parts[0] + '_ordered.' + parts[-1]
            df = pd.read_csv(filename, sep='\t', index_col='sample')
            df = df.reindex(ARM_ORDER, axis=1).dropna(how='all', axis=1)
            df.to_csv(ordered_path, sep='\t')
            print(f"Saved {os.path.basename(ordered_path)}")


# ==============================================================================
# 3.  Build Enrichment Q-Value Matrix from DGE Results
# ==============================================================================
# Reads the per-sample enricher output CSVs (from 01_main_analysis.R) and
# assembles a pathways × samples matrix of q-values (significant at q < 0.25).
# Separate matrices are built for positively and negatively enriched pathways.
#
# The same function is used for both:
#   - Group-based DGE results
#   - Continuous lmFit results
# The difference is only in the glob pattern passed to locate the enricher files.
#
# Input:
#   results_base    — root directory of the DGE/enrichment output tree
#   db              — pathway database prefix, e.g. "HALLMARK", "REACTOME"
#   direction       — "positive" (up in highAS) or "negative" (down in highAS)
#   file_pattern    — glob pattern for enricher CSV files within results_base
#   excluded_studies — list of study names to skip (e.g. low-quality samples)
#   output_file     — path to write the resulting q-value matrix CSV
#
# Output: pathways × samples CSV; values are q-values for significant results

def build_enrichment_matrix(
    results_base,
    db,
    direction,
    file_pattern,
    excluded_studies=None,
    output_file=None,
):
    """
    Assemble a pathways × samples enrichment q-value matrix.
    """
    excluded_studies = excluded_studies or []

    all_files = glob.glob(
        os.path.join(results_base, file_pattern), recursive=True
    )
    db_files = [
        f for f in all_files
        if os.path.basename(f).split('_')[0] == db
    ]
    print(f"Found {len(db_files)} {db} {direction} enricher files")

    matrix = pd.DataFrame()
    for filepath in db_files:
        # Extract study and sample from path; skip excluded studies
        parts = filepath.split('/')
        study  = parts[-5]
        sample = parts[-4]
        if study in excluded_studies:
            print(f"Skipping {study} — excluded")
            continue

        sample_id = f"{study}/{sample}/{parts[-3]}"

        df = pd.read_csv(filepath)
        # Harmonise q-value column: some clusterProfiler versions use
        # 'qvalue', others use 'p.adjust' when BH correction is applied.
        df['qvalues'] = df['qvalue'].fillna(df.pop('p.adjust'))
        df = df.sort_values('qvalues').set_index('Description')
        df = df[df['qvalues'] < 0.25]

        if len(df) > 0:
            col = pd.DataFrame({sample_id: df['qvalues']})
            matrix = pd.concat([matrix, col], axis=1)

    if output_file:
        matrix.to_csv(output_file)
        print(f"Saved enrichment matrix → {output_file}")

    return matrix


# ==============================================================================
# 4.  Enrichment Heatmaps (Supp Fig 3A+B)
# ==============================================================================
# Merges positive and negative enrichment matrices into a single signed
# log10(1/q) matrix, then renders a clustered heatmap (seaborn clustermap).
#
# For each pathway and sample:
#   - If enriched in the positive direction: log10(1/q)  (positive value)
#   - If enriched in the negative direction: -log10(1/q) (negative value)
#   - If enriched in both directions for the same sample: set to NaN
#     (ambiguous cases are excluded)
#
# This function is shared between:
#   - Group-based DGE (highAS vs lowAS) heatmap
#   - Continuous lmFit-based heatmap
# The only difference is the input matrices passed in.
#
# Inputs:
#   pos_matrix      — pathways × samples q-value matrix, positive direction
#   neg_matrix      — pathways × samples q-value matrix, negative direction
#   sample_map      — sampleFilesMap data.frame (CSV) for cancer-type colours
#   samples_filter  — optional data.frame with a 'col_name' column used to
#                     restrict to a validated sample subset
#   nan_threshold   — drop pathways present in fewer than this fraction of
#                     samples (default 0.15)
#   vmin, vmax      — heatmap colour scale limits (default -10, 10)
#   figsize         — figure size tuple (default (30, 25))
#
# Output: matplotlib figure (call plt.show() or plt.savefig() after)

def plot_enrichment_heatmap(
    pos_matrix,
    neg_matrix,
    sample_map_path,
    samples_filter=None,
    nan_threshold=0.15,
    vmin=-10,
    vmax=10,
    figsize=(30, 25),
):
    sample_map = pd.read_csv(sample_map_path)

    # Convert q-values to signed log10 scores
    pos_log = np.log10(1 / pos_matrix)
    neg_log = np.log10(1 / neg_matrix) * -1

    # Remove ambiguous cases: pathway enriched in both directions for same sample
    for pathway in pos_log.index:
        if pathway in neg_log.index:
            common_samples = set(
                pos_log.loc[pathway].dropna().index
            ).intersection(neg_log.loc[pathway].dropna().index)
            for s in common_samples:
                pos_log.loc[pathway, s] = np.nan
                neg_log.loc[pathway, s] = np.nan

    # Merge: keep the non-NaN value per (pathway, sample)
    combined = pd.concat([pos_log, neg_log])
    merged   = pd.DataFrame()
    for col in combined.columns:
        series = combined.groupby(level=0)[col].apply(
            lambda x: x.dropna().iloc[0] if x.dropna().size > 0 else np.nan
        )
        merged = pd.concat([merged, series.rename(col)], axis=1)

    # Optionally restrict to validated samples
    if samples_filter is not None:
        valid_cols = samples_filter['col_name'][
            samples_filter['col_name'].isin(merged.columns)
        ].unique()
        merged = merged[valid_cols]

    # Cancer-type colour bar
    cancer_types = [
        sample_map.loc[
            (sample_map['Sample'] == c.split('/')[1]) &
            (sample_map['Study']  == c.split('/')[0]),
            'Cancer_type'
        ].values[0]
        for c in merged.columns
    ]
    palette     = sns.color_palette('deep', len(np.unique(cancer_types)))
    colour_lut  = dict(zip(np.unique(cancer_types), palette))
    col_colours = pd.Series(cancer_types).map(colour_lut).values

    # Drop sparse pathways and replace inf / NaN for clustering
    df_plot = (
        merged
        .dropna(thresh=int(len(merged.columns) * nan_threshold))
        .replace([np.inf, -np.inf], np.nan)
        .fillna(1e-8)
    )

    g = sns.clustermap(
        df_plot,
        cmap='bwr',
        vmin=vmin, vmax=vmax,
        yticklabels=1, xticklabels=0,
        metric='euclidean',
        figsize=figsize,
        dendrogram_ratio=(0.05, 0.1),
        cbar_kws={'orientation': 'horizontal'},
    )

    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    g.ax_cbar.set_xlabel('log10(1/qvalue) & direction', size=55)
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)
    bbox = g.ax_heatmap.get_position()
    g.gs.update(left=0.04, right=bbox.x1)
    g.ax_cbar.set_position((0.709, 0.98, 0.215, 0.03))
    g.ax_cbar.tick_params(labelsize=50)
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_ymajorticklabels(), fontsize=60
    )

    return g


# ==============================================================================
# 5.  Copy Number Profile Heatmaps
# ==============================================================================
# Three heatmap functions visualise per-cell copy number profiles for a single
# sample.  All share the same core layout logic: cells are rows, genomic
# positions are columns, and cells are grouped by a biological label with
# horizontal dividing lines between groups.
#
# The same _add_group_dividers() helper is used by all three functions.
#
# Shared inputs across all three functions:
#   paper      — study identifier (for plot title and file naming)
#   sample     — sample identifier
#   cells_path — CSV with columns: cell_name, sample, cell_type
#
# ==============================================================================

def _add_group_dividers(ax, sample_cells_filtered, type_col='type'):
    """
    Add horizontal lines between cell-type groups on a heatmap axes object,
    and label each group at its midpoint on the y-axis.

    Input:
        ax                     — matplotlib axes (heatmap)
        sample_cells_filtered  — data.frame sorted in row order of the heatmap,
                                 with a column identifying the group
        type_col               — name of the grouping column (default 'type')
    """
    sizes     = sample_cells_filtered.groupby(type_col, sort=False).size()
    cumsum    = sizes.cumsum()
    midpoints = (cumsum - sizes // 2).values
    labels    = sizes.index.tolist()

    for boundary in cumsum.values[:-1]:
        ax.axhline(boundary, color='black', linewidth=1.5)

    ax.set_yticks(midpoints)
    ax.set_yticklabels(labels, rotation=0, fontsize=15)


# ==============================================================================
# 5a.  Gene-Level CNA Heatmap (Supp Fig 8E)
# ==============================================================================
# Renders a cells × segments heatmap using the raw log2-ratio values from the
# .seg file.  Cells are grouped as Malignant vs Non-malignant.
#
# This provides a higher-resolution CNA view than the arm-level heatmap.
#
# Input:
#   seg_path   — path to the .seg file (columns: chrom, segment_start,
#                segment_end, log2ratio; index: cell_name)
#   cells_path — path to cells CSV
#
# Output: matplotlib figure

def plot_gene_level_CN_heatmap(seg_path, cells_path, paper, sample):
    """
    Plot a gene-level (segment-resolution) CNA heatmap for one sample.

    Cells are split into Malignant and Non-malignant groups.
    Columns are genomic segments, labelled by chromosome number on the x-axis.
    Colour scale: blue (loss) → white (neutral) → red (gain), log2 ratio ±1.
    """
    seg = pd.read_csv(seg_path, sep='\t', index_col='sample')
    # Restrict to autosomal chromosomes
    seg = seg[seg['chrom'].apply(lambda x: str(x).isdigit())]

    cells_df = pd.read_csv(cells_path, encoding='unicode_escape')
    sample_cells = cells_df[cells_df['sample'].astype(str) == sample].copy()
    sample_cells['cell_name'] = sample_cells['cell_name'].astype(str)
    sample_cells['type'] = sample_cells['cell_type'].where(
        sample_cells['cell_type'].str.lower() == 'malignant', 'Non malignant'
    )

    call = seg[seg.index.astype(str).isin(sample_cells['cell_name'])].copy()
    call = call.sort_values(['chrom', 'segment_start'])
    call['segment_id'] = (
        call['chrom'].astype(str) + ':' +
        call['segment_start'].astype(str) + '-' +
        call['segment_end'].astype(str)
    )
    call['cell_id'] = call.index
    call = call.drop_duplicates(subset=['cell_id', 'segment_id'])

    matrix = call.pivot(index='cell_id', columns='segment_id', values='log2ratio')
    matrix = matrix[sorted(
        matrix.columns,
        key=lambda x: (int(x.split(':')[0]), int(x.split(':')[1].split('-')[0]))
    )]
    matrix.index = matrix.index.astype(str)

    # Sort cells: Non-malignant first (alphabetical), then Malignant
    sample_cells_filt = (
        sample_cells[sample_cells['cell_name'].isin(matrix.index)]
        .sort_values('type', ascending=False)
    )
    matrix = matrix.reindex(sample_cells_filt['cell_name'])

    # Build chromosome x-tick positions (centred per chromosome)
    chrom_labels = [s.split(':')[0] for s in matrix.columns]
    xtick_positions, xtick_labels = [], []
    last_chrom, start_idx = None, 0
    for i, chrom in enumerate(chrom_labels):
        if chrom != last_chrom:
            if last_chrom is not None:
                xtick_positions.append((start_idx + i - 1) // 2)
                xtick_labels.append(last_chrom)
            start_idx, last_chrom = i, chrom
    xtick_positions.append((start_idx + len(chrom_labels) - 1) // 2)
    xtick_labels.append(last_chrom)

    # Diverging colour map: strong loss (dark blue) → neutral (white) → strong gain (dark red)
    # A ±0.04 noise gate around zero is rendered white to suppress low-level noise.
    colours = [
        (0.00, '#00008B'),   # log2 = -1.0 : dark navy  (strong loss)
        (0.15, '#4169E1'),   # log2 = -0.6 : royal blue (moderate loss)
        (0.45, '#FFFFFF'),   # log2 = -0.04: start noise gate
        (0.50, '#FFFFFF'),   # log2 =  0.0 : neutral
        (0.55, '#FFFFFF'),   # log2 = +0.04: end noise gate
        (0.85, '#FA8072'),   # log2 = +0.6 : salmon     (moderate gain)
        (1.00, '#B22222'),   # log2 = +1.0 : firebrick  (strong gain)
    ]
    cmap = mcolors.LinearSegmentedColormap.from_list('Goldilocks', colours, N=256)

    fig, ax = plt.subplots(figsize=(8, 5))
    hm = sns.heatmap(
        matrix, ax=ax,
        cmap=cmap, vmin=-1, vmax=1,
        yticklabels=False, xticklabels=False,
        rasterized=True,
        cbar_kws={'ticks': [-1, -0.5, 0, 0.5, 1],
                  'shrink': 0.3, 'anchor': (0, 1),
                  'label': 'inferred CNA\n(log2 ratio)'},
    )

    _add_group_dividers(ax, sample_cells_filt)

    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=13, fontweight='light')
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('inferred CNA\n(log2 ratio)', size=14)
    ax.set_ylabel('Cells — Gene-Level CN', fontsize=20)
    ax.set_xlabel('Chromosome', fontsize=20)
    ax.set_title(f'{paper}  {sample}', fontsize=20)

    return fig


# ==============================================================================
# 5b.  Arm-Level CNA Heatmap — Malignant vs Non-Malignant (Supp Fig 8F)
# ==============================================================================
# Renders a cells × chromosome-arm heatmap using the ASCETS arm-level call
# file.  CNA values (AMP/DEL/NEUTRAL) are mapped to +1 / -1 / 0.
# Cells are split into Malignant and Non-malignant groups.
#
# Input:
#   calls_path  — path to *_arm_level_calls_ordered.txt
#   cells_path  — path to cells CSV
#   as_path     — path to per-cell AS CSV (columns: cell_name, pred, AS)
#
# Output: matplotlib figure

def plot_arm_level_CN_heatmap(calls_path, cells_path, as_path, paper, sample):
    """
    Plot an arm-level CNA heatmap for one sample.

    CNA values are: AMP=+1, DEL=-1, NEUTRAL/NC/LOWCOV=0.
    Cells are split into Malignant (sorted by pred) and Non-malignant groups.
    Sex chromosomes (Xp, Xq, Yp, Yq) are excluded.
    """
    calls = pd.read_csv(calls_path, sep='\t', index_col='sample')
    calls = calls.drop(['Xp', 'Xq', 'Yp', 'Yq'], axis=1, errors='ignore')
    calls = calls.replace(CNA_VALUE_MAP)

    cells_df  = pd.read_csv(cells_path, encoding='unicode_escape')
    sample_cells = cells_df[cells_df['sample'].astype(str) == sample].copy()

    as_df = pd.read_csv(as_path, encoding='unicode_escape')
    as_df.rename(columns={as_df.columns[0]: 'cell_name'}, inplace=True)
    as_df['cell_name'] = as_df['cell_name'].astype(str)

    sample_cells = sample_cells.merge(
        as_df[['cell_name', 'pred', 'AS']], on='cell_name', how='left'
    )
    sample_cells['type'] = sample_cells['cell_type'].where(
        sample_cells['cell_type'].str.lower() == 'malignant', 'Non malignant'
    )

    sample_cells_filt = (
        sample_cells[sample_cells['cell_name'].isin(calls.index)]
        .sort_values('type', ascending=False)
    )
    calls = calls.reindex(sample_cells_filt['cell_name'])

    cmap = mcolors.LinearSegmentedColormap.from_list(
        'cna_3class', ['blue', 'white', 'red'], N=3
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    hm = sns.heatmap(
        calls, ax=ax,
        cmap=cmap, vmin=-1, vmax=1,
        yticklabels=False,
        cbar_kws={'ticks': [-1, 0, 1], 'shrink': 0.3,
                  'anchor': (0, 1), 'label': 'Arm CN'},
    )

    _add_group_dividers(ax, sample_cells_filt)

    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('Arm CN', size=14)
    ax.set_ylabel('Cells CN', fontsize=20)
    ax.set_xlabel('Chromosome Arm', fontsize=20)
    ax.set_title(f'{paper}  {sample}', fontsize=20)
    plt.xticks(rotation=90, fontweight='light', fontsize=11)

    return fig


# ==============================================================================
# 5c.  Arm-Level CNA Heatmap — AS Groups or CN Subclones
# ==============================================================================
# Same layout as 5b, but cells are grouped by either:
#   (a) AS group (lowAS / midAS / highAS / NonMalignant) — for AS-stratified figures
#   (b) CN subclone identity (subclone_1, subclone_2, … / NonMalignant) — for
#       subclone figures
#
# The grouping is controlled by the `group_by` argument:
#   "as_group"  — reads per-cell pred column from as_path
#   "subclone"  — reads MergedCluster column from subclones_path
#
# Input:
#   calls_path     — path to *_arm_level_calls_ordered.txt
#   cells_path     — path to cells CSV
#   as_path        — path to per-cell AS CSV (used when group_by="as_group")
#   subclones_path — path to sample_clusters_louvain_seurat.csv
#                    (used when group_by="subclone")
#   group_by       — "as_group" or "subclone"
#   keep_types     — list of type labels to display (None = all groups)
#
# Output: matplotlib figure

def plot_arm_level_CN_heatmap_groups(
    calls_path,
    cells_path,
    paper,
    sample,
    as_path=None,
    subclones_path=None,
    group_by='as_group',
    keep_types=None,
):
    """
    Plot an arm-level CNA heatmap with cells grouped by AS group or CN subclone.

    group_by='as_group':
        Cells labelled as lowAS / midAS / highAS (from pred column) or NonMalignant.
        Requires as_path.

    group_by='subclone':
        Cells labelled as subclone_<N> (from MergedCluster column) or Non malignant.
        Requires subclones_path.
        Cell names in the subclones file are assumed to use '.' as separator;
        they are converted to match the index of the calls matrix.
    """
    calls = pd.read_csv(calls_path, sep='\t', index_col='sample')
    calls = calls.drop(['Xp', 'Xq', 'Yp', 'Yq'], axis=1, errors='ignore')
    calls = calls.replace(CNA_VALUE_MAP)

    cells_df     = pd.read_csv(cells_path, encoding='unicode_escape')
    sample_cells = cells_df[cells_df['sample'].astype(str) == sample].copy()
    sample_cells['cell_name'] = sample_cells['cell_name'].astype(str)

    if group_by == 'as_group':
        as_df = pd.read_csv(as_path, encoding='unicode_escape')
        as_df.rename(columns={as_df.columns[0]: 'cell_name'}, inplace=True)
        as_df['cell_name'] = as_df['cell_name'].astype(str)
        sample_cells = sample_cells.merge(
            as_df[['cell_name', 'pred', 'AS']], on='cell_name', how='left'
        )
        sample_cells['type'] = sample_cells['pred'].combine_first(
            sample_cells['cell_type']
        )
        sample_cells['type'] = sample_cells['type'].where(
            sample_cells['type'].isin(['lowAS', 'midAS', 'highAS']), 'NonMalignant'
        )

    elif group_by == 'subclone':
        clones = pd.read_csv(subclones_path, encoding='unicode_escape')
        clones.rename(columns={clones.columns[1]: 'cell_name'}, inplace=True)
        clones['cell_name']     = clones['cell_name'].astype(str).str.replace('-', '.', regex=False)
        clones['MergedCluster'] = 'subclone_' + clones['MergedCluster'].astype(str)
        sample_cells = sample_cells.merge(
            clones[['cell_name', 'MergedCluster']], on='cell_name', how='left'
        )
        sample_cells['type'] = sample_cells['MergedCluster'].combine_first(
            sample_cells['cell_type']
        )
        is_subclone = sample_cells['type'].astype(str).str.fullmatch(r'subclone_\d+')
        sample_cells['type'] = sample_cells['type'].where(is_subclone, 'Non malignant')

    else:
        raise ValueError(f"group_by must be 'as_group' or 'subclone', got '{group_by}'")

    # Filter to selected groups and sort
    sample_cells_filt = sample_cells[sample_cells['cell_name'].isin(calls.index)]
    if keep_types is not None:
        sample_cells_filt = sample_cells_filt[
            sample_cells_filt['type'].isin(keep_types)
        ]
    sample_cells_filt = sample_cells_filt.sort_values('type')
    calls = calls.loc[sample_cells_filt['cell_name']]

    cmap = mcolors.LinearSegmentedColormap.from_list(
        'cna_3class', ['blue', 'white', 'red'], N=3
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    hm = sns.heatmap(
        calls, ax=ax,
        cmap=cmap, vmin=-1, vmax=1,
        yticklabels=False,
        cbar_kws={'ticks': [-1, 0, 1], 'shrink': 0.3,
                  'anchor': (0, 1), 'label': 'Arm CN'},
    )

    _add_group_dividers(ax, sample_cells_filt)

    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('Arm CN', size=14)
    ax.set_ylabel('Cells CN', fontsize=20)
    ax.set_xlabel('Chromosome Arm', fontsize=20)
    ax.set_title(f'{paper}  {sample}', fontsize=20)
    plt.xticks(rotation=90, fontweight='light', fontsize=11)

    return fig
