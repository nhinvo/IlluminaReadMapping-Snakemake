import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path 

def plot_cov_hist(plot_data):
    """
    """
    num_subplots = len(plot_data)
    num_cols = int(num_subplots ** 0.5 + 1)
    num_rows = (num_subplots + num_cols - 1) // num_cols

    # make figure
    figsize = (num_cols * 4.5, num_rows * 4.5)
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)
    axes = axes.flatten()  

    for count, ax in enumerate(axes):
        try: 
            sample_tuple = plot_data[count]
        except IndexError: 
            "More subplots specified than samples available. Leaving subplot empty."
            ax.set_axis_off()  # hide unused subplots 
            continue 

        ax.hist(sample_tuple[1], bins=50, alpha=0.9, color='#43a2ca')
        ax.set_title(sample_tuple[0])  
        ax.set_ylabel('counts')
        ax.set_xlabel('coverage')

    # title and axis labels
    fig.suptitle('Coverage Distribution Histogram', fontsize=16)

    plt.tight_layout()
    # plt.subplots_adjust(top=0.90, bottom=0.15, left=0.12, right=0.95)
    plt.savefig(snakemake.output['cov_distribution'])
    plt.close() 

def process_pos(input_paths):
    """
    Returns df with cols: [sample, median genome cov]
        - sample: sample name
        - median genome cov: median coverage across genome 
    """
    df_data = []
    plot_data = []

    for fpath in input_paths:
        sname = Path(fpath).stem.replace('_pos', '')
        df = pd.read_table(fpath, names=['contig_id', 'pos', 'count'])

        # obtain genome median coverage
        median_cov = df['count'].median()
        df_data.append({'sample': sname, 'median genome coverage': median_cov})

        # obtain array of depth (coverage)
        count_array = np.array(df['count'].values.tolist())
        plot_data.append((sname, count_array))

    # plot distribution of coverage 
    plot_cov_hist(plot_data)  # histograms    

    # df of median genome coverage
    df = pd.DataFrame(df_data)

    return df


def process_contig(input_paths):
    """
    Returns df with cols: [contig_name, read_count, percent, sample]
    """
    dfs = []
    for fpath in input_paths:
        sname = Path(fpath).stem.replace('_contig', '')
        df = pd.read_table(fpath, names=['contig_name', 'genome_length', 'read_count', 'unmapped'])

        # move the unmapped value in contig_name row '*' to mapped_reads
        df.loc[df['contig_name'] == '*', 'read_count'] = df.loc[df['contig_name'] == '*', 'unmapped']
        df['contig_name'] = df['contig_name'].str.replace('*', 'unmapped_reads')

        df = df[['contig_name', 'read_count']]

        df['read_count'] = df['read_count'].astype(int)

        total = df['read_count'].sum()  # total read count 
        df['percent'] = (df['read_count'] / total) * 100

        df = df.sort_values(by=['percent'], ascending=False)
        df = df[df['percent'] > 1]

        dfs.append(df)

    df = pd.concat(dfs)
    return df

def process_stats(input_paths):
    """
    Returns df with cols: [sample, total read, read mapped, percent mapped]
    """
    df_data = []
    for fpath in input_paths:
        sname = Path(fpath).stem.replace('_mapping_stats', '')
        df = pd.read_table(fpath, names=['desc', 'stats', 'notes'])
        df['desc'] = df['desc'].str.strip().str[:-1]  # remove the ':' at the end

        total_reads = df.loc[df['desc'] == 'raw total sequences', 'stats'].values[0]
        reads_mapped = df.loc[df['desc'] == 'reads mapped', 'stats'].values[0]

        df_data.append({
            'sample': sname, 
            'total read': total_reads, 
            'read mapped': reads_mapped, 
        })

    df = pd.DataFrame(df_data)
    df['percent_mapped'] = (df['read mapped'] / df['total read']) * 100 

    return df

def plot_summary_bar(df):
    """
    Plot stacked bars of mapped/unmapped reads for each sample. 
    """
    # prep df
    df['percent_unmapped'] = 100 - df['percent_mapped']
    df = df[['sample', 'percent_mapped', 'percent_unmapped']]
    df = df.set_index(['sample'])

    # Plot 
    fig_width = len(df) * 0.5
    fig, ax = plt.subplots(figsize=(fig_width, 5))  
    colors = ['#8cbf88', '#dedede'] 
    df.plot(kind='bar', stacked=True, width=0.9, color=colors, ax=ax)

    # Legend
    ax.legend(labels=['Mapped', 'Un-Mapped'], loc='upper left', bbox_to_anchor=(1, 1))

    # Ticks 
    plt.xticks(rotation=90)

    # Titles and Labels 
    plt.title('Mapped vs Unmapped Reads')
    plt.ylabel('Read Percent')
    plt.xlabel('Sample Name')

    plt.tight_layout()

    plt.savefig(snakemake.output['summary_mapping_stats'])
    plt.close()    

def main():
    median_cov_df = process_pos(snakemake.input['pos'])  # process & plot position coverage 
    contig_df = process_contig(snakemake.input['contig'])  # process contig coverage 
    stats_df = process_stats(snakemake.input['stats'])  # obtain summary stats (i.e. % reads mapped)

    # merge and save data 
    df = pd.merge(median_cov_df, stats_df, on='sample', how='outer')

    # save data as Excel file 
    with pd.ExcelWriter(snakemake.output['final_excel']) as writer:
        df.to_excel(writer, sheet_name='mapping_stats', index=False)
        contig_df.to_excel(writer, sheet_name='contig_mapping_stats', index=False)

    # plot data 
    plot_summary_bar(df)

main()