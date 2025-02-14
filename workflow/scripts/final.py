import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path 

def process_pos(input_paths):
    """
    Returns df with cols: [sample, median genome cov]
        - sample: sample name
        - median genome cov: median coverage across genome 
    """
    df_data = []
    for fpath in input_paths:
        sname = Path(fpath).stem.replace('_pos', '')
        df = pd.read_table(fpath, names=['contig_id', 'pos', 'count'])
        median_cov = df['count'].median()

        df_data.append({'sample': sname, 'median genome coverage': median_cov})

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

    plot_outdir = snakemake.input['plot_outdir']
    plt.savefig(f'{plot_outdir}/SummaryMapping.png')
    plt.close()    

def plot_summary(df):
    """
    """
    plot_summary_bar(df)

def main():
    pos_df = process_pos(snakemake.input['pos'])  # process position coverage 
    contig_df = process_contig(snakemake.input['contig'])  # process contig coverage 
    stats_df = process_stats(snakemake.input['stats'])  # obtain summary stats (i.e. % reads mapped)

    # merge and save data 
    df = pd.merge(pos_df, stats_df, on='sample', how='outer')

    # save data as Excel file 
    with pd.ExcelWriter(snakemake.output[0]) as writer:
        df.to_excel(writer, sheet_name='mapping_stats', index=False)
        contig_df.to_excel(writer, sheet_name='contig_mapping_stats', index=False)

    # plot data 
    plot_summary(df)

main()