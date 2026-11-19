#!/usr/bin/env python3
"""
Centromeric repeat sequence cluster analysis.
Analyzes sequences from out_hor.rawTOP5b_rc.fa to cluster R1-R5 types
based on sequence length and mutations from consensus.
"""

import re
import numpy as np
import pandas as pd
from collections import Counter
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import warnings
warnings.filterwarnings('ignore')


def parse_fasta(filepath):
    """Parse FASTA file and extract sequence type, coordinates, and sequence."""
    sequences = []
    current_header = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    seq = ''.join(current_seq)
                    sequences.append((current_header, seq))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header is not None:
            seq = ''.join(current_seq)
            sequences.append((current_header, seq))
    
    parsed = []
    for header, seq in sequences:
        match = re.match(r'(R[1-5]L\d+)', header)
        seq_type = match.group(1) if match else 'Unknown'
        if len(seq) > 0:
            parsed.append({
                'header': header,
                'seq_type': seq_type,
                'sequence': seq,
                'length': len(seq)
            })
    
    return pd.DataFrame(parsed)


def compute_consensus(sequences):
    """Compute consensus sequence from a list of sequences."""
    if not sequences:
        return ""
    
    max_len = max(len(s) for s in sequences)
    consensus = []
    
    for i in range(max_len):
        bases = [s[i] for s in sequences if i < len(s)]
        if bases:
            consensus.append(Counter(bases).most_common(1)[0][0])
        else:
            consensus.append('N')
    
    return ''.join(consensus)


def count_mutations(seq, consensus):
    """Count mutations (differences) between sequence and consensus."""
    min_len = min(len(seq), len(consensus))
    mutations = sum(1 for i in range(min_len) if seq[i] != consensus[i])
    mutations += abs(len(seq) - len(consensus))
    return mutations


def main():
    fasta_file = 'out_hor.rawTOP5b_rc.fa'
    
    print("Parsing FASTA file...")
    df = parse_fasta(fasta_file)
    print(f"Loaded {len(df)} sequences")
    print(f"Sequence types: {df['seq_type'].value_counts().to_dict()}")
    
    df['R_group'] = df['seq_type'].str.extract(r'(R[1-5])')[0]
    print(f"R groups: {df['R_group'].value_counts().to_dict()}")
    
    print("\nComputing consensus sequences...")
    consensus_by_group = {}
    for group in sorted(df['R_group'].unique()):
        group_seqs = df[df['R_group'] == group]['sequence'].tolist()
        consensus_by_group[group] = compute_consensus(group_seqs)
        print(f"  {group}: consensus length = {len(consensus_by_group[group])}")
    
    overall_consensus = compute_consensus(df['sequence'].tolist())
    print(f"Overall consensus length = {len(overall_consensus)}")
    
    print("\nCalculating mutations from consensus...")
    df['mutations_overall'] = df['sequence'].apply(lambda s: count_mutations(s, overall_consensus))
    df['mutations_group'] = df.apply(
        lambda row: count_mutations(row['sequence'], consensus_by_group[row['R_group']]), axis=1
    )
    
    mutation_rate_overall = df['mutations_overall'] / df['length']
    mutation_rate_group = df['mutations_group'] / df['length']
    df['mutation_rate_overall'] = mutation_rate_overall
    df['mutation_rate_group'] = mutation_rate_group
    
    print("\nMutation statistics (overall consensus):")
    print(df.groupby('R_group')['mutations_overall'].describe())
    print("\nMutation rate statistics (overall consensus):")
    print(df.groupby('R_group')['mutation_rate_overall'].describe())
    
    print("\nMutation statistics (group consensus):")
    print(df.groupby('R_group')['mutations_group'].describe())
    print("\nMutation rate statistics (group consensus):")
    print(df.groupby('R_group')['mutation_rate_group'].describe())
    
    print("\nPerforming cluster analysis...")
    features = df[['length', 'mutations_overall', 'mutation_rate_overall']].copy()
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)
    
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=25)
    df['cluster'] = kmeans.fit_predict(features_scaled)
    
    print(f"\nCluster sizes: {df['cluster'].value_counts().to_dict()}")
    print("\nCluster composition by R_group:")
    cluster_comp = pd.crosstab(df['R_group'], df['cluster'])
    print(cluster_comp)
    
    print("\nCluster composition by seq_type:")
    cluster_comp2 = pd.crosstab(df['seq_type'], df['cluster'])
    print(cluster_comp2)
    
    active_in_cluster = {}
    for c in [0, 1]:
        active = df[(df['cluster'] == c) & (df['R_group'].isin(['R1', 'R2', 'R3']))].shape[0]
        inactive = df[(df['cluster'] == c) & (df['R_group'].isin(['R4', 'R5']))].shape[0]
        active_in_cluster[c] = {'active': active, 'inactive': inactive}
    
    print("\nActive (R1,R2,R3) vs Inactive (R4,R5) in clusters:")
    for c, counts in active_in_cluster.items():
        print(f"  Cluster {c}: Active={counts['active']}, Inactive={counts['inactive']}")
    
    active_majority = max(active_in_cluster, key=lambda c: active_in_cluster[c]['active'] / 
                          (active_in_cluster[c]['active'] + active_in_cluster[c]['inactive'] + 1e-9))
    inactive_majority = max(active_in_cluster, key=lambda c: active_in_cluster[c]['inactive'] / 
                            (active_in_cluster[c]['active'] + active_in_cluster[c]['inactive'] + 1e-9))
    
    if active_majority != inactive_majority:
        print(f"\n>>> HYPOTHESIS CONFIRMED: Active (R1,R2,R3) and Inactive (R4,R5) separate into different clusters!")
        print(f"    Cluster {active_majority}: predominantly active repeats")
        print(f"    Cluster {inactive_majority}: predominantly inactive repeats")
    else:
        print(f"\n>>> HYPOTHESIS REJECTED: Active and inactive repeats do not clearly separate into different clusters.")
    
    print("\nGenerating plots...")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    for ax, (group, color) in zip(axes[0, 0].collections if hasattr(axes[0, 0], 'collections') else [], []):
        pass
    
    scatter = axes[0, 0].scatter(df['length'], df['mutations_overall'], 
                                  c=df['cluster'], cmap='Set1', alpha=0.6, s=30)
    axes[0, 0].set_xlabel('Sequence Length')
    axes[0, 0].set_ylabel('Mutations from Overall Consensus')
    axes[0, 0].set_title('Clusters: Length vs Mutations (Overall Consensus)')
    plt.colorbar(scatter, ax=axes[0, 0], label='Cluster')
    
    for group in sorted(df['R_group'].unique()):
        mask = df['R_group'] == group
        axes[0, 1].scatter(df.loc[mask, 'length'], df.loc[mask, 'mutations_overall'], 
                           label=group, alpha=0.6, s=30)
    axes[0, 1].set_xlabel('Sequence Length')
    axes[0, 1].set_ylabel('Mutations from Overall Consensus')
    axes[0, 1].set_title('R Groups: Length vs Mutations (Overall Consensus)')
    axes[0, 1].legend()
    
    scatter2 = axes[0, 2].scatter(df['length'], df['mutation_rate_overall'], 
                                   c=df['cluster'], cmap='Set1', alpha=0.6, s=30)
    axes[0, 2].set_xlabel('Sequence Length')
    axes[0, 2].set_ylabel('Mutation Rate (Overall Consensus)')
    axes[0, 2].set_title('Clusters: Length vs Mutation Rate')
    plt.colorbar(scatter2, ax=axes[0, 2], label='Cluster')
    
    axes[1, 0].boxplot([df[df['R_group'] == g]['length'].values for g in sorted(df['R_group'].unique())],
                        tick_labels=sorted(df['R_group'].unique()))
    axes[1, 0].set_ylabel('Sequence Length')
    axes[1, 0].set_title('Length Distribution by R Group')
    axes[1, 0].tick_params(axis='x', rotation=45)
    
    axes[1, 1].boxplot([df[df['R_group'] == g]['mutations_overall'].values for g in sorted(df['R_group'].unique())],
                        tick_labels=sorted(df['R_group'].unique()))
    axes[1, 1].set_ylabel('Mutations from Overall Consensus')
    axes[1, 1].set_title('Mutations by R Group')
    axes[1, 1].tick_params(axis='x', rotation=45)
    
    axes[1, 2].boxplot([df[df['R_group'] == g]['mutation_rate_overall'].values for g in sorted(df['R_group'].unique())],
                        tick_labels=sorted(df['R_group'].unique()))
    axes[1, 2].set_ylabel('Mutation Rate')
    axes[1, 2].set_title('Mutation Rate by R Group')
    axes[1, 2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('cluster_analysis_R.png', dpi=150, bbox_inches='tight')
    print("Saved cluster_analysis_R.png")
    
    fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
    
    sns.violinplot(data=df, x='R_group', y='mutations_overall', ax=axes2[0], order=sorted(df['R_group'].unique()))
    axes2[0].set_title('Mutations from Overall Consensus by R Group')
    axes2[0].set_xlabel('R Group')
    axes2[0].set_ylabel('Number of Mutations')
    
    sns.violinplot(data=df, x='R_group', y='mutation_rate_overall', ax=axes2[1], order=sorted(df['R_group'].unique()))
    axes2[1].set_title('Mutation Rate by R Group')
    axes2[1].set_xlabel('R Group')
    axes2[1].set_ylabel('Mutation Rate')
    
    plt.tight_layout()
    plt.savefig('mutation_distributions.png', dpi=150, bbox_inches='tight')
    print("Saved mutation_distributions.png")
    
    df.to_csv('sequence_analysis_results.csv', index=False)
    print("Saved sequence_analysis_results.csv")
    
    print("\nDone!")


if __name__ == '__main__':
    main()