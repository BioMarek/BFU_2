#!/usr/bin/env python3
"""Cluster analysis of Rca3_R3 centromeric repeats from out_hor.rawTOP5b_rc.fa."""

from __future__ import annotations

import collections
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import StandardScaler

FASTA_PATH = Path(__file__).resolve().with_name("out_hor.rawTOP5b_rc.fa")
OUTPUT_DIR = Path(__file__).resolve().parent

REPEAT_ORDER = ["R1", "R2", "R3", "R4", "R5"]
ACTIVE_REPEATS = ("R1", "R2", "R3")
INACTIVE_REPEATS = ("R4", "R5")
ACTIVITY_ORDER = ["active", "inactive"]
ACTIVITY_COLORS = {"active": "#2ca02c", "inactive": "#d62728"}

MIN_SEQUENCE_LENGTH = 50
PERIOD_RANGE = (140, 180)
MAX_CONSENSUS_ITERATIONS = 5
N_CLUSTERS = 2
RANDOM_STATE = 42

BASE_CODE = np.full(256, -1, dtype=np.int16)
for _code, _base in enumerate(b"ACGT"):
    BASE_CODE[_base] = _code


def encode(sequence: str) -> np.ndarray:
    return BASE_CODE[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]


def decode(encoded: np.ndarray) -> str:
    return "".join("ACGT"[int(base)] for base in encoded)


def read_fasta(path: Path) -> tuple[pd.DataFrame, int, int]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.upper())
        if header is not None:
            records.append((header, "".join(chunks)))

    rows = []
    skipped = 0
    for header, sequence in records:
        if len(sequence) < MIN_SEQUENCE_LENGTH or set(sequence) - set("ACGT"):
            skipped += 1
            continue
        repeat_match = re.match(r"(R[1-5])", header)
        if not repeat_match:
            skipped += 1
            continue
        coordinate_match = re.search(r"::(\d+)-(\d+)::([+-])", header)
        locus_match = re.search(r"L(\d+)", header)
        rows.append(
            {
                "header": header,
                "repeat": repeat_match.group(1),
                "locus": f"L{locus_match.group(1)}" if locus_match else "NA",
                "start": int(coordinate_match.group(1)) if coordinate_match else pd.NA,
                "end": int(coordinate_match.group(2)) if coordinate_match else pd.NA,
                "strand": coordinate_match.group(3) if coordinate_match else "NA",
                "length": len(sequence),
                "sequence": sequence,
            }
        )

    return pd.DataFrame(rows), skipped, len(records)


def detect_period(sequences: list[str], low: int, high: int, max_sequences: int = 800) -> tuple[int, dict[int, float]]:
    candidates = [sequence for sequence in sequences if len(sequence) >= 2 * high]
    if not candidates:
        candidates = [sequence for sequence in sequences if len(sequence) >= 2 * (high - 20)]
    rng = np.random.default_rng(RANDOM_STATE)
    if len(candidates) > max_sequences:
        indices = rng.choice(len(candidates), size=max_sequences, replace=False)
        candidates = [candidates[index] for index in indices]

    scores = {}
    for period in range(low, high + 1):
        fractions = []
        for sequence in candidates:
            left = sequence[:-period]
            right = sequence[period:]
            fractions.append(sum(a != b for a, b in zip(left, right)) / len(left))
        scores[period] = float(np.mean(fractions))

    return min(scores, key=scores.get), scores


def phase_mismatches(encoded: np.ndarray, consensus: np.ndarray) -> np.ndarray:
    n = encoded.size
    periods = consensus.size
    phases = np.arange(periods)
    indices = (np.arange(n)[None, :] + phases[:, None]) % periods
    return (consensus[indices] != encoded[None, :]).sum(axis=1)


def build_consensus(encoded_sequences: list[np.ndarray], sequences: list[str], period: int) -> tuple[np.ndarray, str, int]:
    seed_counts = collections.Counter(sequence[:period] for sequence in sequences if len(sequence) >= period)
    seed = seed_counts.most_common(1)[0][0]
    consensus = encode(seed)
    iterations = 0

    for iterations in range(1, MAX_CONSENSUS_ITERATIONS + 1):
        counts = np.zeros((period, 4), dtype=np.int64)
        for encoded in encoded_sequences:
            phase = int(np.argmin(phase_mismatches(encoded, consensus)))
            positions = (np.arange(encoded.size) + phase) % period
            for base in range(4):
                counts[:, base] += np.bincount(positions[encoded == base], minlength=period)
        updated = counts.argmax(axis=1).astype(np.int16)
        if np.array_equal(updated, consensus):
            break
        consensus = updated

    return consensus, decode(consensus), iterations


def annotate_mutations(df: pd.DataFrame, encoded_sequences: list[np.ndarray], consensus: np.ndarray) -> pd.DataFrame:
    mutations = []
    rates = []
    per_monomer = []
    phases = []

    for encoded, length in zip(encoded_sequences, df["length"].to_numpy()):
        mismatch = phase_mismatches(encoded, consensus)
        count = int(mismatch.min())
        mutations.append(count)
        rates.append(count / length)
        per_monomer.append(count / (length / consensus.size))
        phases.append(int(mismatch.argmin()))

    annotated = df.copy()
    annotated["mutations"] = mutations
    annotated["mutation_rate"] = rates
    annotated["mutations_per_monomer"] = per_monomer
    annotated["phase"] = phases
    annotated["activity"] = np.where(annotated["repeat"].isin(ACTIVE_REPEATS), "active", "inactive")
    return annotated


def fit_kmeans(df: pd.DataFrame, features: list[str], n_clusters: int = N_CLUSTERS) -> tuple[np.ndarray, np.ndarray, KMeans]:
    values = df[features].to_numpy(dtype=float)
    scaled = StandardScaler().fit_transform(values)
    model = KMeans(n_clusters=n_clusters, n_init=25, random_state=RANDOM_STATE).fit(scaled)
    labels = model.labels_
    order = (
        df.assign(cluster=labels)
        .groupby("cluster")["length"]
        .mean()
        .sort_values()
        .index.tolist()
    )
    remap = {old: new for new, old in enumerate(order)}
    ordered = np.array([remap[label] for label in labels])
    return ordered, scaled, model


def balanced_accuracy(y_true: np.ndarray, cluster_labels: np.ndarray) -> float:
    best = 0.0
    for active_cluster in np.unique(cluster_labels):
        predicted_active = cluster_labels == active_cluster
        sensitivity = predicted_active[y_true == 1].mean()
        specificity = (~predicted_active)[y_true == 0].mean()
        best = max(best, (sensitivity + specificity) / 2)
    return float(best)


def cluster_summary(df: pd.DataFrame, labels: np.ndarray) -> dict[str, float]:
    activity = np.where(df["activity"].to_numpy() == "active", 1, 0)
    table = pd.crosstab(df["activity"], labels).reindex(ACTIVITY_ORDER)
    ari = adjusted_rand_score(activity, labels)
    accuracy = balanced_accuracy(activity, labels)
    if len(np.unique(labels)) > 1:
        silhouette = silhouette_score(
            StandardScaler().fit_transform(df[["length", "mutations", "mutation_rate"]].to_numpy(dtype=float)),
            labels,
            sample_size=min(3000, len(df)),
            random_state=RANDOM_STATE,
        )
    else:
        silhouette = float("nan")
    try:
        chi2, p_value, _, _ = chi2_contingency(table.to_numpy())
    except ValueError:
        chi2, p_value = float("nan"), float("nan")
    return {
        "ari": float(ari),
        "balanced_accuracy": accuracy,
        "silhouette": float(silhouette),
        "chi2": float(chi2),
        "p_value": float(p_value),
        "table": table,
    }


def make_main_figure(df: pd.DataFrame, verdict: str) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(2, 3, figsize=(19, 11))

    cluster_colors = {0: "#1f77b4", 1: "#ff7f0e"}
    point_colors = df["cluster"].map(cluster_colors)
    activity_colors = df["activity"].map(ACTIVITY_COLORS)
    mutation_cap = df["mutations"].quantile(0.995) * 1.05

    ax = axes[0, 0]
    ax.scatter(df["length"], df["mutations"], c=point_colors, s=14, alpha=0.55, linewidths=0)
    ax.set_ylim(0, mutation_cap)
    ax.set_xlabel("Sequence length (bp)")
    ax.set_ylabel("Mutations vs consensus")
    ax.set_title("K-means clusters: length + mutations")
    ax.legend(
        handles=[Patch(color=color, label=f"Cluster {label}") for label, color in cluster_colors.items()],
        loc="upper left",
    )

    ax = axes[0, 1]
    ax.scatter(df["length"], df["mutations"], c=activity_colors, s=14, alpha=0.55, linewidths=0)
    ax.set_ylim(0, mutation_cap)
    ax.set_xlabel("Sequence length (bp)")
    ax.set_ylabel("Mutations vs consensus")
    ax.set_title("Same data coloured by activity")
    ax.legend(
        handles=[Patch(color=ACTIVITY_COLORS[label], label=label) for label in ACTIVITY_ORDER],
        loc="upper left",
    )

    ax = axes[0, 2]
    sns.boxplot(
        data=df,
        x="repeat",
        y="length",
        hue="activity",
        order=REPEAT_ORDER,
        hue_order=ACTIVITY_ORDER,
        palette=ACTIVITY_COLORS,
        dodge=False,
        ax=ax,
        fliersize=2,
    )
    ax.set_xlabel("Repeat class")
    ax.set_ylabel("Sequence length (bp)")
    ax.set_title("Length by repeat class")
    ax.legend_.remove()

    ax = axes[1, 0]
    rate_cap = df["mutation_rate"].quantile(0.995) * 1.05
    sns.boxplot(
        data=df,
        x="repeat",
        y="mutation_rate",
        hue="activity",
        order=REPEAT_ORDER,
        hue_order=ACTIVITY_ORDER,
        palette=ACTIVITY_COLORS,
        dodge=False,
        ax=ax,
        fliersize=2,
    )
    ax.set_ylim(0, rate_cap)
    ax.set_xlabel("Repeat class")
    ax.set_ylabel("Mutations per bp")
    ax.set_title("Mutation rate by repeat class")
    ax.legend(
        handles=[Patch(color=ACTIVITY_COLORS[label], label=label) for label in ACTIVITY_ORDER],
        loc="upper left",
    )

    ax = axes[1, 1]
    composition = pd.crosstab(df["repeat"], df["cluster"]).reindex(REPEAT_ORDER)
    sns.heatmap(composition, annot=True, fmt="d", cmap="Blues", ax=ax, cbar_kws={"label": "sequences"})
    ax.set_xlabel("K-means cluster")
    ax.set_ylabel("Repeat class")
    ax.set_title("Cluster composition")

    ax = axes[1, 2]
    activity_composition = pd.crosstab(df["cluster"], df["activity"], normalize="index").reindex(
        columns=ACTIVITY_ORDER
    )
    activity_composition.plot(
        kind="bar",
        stacked=True,
        color=[ACTIVITY_COLORS[label] for label in ACTIVITY_ORDER],
        ax=ax,
        rot=0,
        legend=False,
    )
    for position, (_, row) in enumerate(activity_composition.iterrows()):
        ax.text(
            position,
            row["active"] / 2,
            f"{row['active']:.1%}",
            ha="center",
            va="center",
            color="white",
            fontweight="bold",
        )
        ax.text(
            position,
            row["active"] + row["inactive"] / 2,
            f"{row['inactive']:.1%}",
            ha="center",
            va="center",
            color="white",
            fontweight="bold",
        )
    ax.set_xlabel("K-means cluster")
    ax.set_ylabel("Fraction of sequences")
    ax.set_ylim(0, 1)
    ax.set_title("Activity composition per cluster")

    fig.suptitle(f"Rca3_R3 centromeric repeat clustering — hypothesis {verdict}", fontsize=22)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUTPUT_DIR / "cluster_analysis_R.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def make_distribution_figure(df: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 3, figsize=(19, 6))
    rate_cap = df["mutation_rate"].quantile(0.995) * 1.05
    monomer_cap = df["mutations_per_monomer"].quantile(0.995) * 1.05

    ax = axes[0]
    sns.violinplot(
        data=df,
        x="repeat",
        y="mutation_rate",
        hue="activity",
        order=REPEAT_ORDER,
        hue_order=ACTIVITY_ORDER,
        palette=ACTIVITY_COLORS,
        dodge=False,
        cut=0,
        ax=ax,
    )
    ax.set_ylim(0, rate_cap)
    ax.set_xlabel("Repeat class")
    ax.set_ylabel("Mutations per bp")
    ax.set_title("Mutation-rate distributions")
    ax.legend_.remove()

    ax = axes[1]
    sns.violinplot(
        data=df,
        x="repeat",
        y="mutations_per_monomer",
        hue="activity",
        order=REPEAT_ORDER,
        hue_order=ACTIVITY_ORDER,
        palette=ACTIVITY_COLORS,
        dodge=False,
        cut=0,
        ax=ax,
    )
    ax.set_ylim(0, monomer_cap)
    ax.set_xlabel("Repeat class")
    ax.set_ylabel("Mutations per repeat unit")
    ax.set_title("Mutations per repeat unit")
    ax.legend(
        handles=[Patch(color=ACTIVITY_COLORS[label], label=label) for label in ACTIVITY_ORDER],
        loc="upper left",
    )

    ax = axes[2]
    for activity in ACTIVITY_ORDER:
        values = np.sort(df.loc[df["activity"] == activity, "mutation_rate"].to_numpy())
        if values.size:
            cumulative = np.arange(1, values.size + 1) / values.size
            ax.plot(values, cumulative, color=ACTIVITY_COLORS[activity], label=activity, linewidth=2)
    ax.set_xlim(0, rate_cap)
    ax.set_xlabel("Mutations per bp")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("Empirical cumulative distribution")
    ax.legend()

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "mutation_distributions.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not FASTA_PATH.exists():
        raise FileNotFoundError(f"Input FASTA not found: {FASTA_PATH}")

    df, skipped, total_records = read_fasta(FASTA_PATH)
    print(f"Loaded {len(df)} sequences ({skipped} skipped of {total_records} records)")
    print("Sequences per repeat class:")
    print(df["repeat"].value_counts().reindex(REPEAT_ORDER).to_string())

    period, period_scores = detect_period(df["sequence"].tolist(), *PERIOD_RANGE)
    print(f"\nDetected repeat period: {period} bp (mismatch score {period_scores[period]:.3f})")
    neighbouring = {p: round(period_scores[p], 3) for p in range(period - 2, period + 3)}
    print(f"Neighbouring period scores: {neighbouring}")

    encoded_sequences = [encode(sequence) for sequence in df["sequence"]]
    consensus, consensus_string, iterations = build_consensus(
        encoded_sequences, df["sequence"].tolist(), period
    )
    print(f"Consensus built in {iterations} iteration(s), length {len(consensus_string)} bp")
    print(f"Consensus: {consensus_string}")

    df = annotate_mutations(df, encoded_sequences, consensus)

    summary = (
        df.groupby("repeat")[["length", "mutations", "mutation_rate", "mutations_per_monomer"]]
        .agg(["median", "mean"])
        .reindex(REPEAT_ORDER)
    )
    print("\nSequence statistics by repeat class:")
    print(summary.round(4).to_string())

    print("\nLength and mutation summary:")
    overview = df.groupby("repeat").agg(
        n=("length", "size"),
        length_median=("length", "median"),
        length_iqr=("length", lambda values: values.quantile(0.75) - values.quantile(0.25)),
        mutations_median=("mutations", "median"),
        mutation_rate_median=("mutation_rate", "median"),
        mutation_rate_mean=("mutation_rate", "mean"),
    )
    print(overview.reindex(REPEAT_ORDER).round(4).to_string())

    labels, scaled, model = fit_kmeans(df, ["length", "mutations"])
    df["cluster"] = labels

    print("\nK-means (k=2) on length + mutations, cluster sizes:")
    print(pd.Series(labels).value_counts().sort_index().to_string())
    print("\nRepeat class by cluster:")
    print(pd.crosstab(df["repeat"], df["cluster"]).reindex(REPEAT_ORDER).to_string())
    print("\nActivity by cluster:")
    activity_table = pd.crosstab(df["activity"], df["cluster"]).reindex(ACTIVITY_ORDER)
    print(activity_table.to_string())

    metrics = cluster_summary(df, labels)
    print(f"\nAdjusted Rand index (cluster vs activity): {metrics['ari']:.4f}")
    print(f"Balanced accuracy of best cluster-to-activity mapping: {metrics['balanced_accuracy']:.4f}")
    print(f"Silhouette score: {metrics['silhouette']:.4f}")
    print(f"Chi-square test of cluster vs activity: chi2={metrics['chi2']:.2f}, p={metrics['p_value']:.3e}")

    verdict = (
        "CONFIRMED"
        if metrics["balanced_accuracy"] >= 0.85 and metrics["ari"] >= 0.5
        else "REJECTED"
    )
    print(f"\n>>> HYPOTHESIS {verdict}")
    if verdict == "REJECTED":
        print(
            "    The two clusters do not separate active (R1-R3) from inactive (R4-R5) repeats; "
            "both clusters contain a mixture of activity states."
        )
        print(
            "    The dominant structure is sequence length: R2/R5 (~159 bp) fall into one cluster "
            "while R1/R3/R4 (~318-477 bp) fall into the other."
        )
    else:
        print("    Active (R1-R3) and inactive (R4-R5) repeats occupy different clusters.")

    print("\nSensitivity of the two-cluster solution to feature choice:")
    feature_sets = {
        "length + mutations": ["length", "mutations"],
        "length + mutation_rate": ["length", "mutation_rate"],
        "mutations only": ["mutations"],
        "mutation_rate only": ["mutation_rate"],
        "length only": ["length"],
        "length + mutations + mutation_rate": ["length", "mutations", "mutation_rate"],
    }
    sensitivity_rows = []
    for name, features in feature_sets.items():
        feature_labels, feature_scaled, _ = fit_kmeans(df, features)
        activity = np.where(df["activity"].to_numpy() == "active", 1, 0)
        ari = adjusted_rand_score(activity, feature_labels)
        accuracy = balanced_accuracy(activity, feature_labels)
        clusters = pd.crosstab(df["repeat"], feature_labels).reindex(REPEAT_ORDER)
        active_fraction = (
            df.assign(cluster=feature_labels)
            .groupby("cluster")["activity"]
            .apply(lambda values: float((values == "active").mean()))
        )
        dominant_repeats = {}
        for cluster in [0, 1]:
            column = clusters.get(cluster, pd.Series(0, index=clusters.index))
            dominant_repeat = column.idxmax()
            dominant_share = column.max() / column.sum()
            dominant_repeats[cluster] = f"{dominant_repeat} ({dominant_share:.0%})"
        sensitivity_rows.append(
            {
                "features": name,
                "ari": ari,
                "balanced_accuracy": accuracy,
                "active_fraction_cluster_0": active_fraction.get(0, float("nan")),
                "active_fraction_cluster_1": active_fraction.get(1, float("nan")),
                "dominant_cluster_0": dominant_repeats[0],
                "dominant_cluster_1": dominant_repeats[1],
            }
        )
    sensitivity = pd.DataFrame(sensitivity_rows)
    print(sensitivity.round(4).to_string(index=False))
    sensitivity.to_csv(OUTPUT_DIR / "clustering_sensitivity.csv", index=False)

    labels3, _, _ = fit_kmeans(df, ["length", "mutations"], n_clusters=3)
    print("\nLength + mutation structure with k=3 (sanity check for length classes):")
    print(pd.crosstab(df["repeat"], labels3).reindex(REPEAT_ORDER).to_string())

    print("\nLength-matched comparisons of mutation rate (Mann-Whitney U):")
    for left, right in [("R2", "R5"), ("R3", "R4")]:
        left_values = df.loc[df["repeat"] == left, "mutation_rate"]
        right_values = df.loc[df["repeat"] == right, "mutation_rate"]
        statistic, p_value = mannwhitneyu(left_values, right_values, alternative="two-sided")
        print(
            f"  {left} vs {right}: median {left_values.median():.4f} vs {right_values.median():.4f}, "
            f"mean {left_values.mean():.4f} vs {right_values.mean():.4f}, p={p_value:.3e}"
        )

    statistic, p_value = kruskal(*[df.loc[df["repeat"] == repeat, "mutation_rate"] for repeat in REPEAT_ORDER])
    print(f"\nKruskal-Wallis across repeat classes: H={statistic:.2f}, p={p_value:.3e}")

    make_main_figure(df, verdict)
    make_distribution_figure(df)

    output_columns = [
        "header",
        "repeat",
        "locus",
        "start",
        "end",
        "strand",
        "length",
        "mutations",
        "mutation_rate",
        "mutations_per_monomer",
        "cluster",
        "activity",
        "sequence",
    ]
    df[output_columns].to_csv(OUTPUT_DIR / "sequence_analysis_results.csv", index=False)

    print("\nSaved outputs:")
    print(f"  {OUTPUT_DIR / 'cluster_analysis_R.png'}")
    print(f"  {OUTPUT_DIR / 'mutation_distributions.png'}")
    print(f"  {OUTPUT_DIR / 'sequence_analysis_results.csv'}")
    print(f"  {OUTPUT_DIR / 'clustering_sensitivity.csv'}")


if __name__ == "__main__":
    main()
