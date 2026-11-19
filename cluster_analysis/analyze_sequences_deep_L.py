#!/usr/bin/env python3
"""Cluster analysis of Rca3_R3 centromeric repeats by locus (L) from out_hor.rawTOP5b_rc.fa."""

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

LOCUS_ORDER = ["L1", "L2", "L3"]
REPEAT_ORDER = ["R1", "R2", "R3", "R4", "R5"]
ACTIVE_REPEATS = ("R1", "R2", "R3")
INACTIVE_REPEATS = ("R4", "R5")
ACTIVITY_ORDER = ["active", "inactive"]
ACTIVITY_COLORS = {"active": "#2ca02c", "inactive": "#d62728"}
LOCUS_COLORS = {"L1": "#1f77b4", "L2": "#ff7f0e", "L3": "#9467bd"}

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
        name_match = re.match(r"(R[1-5])(L\d+)", header)
        if not name_match:
            skipped += 1
            continue
        coordinate_match = re.search(r"::(\d+)-(\d+)::([+-])", header)
        rows.append(
            {
                "header": header,
                "repeat": name_match.group(1),
                "locus": name_match.group(2),
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


def mutation_metrics(encoded: np.ndarray, length: int, consensus: np.ndarray) -> tuple[int, float, float, int]:
    mismatch = phase_mismatches(encoded, consensus)
    count = int(mismatch.min())
    return count, count / length, count / (length / consensus.size), int(mismatch.argmin())


def annotate_mutations(df: pd.DataFrame, encoded_sequences: list[np.ndarray], consensus: np.ndarray) -> pd.DataFrame:
    mutations = []
    rates = []
    per_monomer = []
    phases = []

    for encoded, length in zip(encoded_sequences, df["length"].to_numpy()):
        count, rate, per_unit, phase = mutation_metrics(encoded, length, consensus)
        mutations.append(count)
        rates.append(rate)
        per_monomer.append(per_unit)
        phases.append(phase)

    annotated = df.copy()
    annotated["mutations"] = mutations
    annotated["mutation_rate"] = rates
    annotated["mutations_per_monomer"] = per_monomer
    annotated["phase"] = phases
    annotated["activity"] = np.where(annotated["repeat"].isin(ACTIVE_REPEATS), "active", "inactive")
    return annotated


def annotate_locus_mutations(
    df: pd.DataFrame, encoded_sequences: list[np.ndarray], locus_consensuses: dict[str, np.ndarray]
) -> pd.DataFrame:
    counts = []
    rates = []
    per_monomer = []
    for encoded, length, locus in zip(encoded_sequences, df["length"].to_numpy(), df["locus"]):
        count, rate, per_unit, _ = mutation_metrics(encoded, length, locus_consensuses[locus])
        counts.append(count)
        rates.append(rate)
        per_monomer.append(per_unit)

    annotated = df.copy()
    annotated["mutations_locus"] = counts
    annotated["mutation_rate_locus"] = rates
    annotated["mutations_per_monomer_locus"] = per_monomer
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
    for positive_cluster in np.unique(cluster_labels):
        predicted_positive = cluster_labels == positive_cluster
        sensitivity = predicted_positive[y_true == 1].mean()
        specificity = (~predicted_positive)[y_true == 0].mean()
        best = max(best, (sensitivity + specificity) / 2)
    return float(best)


def cluster_summary(df: pd.DataFrame, labels: np.ndarray) -> dict[str, object]:
    activity = np.where(df["activity"].to_numpy() == "active", 1, 0)
    activity_table = pd.crosstab(df["activity"], labels).reindex(ACTIVITY_ORDER)
    activity_ari = adjusted_rand_score(activity, labels)
    activity_accuracy = balanced_accuracy(activity, labels)

    is_l1 = np.where(df["locus"].to_numpy() == "L1", 1, 0)
    locus_table = pd.crosstab(df["locus"], labels).reindex([loc for loc in LOCUS_ORDER if loc in df["locus"].unique()])
    locus_ari = adjusted_rand_score(df["locus"], labels)
    locus_accuracy = balanced_accuracy(is_l1, labels)

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
        activity_chi2, activity_p, _, _ = chi2_contingency(activity_table.to_numpy())
    except ValueError:
        activity_chi2, activity_p = float("nan"), float("nan")
    try:
        locus_chi2, locus_p, _, _ = chi2_contingency(locus_table.to_numpy())
    except ValueError:
        locus_chi2, locus_p = float("nan"), float("nan")

    return {
        "activity_ari": float(activity_ari),
        "activity_accuracy": activity_accuracy,
        "activity_chi2": float(activity_chi2),
        "activity_p_value": float(activity_p),
        "locus_ari": float(locus_ari),
        "locus_l1_accuracy": locus_accuracy,
        "locus_chi2": float(locus_chi2),
        "locus_p_value": float(locus_p),
        "silhouette": float(silhouette),
        "activity_table": activity_table,
        "locus_table": locus_table,
    }


def make_main_figure(df: pd.DataFrame, activity_verdict: str, locus_verdict: str) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(2, 3, figsize=(19, 11))
    locuses = [loc for loc in LOCUS_ORDER if loc in df["locus"].unique()]

    cluster_colors = {0: "#1f77b4", 1: "#ff7f0e"}
    point_colors = df["cluster"].map(cluster_colors)
    locus_colors = df["locus"].map(LOCUS_COLORS)
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
    ax.scatter(df["length"], df["mutations"], c=locus_colors, s=14, alpha=0.55, linewidths=0)
    ax.set_ylim(0, mutation_cap)
    ax.set_xlabel("Sequence length (bp)")
    ax.set_ylabel("Mutations vs consensus")
    ax.set_title("Same data coloured by locus (L)")
    ax.legend(
        handles=[Patch(color=LOCUS_COLORS[loc], label=loc) for loc in locuses],
        loc="upper left",
    )

    ax = axes[0, 2]
    sns.boxplot(
        data=df,
        x="locus",
        y="length",
        hue="activity",
        order=locuses,
        hue_order=ACTIVITY_ORDER,
        palette=ACTIVITY_COLORS,
        dodge=False,
        ax=ax,
        fliersize=2,
    )
    ax.set_xlabel("Locus")
    ax.set_ylabel("Sequence length (bp)")
    ax.set_title("Length by locus")
    if ax.legend_ is not None:
        ax.legend_.remove()

    ax = axes[1, 0]
    rate_cap = df["mutation_rate"].quantile(0.995) * 1.05
    sns.boxplot(
        data=df,
        x="locus",
        y="mutation_rate",
        hue="activity",
        order=locuses,
        hue_order=ACTIVITY_ORDER,
        palette=ACTIVITY_COLORS,
        dodge=False,
        ax=ax,
        fliersize=2,
    )
    ax.set_ylim(0, rate_cap)
    ax.set_xlabel("Locus")
    ax.set_ylabel("Mutations per bp")
    ax.set_title("Mutation rate by locus")
    ax.legend(
        handles=[Patch(color=ACTIVITY_COLORS[label], label=label) for label in ACTIVITY_ORDER],
        loc="upper left",
    )

    ax = axes[1, 1]
    composition = pd.crosstab(df["locus"], df["cluster"]).reindex(locuses)
    sns.heatmap(composition, annot=True, fmt="d", cmap="Blues", ax=ax, cbar_kws={"label": "sequences"})
    ax.set_xlabel("K-means cluster")
    ax.set_ylabel("Locus")
    ax.set_title("Cluster composition by locus")

    ax = axes[1, 2]
    locus_composition = (
        pd.crosstab(df["cluster"], df["locus"], normalize="index")
        .reindex(columns=locuses)
        .fillna(0.0)
    )
    locus_composition.plot(
        kind="bar",
        stacked=True,
        color=[LOCUS_COLORS[loc] for loc in locuses],
        ax=ax,
        rot=0,
        legend=False,
    )
    for position, (_, row) in enumerate(locus_composition.iterrows()):
        baseline = 0.0
        for loc in locuses:
            value = row[loc]
            if value > 0.05:
                ax.text(
                    position,
                    baseline + value / 2,
                    f"{value:.0%}",
                    ha="center",
                    va="center",
                    color="white",
                    fontweight="bold",
                )
            baseline += value
    ax.set_xlabel("K-means cluster")
    ax.set_ylabel("Fraction of sequences")
    ax.set_ylim(0, 1)
    ax.set_title("Locus composition per cluster")
    ax.legend(handles=[Patch(color=LOCUS_COLORS[loc], label=loc) for loc in locuses], loc="upper right")

    fig.suptitle(
        f"Rca3_R3 centromeric repeat clustering by locus (L)\n"
        f"AGENTS activity hypothesis {activity_verdict} | locus structure {locus_verdict}",
        fontsize=20,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUTPUT_DIR / "cluster_analysis_L.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def make_distribution_figure(df: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 3, figsize=(19, 6))
    locuses = [loc for loc in LOCUS_ORDER if loc in df["locus"].unique()]
    rate_cap = df["mutation_rate"].quantile(0.995) * 1.05
    monomer_cap = df["mutations_per_monomer"].quantile(0.995) * 1.05

    ax = axes[0]
    sns.violinplot(
        data=df,
        x="locus",
        y="mutation_rate",
        hue="activity",
        order=locuses,
        hue_order=ACTIVITY_ORDER,
        palette=ACTIVITY_COLORS,
        dodge=False,
        cut=0,
        ax=ax,
    )
    ax.set_ylim(0, rate_cap)
    ax.set_xlabel("Locus")
    ax.set_ylabel("Mutations per bp")
    ax.set_title("Mutation-rate distributions by locus")
    if ax.legend_ is not None:
        ax.legend_.remove()

    ax = axes[1]
    sns.violinplot(
        data=df,
        x="locus",
        y="mutations_per_monomer",
        hue="activity",
        order=locuses,
        hue_order=ACTIVITY_ORDER,
        palette=ACTIVITY_COLORS,
        dodge=False,
        cut=0,
        ax=ax,
    )
    ax.set_ylim(0, monomer_cap)
    ax.set_xlabel("Locus")
    ax.set_ylabel("Mutations per repeat unit")
    ax.set_title("Mutations per repeat unit by locus")
    ax.legend(
        handles=[Patch(color=ACTIVITY_COLORS[label], label=label) for label in ACTIVITY_ORDER],
        loc="upper left",
    )

    ax = axes[2]
    for loc in locuses:
        values = np.sort(df.loc[df["locus"] == loc, "mutation_rate"].to_numpy())
        if values.size:
            cumulative = np.arange(1, values.size + 1) / values.size
            ax.plot(values, cumulative, color=LOCUS_COLORS[loc], label=loc, linewidth=2)
    ax.set_xlim(0, rate_cap)
    ax.set_xlabel("Mutations per bp")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("Empirical cumulative distribution")
    ax.legend()

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "mutation_distributions_L.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not FASTA_PATH.exists():
        raise FileNotFoundError(f"Input FASTA not found: {FASTA_PATH}")

    df, skipped, total_records = read_fasta(FASTA_PATH)
    locuses = [loc for loc in LOCUS_ORDER if loc in df["locus"].unique()]
    print(f"Loaded {len(df)} sequences ({skipped} skipped of {total_records} records)")
    print("Sequences per locus:")
    print(df["locus"].value_counts().reindex(locuses).to_string())
    print("\nRepeat composition within each locus:")
    print(pd.crosstab(df["locus"], df["repeat"]).reindex(locuses).reindex(columns=REPEAT_ORDER).fillna(0).astype(int).to_string())

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

    locus_consensuses = {}
    print("\nLocus-specific consensus sequences:")
    for loc in locuses:
        mask = (df["locus"] == loc).to_numpy()
        group_sequences = df.loc[mask, "sequence"].tolist()
        group_encoded = [encoded for encoded, keep in zip(encoded_sequences, mask) if keep]
        locus_consensus, locus_string, locus_iterations = build_consensus(
            group_encoded, group_sequences, period
        )
        locus_consensuses[loc] = locus_consensus
        print(f"  {loc}: n={int(mask.sum())}, built in {locus_iterations} iteration(s)")
        print(f"    {locus_string}")
    df = annotate_locus_mutations(df, encoded_sequences, locus_consensuses)

    print("\nSequence statistics by locus (global consensus):")
    summary = (
        df.groupby("locus")[
            [
                "length",
                "mutations",
                "mutation_rate",
                "mutations_per_monomer",
                "mutations_locus",
                "mutation_rate_locus",
                "mutations_per_monomer_locus",
            ]
        ]
        .agg(["median", "mean"])
        .reindex(locuses)
    )
    print(summary.round(4).to_string())

    print("\nLength and mutation summary by locus:")
    overview = df.groupby("locus").agg(
        n=("length", "size"),
        length_median=("length", "median"),
        length_iqr=("length", lambda values: values.quantile(0.75) - values.quantile(0.25)),
        mutations_median=("mutations", "median"),
        mutation_rate_median=("mutation_rate", "median"),
        mutation_rate_mean=("mutation_rate", "mean"),
        mutations_locus_median=("mutations_locus", "median"),
        mutation_rate_locus_median=("mutation_rate_locus", "median"),
    )
    print(overview.reindex(locuses).round(4).to_string())

    labels, scaled, model = fit_kmeans(df, ["length", "mutations"])
    df["cluster"] = labels

    print("\nK-means (k=2) on length + mutations, cluster sizes:")
    print(pd.Series(labels).value_counts().sort_index().to_string())
    print("\nLocus by cluster:")
    print(pd.crosstab(df["locus"], df["cluster"]).reindex(locuses).to_string())
    print("\nRepeat by cluster:")
    print(pd.crosstab(df["repeat"], df["cluster"]).reindex(REPEAT_ORDER).to_string())
    print("\nActivity by cluster:")
    activity_table = pd.crosstab(df["activity"], df["cluster"]).reindex(ACTIVITY_ORDER)
    print(activity_table.to_string())

    metrics = cluster_summary(df, labels)
    print(f"\nAdjusted Rand index (cluster vs activity): {metrics['activity_ari']:.4f}")
    print(f"Balanced accuracy of best cluster-to-activity mapping: {metrics['activity_accuracy']:.4f}")
    print(
        f"Chi-square test of cluster vs activity: chi2={metrics['activity_chi2']:.2f}, "
        f"p={metrics['activity_p_value']:.3e}"
    )
    print(f"\nAdjusted Rand index (cluster vs locus): {metrics['locus_ari']:.4f}")
    print(f"Balanced accuracy separating L1 from L2/L3: {metrics['locus_l1_accuracy']:.4f}")
    print(
        f"Chi-square test of cluster vs locus: chi2={metrics['locus_chi2']:.2f}, "
        f"p={metrics['locus_p_value']:.3e}"
    )
    print(f"Silhouette score: {metrics['silhouette']:.4f}")

    activity_verdict = (
        "CONFIRMED"
        if metrics["activity_accuracy"] >= 0.85 and metrics["activity_ari"] >= 0.5
        else "REJECTED"
    )
    locus_verdict = (
        "CONFIRMED"
        if metrics["locus_l1_accuracy"] >= 0.85 and metrics["locus_ari"] >= 0.5
        else "REJECTED"
    )

    print(f"\n>>> AGENTS HYPOTHESIS (active R1-R3 vs inactive R4-R5) {activity_verdict}")
    if activity_verdict == "REJECTED":
        print(
            "    The two clusters do not separate active (R1-R3) from inactive (R4-R5) repeats; "
            "both clusters contain a mixture of activity states."
        )
    else:
        print("    Active (R1-R3) and inactive (R4-R5) repeats occupy different clusters.")

    print(f"\n>>> LOCUS (L) STRUCTURE HYPOTHESIS {locus_verdict}")
    if locus_verdict == "CONFIRMED":
        print(
            "    Sequences separate by locus: L1 (single repeat unit, R2+R5) forms one cluster "
            "while L2 (R3+R4) and L3 (R1) form the other."
        )
    else:
        print("    Locus groups are not cleanly separated by the two-cluster solution.")

    print("\nSensitivity of the two-cluster solution to feature choice (grouped by locus):")
    feature_sets = {
        "length + mutations": ["length", "mutations"],
        "length + mutation_rate": ["length", "mutation_rate"],
        "mutations only": ["mutations"],
        "mutation_rate only": ["mutation_rate"],
        "length only": ["length"],
        "length + mutation_rate_locus": ["length", "mutation_rate_locus"],
        "mutation_rate_locus only": ["mutation_rate_locus"],
    }
    sensitivity_rows = []
    for name, features in feature_sets.items():
        feature_labels, feature_scaled, _ = fit_kmeans(df, features)
        activity = np.where(df["activity"].to_numpy() == "active", 1, 0)
        is_l1 = np.where(df["locus"].to_numpy() == "L1", 1, 0)
        activity_ari = adjusted_rand_score(activity, feature_labels)
        activity_accuracy = balanced_accuracy(activity, feature_labels)
        locus_ari = adjusted_rand_score(df["locus"], feature_labels)
        locus_accuracy = balanced_accuracy(is_l1, feature_labels)
        clusters = pd.crosstab(df["locus"], feature_labels).reindex(locuses)
        active_fraction = (
            df.assign(cluster=feature_labels)
            .groupby("cluster")["activity"]
            .apply(lambda values: float((values == "active").mean()))
        )
        dominant_loci = {}
        for cluster in [0, 1]:
            column = clusters.get(cluster, pd.Series(0, index=clusters.index))
            dominant_loci[cluster] = f"{column.idxmax()} ({column.max() / column.sum():.0%})" if column.sum() else "NA"
        sensitivity_rows.append(
            {
                "features": name,
                "activity_ari": activity_ari,
                "activity_balanced_accuracy": activity_accuracy,
                "locus_ari": locus_ari,
                "locus_l1_balanced_accuracy": locus_accuracy,
                "active_fraction_cluster_0": active_fraction.get(0, float("nan")),
                "active_fraction_cluster_1": active_fraction.get(1, float("nan")),
                "dominant_locus_cluster_0": dominant_loci[0],
                "dominant_locus_cluster_1": dominant_loci[1],
            }
        )
    sensitivity = pd.DataFrame(sensitivity_rows)
    print(sensitivity.round(4).to_string(index=False))
    sensitivity.to_csv(OUTPUT_DIR / "clustering_sensitivity_L.csv", index=False)

    labels3, _, _ = fit_kmeans(df, ["length", "mutations"], n_clusters=3)
    print("\nLength + mutation structure with k=3 (sanity check for locus/length classes):")
    print(pd.crosstab(df["locus"], labels3).reindex(locuses).to_string())

    print("\nLocus-matched comparisons of mutation rate (Mann-Whitney U):")
    for left, right in [("L1", "L2"), ("L2", "L3"), ("L1", "L3")]:
        left_values = df.loc[df["locus"] == left, "mutation_rate"]
        right_values = df.loc[df["locus"] == right, "mutation_rate"]
        statistic, p_value = mannwhitneyu(left_values, right_values, alternative="two-sided")
        print(
            f"  {left} vs {right}: median {left_values.median():.4f} vs {right_values.median():.4f}, "
            f"mean {left_values.mean():.4f} vs {right_values.mean():.4f}, p={p_value:.3e}"
        )

    statistic, p_value = kruskal(*[df.loc[df["locus"] == loc, "mutation_rate"] for loc in locuses])
    print(f"\nKruskal-Wallis across loci: H={statistic:.2f}, p={p_value:.3e}")

    print("\nRepeat-matched comparisons of mutation rate (Mann-Whitney U):")
    for left, right in [("R2", "R5"), ("R3", "R4")]:
        left_values = df.loc[df["repeat"] == left, "mutation_rate"]
        right_values = df.loc[df["repeat"] == right, "mutation_rate"]
        statistic, p_value = mannwhitneyu(left_values, right_values, alternative="two-sided")
        print(
            f"  {left} vs {right}: median {left_values.median():.4f} vs {right_values.median():.4f}, "
            f"p={p_value:.3e}"
        )

    make_main_figure(df, activity_verdict, locus_verdict)
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
        "mutations_locus",
        "mutation_rate_locus",
        "mutations_per_monomer_locus",
        "cluster",
        "activity",
        "sequence",
    ]
    df[output_columns].to_csv(OUTPUT_DIR / "sequence_analysis_results_L.csv", index=False)

    print("\nSaved outputs:")
    print(f"  {OUTPUT_DIR / 'cluster_analysis_L.png'}")
    print(f"  {OUTPUT_DIR / 'mutation_distributions_L.png'}")
    print(f"  {OUTPUT_DIR / 'sequence_analysis_results_L.csv'}")
    print(f"  {OUTPUT_DIR / 'clustering_sensitivity_L.csv'}")


if __name__ == "__main__":
    main()
