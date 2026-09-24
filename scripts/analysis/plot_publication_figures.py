#!/usr/bin/env python3
"""Create publication figures from terminal-current scores and summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MultipleLocator


GROUP_ORDER = ["mock", "MHV-wt", "MHV-mut"]
VIRUS_GROUP_ORDER = ["MHV-wt", "MHV-mut"]
GROUP_LABEL = {"mock": "Mock", "MHV-wt": "WT", "MHV-mut": "Mut"}
GROUP_COLOR = {"mock": "#7A7A7A", "MHV-wt": "#0072B2", "MHV-mut": "#D55E00"}
BARCODE_ORDER = ["bc04", "bc05", "bc06", "bc07", "bc11"]
INFECTED_BARCODES = ["bc05", "bc06", "bc07", "bc11"]
BARCODE_GROUP = {
    "bc04": "mock",
    "bc05": "MHV-wt",
    "bc06": "MHV-wt",
    "bc07": "MHV-mut",
    "bc11": "MHV-mut",
}
GENE_ORDER = ["sgRNA4_ORF4", "sgRNA5_ORF5_E", "sgRNA6_M", "sgRNA7_N"]
GENE_LABEL = {
    "sgRNA4_ORF4": "ORF4",
    "sgRNA5_ORF5_E": "ORF5/E",
    "sgRNA6_M": "M",
    "sgRNA7_N": "N",
}
GENE_COLOR = {
    "sgRNA4_ORF4": "#CC79A7",
    "sgRNA5_ORF5_E": "#E69F00",
    "sgRNA6_M": "#009E73",
    "sgRNA7_N": "#56B4E9",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--package-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Publication package root (default: parent of scripts/).",
    )
    parser.add_argument(
        "--formats",
        default="png,pdf,svg",
        help="Comma-separated output formats (default: png,pdf,svg).",
    )
    parser.add_argument("--dpi", type=int, default=400, help="PNG resolution.")
    return parser.parse_args()


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 9.5,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8,
            "legend.frameon": False,
            "figure.dpi": 120,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def polish_axis(ax: plt.Axes, grid_axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis=grid_axis, color="#D9D9D9", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )


def save_figure(
    fig: plt.Figure, output_dir: Path, stem: str, formats: list[str], dpi: int
) -> None:
    for extension in formats:
        kwargs = {"dpi": dpi} if extension == "png" else {}
        fig.savefig(output_dir / f"{stem}.{extension}", **kwargs)
    plt.close(fig)


def figure_1(
    package_root: Path, output_dir: Path, data_dir: Path, formats: list[str], dpi: int
) -> None:
    primary_dir = package_root / "results" / "final_primary"
    frames = []
    for origin, filename in [
        ("MHV-mapped", "mhv_reads_terminal_nona_direct.parquet"),
        ("Host-mapped", "host_terminal_nona_direct_scores.parquet"),
    ]:
        frame = pd.read_parquet(
            primary_dir / filename,
            columns=[
                "mapping_barcode",
                "sample_group",
                "eligible_terminal_nona_direct",
                "terminal_nona_signed_difference_pa",
            ],
        )
        frame = frame[frame["eligible_terminal_nona_direct"].eq(1)].copy()
        if origin == "MHV-mapped":
            frame = frame[frame["mapping_barcode"].isin(INFECTED_BARCODES)].copy()
        frame["read_origin"] = origin
        frames.append(frame)
    signed = pd.concat(frames, ignore_index=True)

    bin_edges = np.arange(-15.0, 15.0001, 0.25)
    distribution_rows = []
    summary_rows = []
    for (barcode, origin), subset in signed.groupby(
        ["mapping_barcode", "read_origin"], sort=False
    ):
        values = subset["terminal_nona_signed_difference_pa"].to_numpy(dtype=float)
        in_range = values[(values >= bin_edges[0]) & (values <= bin_edges[-1])]
        counts, _ = np.histogram(in_range, bins=bin_edges)
        for left, right, count in zip(bin_edges[:-1], bin_edges[1:], counts):
            distribution_rows.append(
                {
                    "mapping_barcode": barcode,
                    "sample_group": BARCODE_GROUP[barcode],
                    "read_origin": origin,
                    "bin_left_pa": left,
                    "bin_right_pa": right,
                    "bin_mid_pa": (left + right) / 2,
                    "n_reads": int(count),
                    "percent_of_all_callable": count / len(values) * 100,
                }
            )
        q25, median, q75 = np.quantile(values, [0.25, 0.5, 0.75])
        n_negative = int((values <= -5.0).sum())
        n_positive = int((values >= 5.0).sum())
        summary_rows.append(
            {
                "mapping_barcode": barcode,
                "sample_group": BARCODE_GROUP[barcode],
                "read_origin": origin,
                "n_callable": len(values),
                "q25_signed_deviation_pa": q25,
                "median_signed_deviation_pa": median,
                "q75_signed_deviation_pa": q75,
                "iqr_signed_deviation_pa": q75 - q25,
                "n_negative_tail": n_negative,
                "negative_tail_fraction": n_negative / len(values),
                "n_positive_tail": n_positive,
                "positive_tail_fraction": n_positive / len(values),
                "n_total_beyond_5pa": n_negative + n_positive,
                "total_beyond_5pa_fraction": (n_negative + n_positive) / len(values),
                "n_below_plot_range": int((values < bin_edges[0]).sum()),
                "n_above_plot_range": int((values > bin_edges[-1]).sum()),
            }
        )
    distribution = pd.DataFrame(distribution_rows)
    summary = pd.DataFrame(summary_rows)
    distribution.to_csv(
        data_dir / "figure1_current_distribution.tsv", sep="\t", index=False
    )
    summary.to_csv(
        data_dir / "figure1_current_distribution_summary.tsv", sep="\t", index=False
    )

    fig = plt.figure(figsize=(13.0, 7.45), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 1])
    axes = [fig.add_subplot(grid[row, column]) for row in range(2) for column in range(3)]
    origin_style = {
        "Host-mapped": {"color": "#6F777C", "label": "Host", "linewidth": 1.45},
        "MHV-mapped": {"color": None, "label": "MHV", "linewidth": 1.8},
    }
    ymax = distribution["percent_of_all_callable"].max() * 1.32

    for panel_index, (barcode, ax) in enumerate(zip(BARCODE_ORDER, axes[:5])):
        group = BARCODE_GROUP[barcode]
        ax.axvspan(-15, -5, color="#EAF3F7", alpha=0.9, zorder=0)
        ax.axvspan(-5, 5, color="#F5F6F7", alpha=0.95, zorder=0)
        ax.axvspan(5, 15, color="#FFF0E3", alpha=0.9, zorder=0)
        barcode_summary = summary[summary["mapping_barcode"] == barcode]
        annotation_lines = []
        for origin in ["Host-mapped", "MHV-mapped"]:
            curve = distribution[
                (distribution["mapping_barcode"] == barcode)
                & (distribution["read_origin"] == origin)
            ]
            if curve.empty:
                continue
            style = origin_style[origin].copy()
            color = style["color"] or GROUP_COLOR[group]
            ax.stairs(
                curve["percent_of_all_callable"],
                bin_edges,
                color=color,
                linewidth=style["linewidth"],
                fill=False,
                label=style["label"],
                zorder=3,
            )
            row = barcode_summary[barcode_summary["read_origin"] == origin].iloc[0]
            annotation_lines.append(
                f"{style['label']}: n={int(row['n_callable']):,}; "
                f"median {row['median_signed_deviation_pa']:+.2f}; "
                f"|deviation|>=5: {row['total_beyond_5pa_fraction'] * 100:.2f}%"
            )
        ax.axvline(0, color="#242424", linewidth=0.9, zorder=2)
        ax.axvline(-5, color="#2E6F95", linestyle="--", linewidth=1.0, zorder=2)
        ax.axvline(5, color="#D55E00", linestyle="--", linewidth=1.0, zorder=2)
        ax.text(
            0.02,
            0.97,
            "\n".join(annotation_lines),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.25,
            color="#27333B",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 2.6},
        )
        title_suffix = "Mock (host only)" if barcode == "bc04" else GROUP_LABEL[group]
        ax.set_title(f"{barcode} · {title_suffix}")
        ax.set_xlim(-15, 15)
        ax.set_ylim(0, ymax)
        ax.set_xticks([-15, -10, -5, 0, 5, 10, 15])
        if panel_index in [0, 3]:
            ax.set_ylabel("Callable reads per 0.25-pA bin (%)")
        if panel_index >= 2:
            ax.set_xlabel("Signed deviation from host reference (pA)")
        if panel_index == 0:
            ax.legend(
                handles=[
                    Line2D([], [], color="#6F777C", linewidth=1.45, label="Host-mapped"),
                    Line2D([], [], color=GROUP_COLOR["MHV-wt"], linewidth=1.8, label="MHV-mapped (WT)"),
                    Line2D([], [], color=GROUP_COLOR["MHV-mut"], linewidth=1.8, label="MHV-mapped (Mut)"),
                ],
                loc="lower right",
                fontsize=7.2,
            )
        polish_axis(ax)
        panel_label(ax, chr(ord("A") + panel_index))

    scale_ax = axes[5]
    scale_ax.set_xlim(-15, 15)
    scale_ax.set_ylim(0, 1)
    scale_ax.axvspan(-15, -5, ymin=0.34, ymax=0.62, color="#EAF3F7")
    scale_ax.axvspan(-5, 5, ymin=0.34, ymax=0.62, color="#F0F2F3")
    scale_ax.axvspan(5, 15, ymin=0.34, ymax=0.62, color="#FFF0E3")
    scale_ax.axvline(-5, ymin=0.31, ymax=0.66, color="#2E6F95", linestyle="--", linewidth=1.2)
    scale_ax.axvline(0, ymin=0.31, ymax=0.66, color="#242424", linewidth=1.0)
    scale_ax.axvline(5, ymin=0.31, ymax=0.66, color="#D55E00", linestyle="--", linewidth=1.2)
    scale_ax.annotate(
        "5 pA",
        xy=(5, 0.70),
        xytext=(0, 0.70),
        arrowprops={"arrowstyle": "<->", "color": "#242424", "linewidth": 1.0},
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )
    scale_ax.text(-10, 0.48, "negative tail", ha="center", va="center", fontsize=8)
    scale_ax.text(0, 0.48, "reference-like", ha="center", va="center", fontsize=8)
    scale_ax.text(10, 0.48, "positive tail", ha="center", va="center", fontsize=8)
    scale_ax.text(
        0.02,
        0.97,
        "How to read the threshold",
        transform=scale_ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.2,
        fontweight="bold",
        color="#14324A",
    )
    scale_ax.text(
        0.02,
        0.24,
        (
            "5 pA = 5 × 10⁻¹² A. It is applied to the difference between two "
            "31-sample window means (7.75 ms each), after centering on the "
            "−10.19-pA host contrast. The equivalent raw-contrast boundaries "
            "are −15.19 and −5.19 pA."
        ),
        transform=scale_ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.1,
        color="#3D4B55",
        wrap=True,
    )
    scale_ax.set_xticks([-15, -10, -5, 0, 5, 10, 15])
    scale_ax.set_xlabel("Signed deviation from host reference (pA)")
    scale_ax.set_yticks([])
    for spine in ["left", "right", "top"]:
        scale_ax.spines[spine].set_visible(False)
    panel_label(scale_ax, "F")

    fig.suptitle(
        "Figure 1. Per-sample distributions of terminal-junction current deviation",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.022,
        (
            "Histograms share 0.25-pA bins and are normalized within origin and barcode. "
            "The displayed -15 to +15 pA range contains more than 99.9% of reads; "
            "tail percentages use all callable reads."
        ),
        ha="center",
        va="top",
        fontsize=8.2,
        color="#4A4A4A",
    )
    save_figure(fig, output_dir, "figure1_sample_current_distributions", formats, dpi)


def figure_2(
    package_root: Path, output_dir: Path, data_dir: Path, formats: list[str], dpi: int
) -> None:
    result_dir = package_root / "results" / "final_host_virus_sgrna"
    by_barcode = pd.read_csv(result_dir / "virus_high_conf_sgrna_by_barcode_summary.tsv", sep="\t")
    by_group = pd.read_csv(result_dir / "virus_high_conf_sgrna_by_group_summary.tsv", sep="\t")
    by_barcode = by_barcode[
        by_barcode["sample_group"].isin(["MHV-wt", "MHV-mut"])
        & by_barcode["sgrna_primary_assignment"].isin(GENE_ORDER)
    ].copy()
    by_group = by_group[
        by_group["sample_group"].isin(["MHV-wt", "MHV-mut"])
        & by_group["sgrna_primary_assignment"].isin(GENE_ORDER)
    ].copy()
    by_barcode["gene_label"] = by_barcode["sgrna_primary_assignment"].map(GENE_LABEL)
    by_group["gene_label"] = by_group["sgrna_primary_assignment"].map(GENE_LABEL)
    by_group["mut_minus_wt_pp"] = by_group["sgrna_primary_assignment"].map(
        by_group.pivot(
            index="sgrna_primary_assignment",
            columns="sample_group",
            values="nonA_like_fraction_callable",
        ).assign(delta=lambda frame: (frame["MHV-mut"] - frame["MHV-wt"]) * 100)["delta"]
    )
    by_barcode.to_csv(data_dir / "figure2_sgrna_by_barcode.tsv", sep="\t", index=False)
    by_group.to_csv(data_dir / "figure2_sgrna_by_group.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(13.0, 4.35), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.28, 0.72])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])

    infected_barcodes = ["bc05", "bc06", "bc07", "bc11"]
    y = np.arange(len(infected_barcodes))
    left = np.zeros(len(infected_barcodes), dtype=float)
    for gene in GENE_ORDER:
        counts = (
            by_barcode[by_barcode["sgrna_primary_assignment"] == gene]
            .set_index("mapping_barcode")
            .reindex(infected_barcodes)["n_callable"]
            .fillna(0)
            .to_numpy()
        )
        ax_a.barh(
            y,
            counts,
            left=left,
            height=0.62,
            color=GENE_COLOR[gene],
            edgecolor="white",
            linewidth=0.5,
            label=GENE_LABEL[gene],
        )
        left += counts
    for ypos, total in zip(y, left):
        ax_a.text(total + max(left) * 0.015, ypos, f"n={int(total):,}", va="center", fontsize=8)
    barcode_labels = [
        f"{barcode}  ({GROUP_LABEL[BARCODE_GROUP[barcode]]})" for barcode in infected_barcodes
    ]
    ax_a.set_yticks(y, barcode_labels)
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, max(left) * 1.19)
    ax_a.set_xlabel("Callable high-confidence sgRNA reads")
    ax_a.set_title("Transcript-resolved read depth")
    ax_a.legend(ncol=2, loc="lower right")
    polish_axis(ax_a, grid_axis="x")
    panel_label(ax_a, "A")

    gene_y = np.arange(len(GENE_ORDER))
    offsets = {"MHV-wt": -0.13, "MHV-mut": 0.13}
    for group in ["MHV-wt", "MHV-mut"]:
        for gene_index, gene in enumerate(GENE_ORDER):
            row = by_group[
                (by_group["sample_group"] == group)
                & (by_group["sgrna_primary_assignment"] == gene)
            ].iloc[0]
            value = row["nonA_like_fraction_callable"] * 100
            low = row["wilson95_low_descriptive"] * 100
            high = row["wilson95_high_descriptive"] * 100
            ypos = gene_index + offsets[group]
            ax_b.errorbar(
                value,
                ypos,
                xerr=[[value - low], [high - value]],
                fmt="D",
                markersize=6,
                color=GROUP_COLOR[group],
                markeredgecolor="white",
                markeredgewidth=0.5,
                capsize=2.5,
                elinewidth=1.3,
                zorder=3,
            )
            barcode_rows = by_barcode[
                (by_barcode["sample_group"] == group)
                & (by_barcode["sgrna_primary_assignment"] == gene)
            ]
            barcode_jitter = np.linspace(-0.04, 0.04, len(barcode_rows))
            for jitter, (_, barcode_row) in zip(barcode_jitter, barcode_rows.iterrows()):
                ax_b.scatter(
                    barcode_row["nonA_like_fraction_callable"] * 100,
                    ypos + jitter,
                    s=21,
                    facecolor="white",
                    edgecolor=GROUP_COLOR[group],
                    linewidth=1.0,
                    zorder=4,
                )
            ax_b.text(
                high + 0.7,
                ypos,
                f"{int(row['n_terminal_junction_nonA_like'])}/{int(row['n_callable'])}",
                va="center",
                fontsize=7.5,
                color=GROUP_COLOR[group],
            )
    ax_b.set_yticks(gene_y, [GENE_LABEL[gene] for gene in GENE_ORDER])
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 27)
    ax_b.set_xlabel("Terminal-junction non-A-like (%)")
    ax_b.set_title("WT and Mut estimates within each sgRNA")
    ax_b.legend(
        handles=[
            Line2D([], [], marker="D", linestyle="", color=GROUP_COLOR["MHV-wt"], label="WT pooled"),
            Line2D([], [], marker="D", linestyle="", color=GROUP_COLOR["MHV-mut"], label="Mut pooled"),
            Line2D([], [], marker="o", linestyle="", markerfacecolor="white", color="#555555", label="Barcode"),
        ],
        loc="lower right",
    )
    polish_axis(ax_b, grid_axis="x")
    panel_label(ax_b, "B")

    pivot = by_group.pivot(
        index="sgrna_primary_assignment",
        columns="sample_group",
        values="nonA_like_fraction_callable",
    ).reindex(GENE_ORDER)
    delta = (pivot["MHV-mut"] - pivot["MHV-wt"]) * 100
    bars = ax_c.barh(
        gene_y,
        delta,
        height=0.56,
        color=[GENE_COLOR[gene] for gene in GENE_ORDER],
        edgecolor="white",
    )
    for bar, value in zip(bars, delta):
        ax_c.text(value + 0.18, bar.get_y() + bar.get_height() / 2, f"+{value:.2f}", va="center", fontsize=8)
    ax_c.axvline(0, color="#333333", linewidth=0.8)
    ax_c.set_yticks(gene_y, [GENE_LABEL[gene] for gene in GENE_ORDER])
    ax_c.invert_yaxis()
    ax_c.set_xlim(0, max(delta) + 1.4)
    ax_c.set_xlabel("Mut − WT (percentage points)")
    ax_c.set_title("Descriptive difference")
    polish_axis(ax_c, grid_axis="x")
    panel_label(ax_c, "C")

    fig.suptitle(
        "Figure 2. Terminal-junction non-A-like signal across high-confidence MHV sgRNAs",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.035,
        "sgRNA2 and S are omitted because each genotype has ≤10 callable reads; no genotype p-values are shown.",
        ha="center",
        va="top",
        fontsize=8.2,
        color="#4A4A4A",
    )
    save_figure(fig, output_dir, "figure2_sgrna_wt_mut", formats, dpi)


def empirical_survival(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    values = np.sort(values[np.isfinite(values)])
    indices = np.searchsorted(values, grid, side="left")
    return (values.size - indices) / values.size


def figure_3(
    package_root: Path, output_dir: Path, data_dir: Path, formats: list[str], dpi: int
) -> None:
    primary_dir = package_root / "results" / "final_primary"
    reads = pd.read_parquet(
        primary_dir / "mhv_reads_terminal_nona_direct.parquet",
        columns=[
            "sample_group",
            "eligible_terminal_nona_direct",
            "terminal_nona_abs_difference_pa",
        ],
    )
    reads = reads[
        reads["eligible_terminal_nona_direct"].eq(1)
        & reads["sample_group"].isin(VIRUS_GROUP_ORDER)
    ].copy()
    lengths = pd.read_csv(primary_dir / "mhv_sgrna_nona_direct_length_summary.tsv", sep="\t")
    lengths = lengths[
        lengths["sample_group"].isin(["MHV-wt", "MHV-mut"])
        & lengths["sgrna_primary_assignment"].isin(GENE_ORDER)
        & lengths["terminal_nona_direct_class"].eq("terminal_junction_nonA_like")
    ].copy()
    lengths["gene_label"] = lengths["sgrna_primary_assignment"].map(GENE_LABEL)
    positive_reads = pd.read_csv(
        package_root
        / "results"
        / "final_host_virus_sgrna"
        / "virus_high_conf_sgrna_nonA_like_reads.tsv.gz",
        sep="\t",
        usecols=[
            "sample_group",
            "sgrna_primary_assignment",
            "polya_length_source",
        ],
    )
    positive_reads = positive_reads[
        positive_reads["sample_group"].isin(["MHV-wt", "MHV-mut"])
        & positive_reads["sgrna_primary_assignment"].isin(GENE_ORDER)
    ]
    source_mix = (
        positive_reads.groupby(["sample_group", "polya_length_source"], observed=True)
        .size()
        .rename("n_reads")
        .reset_index()
    )
    source_mix["fraction_within_group"] = source_mix["n_reads"] / source_mix.groupby(
        "sample_group"
    )["n_reads"].transform("sum")

    score_grid = np.linspace(0, 15, 601)
    curves = []
    for group in VIRUS_GROUP_ORDER:
        values = reads.loc[
            reads["sample_group"] == group, "terminal_nona_abs_difference_pa"
        ].to_numpy(dtype=float)
        survival = empirical_survival(values, score_grid)
        curves.append(
            pd.DataFrame(
                {
                    "sample_group": group,
                    "score_threshold_pa": score_grid,
                    "fraction_at_or_above": survival,
                    "n_callable": len(values),
                }
            )
        )
    curve_data = pd.concat(curves, ignore_index=True)
    curve_data.to_csv(data_dir / "figure3_signal_survival.tsv", sep="\t", index=False)
    lengths.to_csv(data_dir / "figure3_positive_length_summary.tsv", sep="\t", index=False)
    source_mix.to_csv(data_dir / "figure3_polya_length_source_mix.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(13.0, 4.4), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[1.18, 1.0, 1.08])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])

    for group in VIRUS_GROUP_ORDER:
        curve = curve_data[curve_data["sample_group"] == group]
        n_callable = int(curve["n_callable"].iloc[0])
        cutoff_fraction = reads.loc[
            reads["sample_group"] == group, "terminal_nona_abs_difference_pa"
        ].ge(5.0).mean()
        ax_a.plot(
            curve["score_threshold_pa"],
            curve["fraction_at_or_above"] * 100,
            color=GROUP_COLOR[group],
            linewidth=2.0,
            label=f"{GROUP_LABEL[group]}: {cutoff_fraction * 100:.2f}% (n={n_callable:,})",
        )
        ax_a.scatter(
            [5],
            [cutoff_fraction * 100],
            s=30,
            color=GROUP_COLOR[group],
            edgecolor="white",
            linewidth=0.5,
            zorder=4,
        )
    ax_a.axvline(5, color="#222222", linestyle="--", linewidth=1.1, label="Fixed cutoff: 5 pA")
    ax_a.set_xlim(0, 12)
    ax_a.set_ylim(0, 100)
    ax_a.set_xlabel("Absolute deviation from host A|DNA reference (pA)")
    ax_a.set_ylabel("Callable reads at or above threshold (%)")
    ax_a.set_title("Observed score distributions")
    ax_a.legend(loc="upper right")
    polish_axis(ax_a)
    panel_label(ax_a, "A")

    gene_y = np.arange(len(GENE_ORDER))
    offsets = {"MHV-wt": -0.13, "MHV-mut": 0.13}
    for group in ["MHV-wt", "MHV-mut"]:
        subset = lengths[lengths["sample_group"] == group].set_index("sgrna_primary_assignment")
        for gene_index, gene in enumerate(GENE_ORDER):
            row = subset.loc[gene]
            ypos = gene_index + offsets[group]
            median = row["polya_nt_median"]
            low = row["polya_nt_q25"]
            high = row["polya_nt_q75"]
            ax_b.errorbar(
                median,
                ypos,
                xerr=[[median - low], [high - median]],
                fmt="o",
                markersize=6,
                color=GROUP_COLOR[group],
                capsize=3,
                elinewidth=1.5,
                markeredgecolor="white",
                markeredgewidth=0.5,
            )
            ax_b.text(high + 1.2, ypos, f"n={int(row['n_reads'])}", va="center", fontsize=7.4)
    ax_b.set_yticks(gene_y, [GENE_LABEL[gene] for gene in GENE_ORDER])
    ax_b.invert_yaxis()
    ax_b.set_xlim(25, 95)
    ax_b.set_xlabel("Poly(A)-tail length (nt), median and IQR")
    ax_b.set_title("Positive-read poly(A) lengths*")
    ax_b.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="", color=GROUP_COLOR["MHV-wt"], label="WT"),
            Line2D([], [], marker="o", linestyle="", color=GROUP_COLOR["MHV-mut"], label="Mut"),
        ],
        loc="lower right",
    )
    polish_axis(ax_b, grid_axis="x")
    panel_label(ax_b, "B")

    for group in ["MHV-wt", "MHV-mut"]:
        subset = lengths[lengths["sample_group"] == group].set_index("sgrna_primary_assignment")
        for gene_index, gene in enumerate(GENE_ORDER):
            row = subset.loc[gene]
            ypos = gene_index + offsets[group]
            median = row["read_nt_median"]
            low = row["read_nt_q25"]
            high = row["read_nt_q75"]
            ax_c.errorbar(
                median,
                ypos,
                xerr=[[median - low], [high - median]],
                fmt="o",
                markersize=6,
                color=GROUP_COLOR[group],
                capsize=3,
                elinewidth=1.5,
                markeredgecolor="white",
                markeredgewidth=0.5,
            )
            ax_c.text(high + 45, ypos, f"n={int(row['n_reads'])}", va="center", fontsize=7.4)
    ax_c.set_yticks(gene_y, [GENE_LABEL[gene] for gene in GENE_ORDER])
    ax_c.invert_yaxis()
    ax_c.set_xlim(1500, 3750)
    ax_c.xaxis.set_major_locator(MultipleLocator(500))
    ax_c.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1000:.1f}k"))
    ax_c.set_xlabel("Read length (nt), median and IQR")
    ax_c.set_title("Positive-read lengths by sgRNA")
    polish_axis(ax_c, grid_axis="x")
    panel_label(ax_c, "C")

    fig.suptitle(
        "Figure 3. Signal-threshold behavior and lengths of positive sgRNA reads",
        fontsize=13,
        fontweight="bold",
    )
    wt_dorado = source_mix.loc[
        (source_mix["sample_group"] == "MHV-wt")
        & (source_mix["polya_length_source"] == "dorado_pt"),
        "fraction_within_group",
    ].iloc[0]
    mut_dorado = source_mix.loc[
        (source_mix["sample_group"] == "MHV-mut")
        & (source_mix["polya_length_source"] == "dorado_pt"),
        "fraction_within_group",
    ].iloc[0]
    fig.text(
        0.5,
        -0.035,
        (
            f"Positive-read poly(A) source mix: WT {wt_dorado * 100:.1f}% Dorado; "
            f"Mut {mut_dorado * 100:.1f}% Dorado (remainder: WarpDemuX fallback). "
            "Do not interpret the pooled length contrast as genotype-specific."
        ),
        ha="center",
        va="top",
        fontsize=8.2,
        color="#4A4A4A",
    )
    save_figure(fig, output_dir, "figure3_signal_and_positive_read_lengths", formats, dpi)


def figure_4(
    package_root: Path, output_dir: Path, data_dir: Path, formats: list[str], dpi: int
) -> None:
    primary_dir = package_root / "results" / "final_primary"
    frames = []
    for origin, filename in [
        ("MHV-mapped", "mhv_reads_terminal_nona_direct.parquet"),
        ("Host-mapped", "host_terminal_nona_direct_scores.parquet"),
    ]:
        frame = pd.read_parquet(
            primary_dir / filename,
            columns=[
                "mapping_barcode",
                "eligible_terminal_nona_direct",
                "terminal_nona_signed_difference_pa",
            ],
        )
        frame = frame[frame["eligible_terminal_nona_direct"].eq(1)].copy()
        if origin == "MHV-mapped":
            frame = frame[frame["mapping_barcode"].isin(INFECTED_BARCODES)].copy()
        frame["read_origin"] = origin
        frames.append(frame)
    signed = pd.concat(frames, ignore_index=True)

    polarity_rows = []
    for (barcode, origin), subset in signed.groupby(["mapping_barcode", "read_origin"]):
        values = subset["terminal_nona_signed_difference_pa"].to_numpy(dtype=float)
        n_callable = len(values)
        n_negative = int((values <= -5.0).sum())
        n_positive = int((values >= 5.0).sum())
        polarity_rows.append(
            {
                "mapping_barcode": barcode,
                "sample_group": BARCODE_GROUP[barcode],
                "read_origin": origin,
                "n_callable": n_callable,
                "n_negative_direction_nonA_like": n_negative,
                "n_within_cutoff": n_callable - n_negative - n_positive,
                "n_positive_direction_nonA_like": n_positive,
                "negative_direction_fraction_callable": n_negative / n_callable,
                "positive_direction_fraction_callable": n_positive / n_callable,
                "total_nonA_like_fraction_callable": (n_negative + n_positive) / n_callable,
            }
        )
    polarity = pd.DataFrame(polarity_rows)
    polarity["barcode_order"] = polarity["mapping_barcode"].map(
        {barcode: index for index, barcode in enumerate(BARCODE_ORDER)}
    )
    polarity["origin_order"] = polarity["read_origin"].map(
        {"Host-mapped": 0, "MHV-mapped": 1}
    )
    polarity = polarity.sort_values(["barcode_order", "origin_order"]).reset_index(drop=True)
    polarity.to_csv(data_dir / "figure4_signed_polarity.tsv", sep="\t", index=False)

    sensitivity = pd.read_csv(
        primary_dir / "mhv_barcode_nona_direct_threshold_sensitivity.tsv", sep="\t"
    )
    sensitivity = sensitivity[
        sensitivity["mapping_barcode"].isin(INFECTED_BARCODES)
    ].copy()
    sensitivity["barcode_order"] = sensitivity["mapping_barcode"].map(
        {barcode: index for index, barcode in enumerate(INFECTED_BARCODES)}
    )
    sensitivity = sensitivity.sort_values(["cutoff_pa", "barcode_order"]).reset_index(drop=True)
    pooled = (
        sensitivity.groupby(["cutoff_pa", "sample_group"], observed=True)[
            ["n_callable", "n_terminal_junction_nonA_like"]
        ]
        .sum()
        .reset_index()
    )
    pooled["terminal_junction_nonA_like_fraction_callable"] = (
        pooled["n_terminal_junction_nonA_like"] / pooled["n_callable"]
    )
    difference = pooled.pivot(
        index="cutoff_pa",
        columns="sample_group",
        values="terminal_junction_nonA_like_fraction_callable",
    ).reset_index()
    difference["mut_minus_wt_pp"] = (
        difference["MHV-mut"] - difference["MHV-wt"]
    ) * 100
    sensitivity.to_csv(data_dir / "figure4_cutoff_sensitivity_by_barcode.tsv", sep="\t", index=False)
    pooled.to_csv(data_dir / "figure4_cutoff_sensitivity_pooled.tsv", sep="\t", index=False)
    difference.to_csv(data_dir / "figure4_cutoff_difference.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(13.0, 5.1), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.05, 0.8])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])

    y_positions: list[float] = []
    y_labels: list[str] = []
    current_y = 0.0
    for barcode in BARCODE_ORDER:
        origins = ["Host-mapped"] if barcode == "bc04" else ["Host-mapped", "MHV-mapped"]
        for origin in origins:
            row = polarity[
                (polarity["mapping_barcode"] == barcode)
                & (polarity["read_origin"] == origin)
            ].iloc[0]
            y_positions.append(current_y)
            y_labels.append(f"{barcode} {origin.replace('-mapped', '')}")
            negative = row["negative_direction_fraction_callable"] * 100
            positive = row["positive_direction_fraction_callable"] * 100
            ax_a.barh(current_y, -negative, height=0.58, color="#4477AA", edgecolor="white")
            ax_a.barh(current_y, positive, height=0.58, color="#EE7733", edgecolor="white")
            ax_a.text(
                -negative - 0.35,
                current_y,
                f"{negative:.1f}",
                va="center",
                ha="right",
                fontsize=7.4,
            )
            ax_a.text(
                positive + 0.35,
                current_y,
                f"{positive:.1f}",
                va="center",
                ha="left",
                fontsize=7.4,
            )
            current_y += 0.82
        current_y += 0.48
    ax_a.axvline(0, color="#2B2B2B", linewidth=0.9)
    ax_a.set_yticks(y_positions, y_labels)
    ax_a.invert_yaxis()
    ax_a.set_xlim(-22, 22)
    ax_a.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{abs(value):.0f}"))
    ax_a.set_xlabel("Callable reads beyond 5-pA cutoff (%)")
    ax_a.set_title("Signed direction exposes barcode structure")
    ax_a.legend(
        handles=[
            Patch(facecolor="#4477AA", label="Negative direction"),
            Patch(facecolor="#EE7733", label="Positive direction"),
        ],
        loc="upper right",
    )
    polish_axis(ax_a, grid_axis="x")
    panel_label(ax_a, "A")

    for barcode in INFECTED_BARCODES:
        group = BARCODE_GROUP[barcode]
        subset = sensitivity[sensitivity["mapping_barcode"] == barcode]
        ax_b.plot(
            subset["cutoff_pa"],
            subset["terminal_junction_nonA_like_fraction_callable"] * 100,
            color=GROUP_COLOR[group],
            linewidth=1.0,
            alpha=0.42,
            marker="o",
            markersize=3.5,
        )
    for group in VIRUS_GROUP_ORDER:
        subset = pooled[pooled["sample_group"] == group]
        ax_b.plot(
            subset["cutoff_pa"],
            subset["terminal_junction_nonA_like_fraction_callable"] * 100,
            color=GROUP_COLOR[group],
            linewidth=2.4,
            marker="D",
            markersize=5.5,
            label=f"{GROUP_LABEL[group]} pooled",
        )
    ax_b.axvline(5, color="#333333", linestyle="--", linewidth=1.0)
    ax_b.set_xticks([3, 5, 7])
    ax_b.set_ylim(0, 50)
    ax_b.set_xlabel("Fixed absolute cutoff (pA)")
    ax_b.set_ylabel("MHV terminal-junction non-A-like (%)")
    ax_b.set_title("Absolute rates depend on the fixed cutoff")
    ax_b.legend(loc="upper right")
    polish_axis(ax_b)
    panel_label(ax_b, "B")

    ax_c.plot(
        difference["cutoff_pa"],
        difference["mut_minus_wt_pp"],
        color=GROUP_COLOR["MHV-mut"],
        linewidth=2.2,
        marker="o",
        markersize=6,
    )
    ax_c.fill_between(
        difference["cutoff_pa"],
        0,
        difference["mut_minus_wt_pp"],
        color=GROUP_COLOR["MHV-mut"],
        alpha=0.12,
    )
    for _, row in difference.iterrows():
        ax_c.text(
            row["cutoff_pa"],
            row["mut_minus_wt_pp"] + 0.28,
            f"+{row['mut_minus_wt_pp']:.2f}",
            ha="center",
            va="bottom",
            fontsize=8.4,
        )
    ax_c.axhline(0, color="#333333", linewidth=0.8)
    ax_c.axvline(5, color="#333333", linestyle="--", linewidth=1.0)
    ax_c.set_xticks([3, 5, 7])
    ax_c.set_ylim(-0.2, 6.0)
    ax_c.set_xlabel("Fixed absolute cutoff (pA)")
    ax_c.set_ylabel("Pooled Mut − WT (percentage points)")
    ax_c.set_title("The pooled contrast is cutoff-dependent")
    polish_axis(ax_c)
    panel_label(ax_c, "C")

    fig.suptitle(
        "Figure 4. Directional sample effects and fixed-cutoff sensitivity",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.028,
        (
            "The primary analysis is the two-sided 5-pA call. Polarity and 3/7-pA results are technical context, "
            "not alternative selected endpoints."
        ),
        ha="center",
        va="top",
        fontsize=8.2,
        color="#4A4A4A",
    )
    save_figure(fig, output_dir, "figure4_direction_and_cutoff_sensitivity", formats, dpi)


def main() -> None:
    args = parse_args()
    package_root = args.package_root.resolve()
    output_dir = package_root / "figures"
    data_dir = output_dir / "source_data"
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    allowed = {"png", "pdf", "svg"}
    unexpected = sorted(set(formats) - allowed)
    if unexpected:
        raise ValueError(f"Unsupported output format(s): {', '.join(unexpected)}")

    configure_style()
    figure_1(package_root, output_dir, data_dir, formats, args.dpi)
    figure_2(package_root, output_dir, data_dir, formats, args.dpi)
    figure_3(package_root, output_dir, data_dir, formats, args.dpi)
    figure_4(package_root, output_dir, data_dir, formats, args.dpi)
    print(f"Wrote publication figures and source data to {output_dir}")


if __name__ == "__main__":
    main()
