#!/usr/bin/env python3
"""Build the overall WT-versus-Mut MHV terminal-junction figure and source data.

The figure adapts the full-composition plus magnified-tail layout of the
provided Figure 6A. The source-data export contains all MHV-mapped reads from
infected libraries (bc05/bc06, WT; bc07/bc11, Mut), including uncallable rows
and their QC fields. Plotting and endpoint summaries use only callable rows.
Mock-library MHV mappings are excluded.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from statistics import NormalDist

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import ConnectionPatch, Patch, Rectangle


INFECTED_BARCODES = ["bc05", "bc06", "bc07", "bc11"]
GROUP_ORDER = ["MHV-wt", "MHV-mut"]
DISPLAY_GROUP = {"MHV-wt": "WT", "MHV-mut": "Mut"}
CLASS_ORDER = [
    "A_junction_like",
    "terminal_junction_nonA_like",
]
CLASS_LABEL = {
    "A_junction_like": "A-junction-like",
    "terminal_junction_nonA_like": "Non-A-like",
}
CLASS_COLOR = {
    "A_junction_like": "#D9D9D9",
    "terminal_junction_nonA_like": "#2B6CB0",
}


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    publication_root = script_path.parents[2]
    output_root = script_path.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-parquet",
        type=Path,
        default=(
            publication_root
            / "results"
            / "final_primary"
            / "mhv_reads_terminal_nona_direct.parquet"
        ),
        help="Primary per-read MHV result table.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=output_root,
        help="Directory containing data/ and output/ subdirectories.",
    )
    parser.add_argument(
        "--figure-mode",
        choices=["all", "pooled", "barcode"],
        default="all",
        help="Which figure variant to write (default: all).",
    )
    parser.add_argument("--dpi", type=int, default=400, help="PNG resolution.")
    return parser.parse_args()


def wilson_interval(successes: int, total: int, confidence: float = 0.95) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if total <= 0:
        return (float("nan"), float("nan"))
    z = NormalDist().inv_cdf(0.5 + confidence / 2)
    proportion = successes / total
    denominator = 1 + z**2 / total
    center = (proportion + z**2 / (2 * total)) / denominator
    half_width = (
        z
        / denominator
        * math.sqrt(
            proportion * (1 - proportion) / total + z**2 / (4 * total**2)
        )
    )
    return center - half_width, center + half_width


def load_and_filter(source_parquet: Path) -> pd.DataFrame:
    columns = [
        "read_id",
        "mapping_barcode",
        "sample_group",
        "mapping_group",
        "barcode_mapping_consistent",
        "boundary_method",
        "eligible_terminal_nona_direct",
        "terminal_nona_direct_qc_status",
        "terminal_window_mean_centered_pa",
        "stableA_window_mean_centered_pa",
        "terminal_minus_stableA_pa",
        "pooled_long_host_A_DNA_reference_pa",
        "terminal_nona_fixed_cutoff_pa",
        "terminal_nona_signed_difference_pa",
        "terminal_nona_abs_difference_pa",
        "terminal_junction_nonA_like",
        "terminal_nona_direct_class",
    ]
    reads = pd.read_parquet(source_parquet, columns=columns)
    reads = reads[
        reads["mapping_barcode"].isin(INFECTED_BARCODES)
        & reads["sample_group"].isin(GROUP_ORDER)
        & reads["mapping_group"].eq("virus")
    ].copy()

    if reads["read_id"].duplicated().any():
        raise ValueError("Filtered infected-virus table contains duplicate read IDs")
    if len(reads) != 84_900:
        raise ValueError(f"Unexpected infected-virus row count: {len(reads):,}")

    callable_mask = reads["eligible_terminal_nona_direct"].eq(1)
    callable_reads = reads.loc[callable_mask]
    if len(callable_reads) != 69_799:
        raise ValueError(f"Unexpected callable infected-virus count: {len(callable_reads):,}")
    if callable_reads["terminal_nona_signed_difference_pa"].isna().any():
        raise ValueError("Callable reads contain missing signed deviations")

    cutoff_values = callable_reads["terminal_nona_fixed_cutoff_pa"].dropna().unique()
    reference_values = callable_reads[
        "pooled_long_host_A_DNA_reference_pa"
    ].dropna().unique()
    if len(cutoff_values) != 1 or not math.isclose(float(cutoff_values[0]), 5.0):
        raise ValueError(f"Unexpected cutoff values: {cutoff_values}")
    if len(reference_values) != 1:
        raise ValueError("Expected one pooled host reference value")

    signed = reads["terminal_nona_signed_difference_pa"].astype(float)
    cutoff = float(cutoff_values[0])
    reads["figure_signal_class"] = "uncallable"
    reads.loc[callable_mask, "figure_signal_class"] = np.where(
        signed.loc[callable_mask].abs() >= cutoff,
        "terminal_junction_nonA_like",
        "A_junction_like",
    )
    derived_positive = reads.loc[callable_mask, "figure_signal_class"].eq(
        "terminal_junction_nonA_like"
    )
    stored_positive = reads.loc[
        callable_mask, "terminal_junction_nonA_like"
    ].astype(bool)
    if not np.array_equal(derived_positive.to_numpy(), stored_positive.to_numpy()):
        raise ValueError("Derived tail classes do not reproduce the primary endpoint")

    reads["display_group"] = reads["sample_group"].map(DISPLAY_GROUP)
    return reads


def summarize(reads: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_key: str | list[str] = group_columns[0] if len(group_columns) == 1 else group_columns
    for keys, subset in reads.groupby(group_key, observed=True, sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_columns, keys))
        n_total_mhv = len(subset)
        callable_subset = subset[subset["eligible_terminal_nona_direct"].eq(1)]
        counts = callable_subset["figure_signal_class"].value_counts()
        total = len(callable_subset)
        if total == 0:
            raise ValueError(f"No callable reads for group {keys}")
        n_a_like = int(counts.get("A_junction_like", 0))
        n_non_a = int(counts.get("terminal_junction_nonA_like", 0))
        ci_low, ci_high = wilson_interval(n_non_a, total)
        row.update(
            {
                "n_total_mhv": n_total_mhv,
                "n_callable": total,
                "n_uncallable": n_total_mhv - total,
                "callable_fraction": total / n_total_mhv,
                "n_A_junction_like": n_a_like,
                "n_terminal_junction_nonA_like": n_non_a,
                "A_junction_like_fraction": n_a_like / total,
                "terminal_junction_nonA_like_fraction": n_non_a / total,
                "wilson95_low_descriptive": ci_low,
                "wilson95_high_descriptive": ci_high,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def effect_summary(group_summary: pd.DataFrame) -> pd.DataFrame:
    indexed = group_summary.set_index("sample_group")
    wt = indexed.loc["MHV-wt"]
    mut = indexed.loc["MHV-mut"]
    wt_positive = int(wt["n_terminal_junction_nonA_like"])
    wt_total = int(wt["n_callable"])
    mut_positive = int(mut["n_terminal_junction_nonA_like"])
    mut_total = int(mut["n_callable"])
    wt_fraction = wt_positive / wt_total
    mut_fraction = mut_positive / mut_total
    difference = mut_fraction - wt_fraction
    risk_ratio = mut_fraction / wt_fraction
    wt_negative = wt_total - wt_positive
    mut_negative = mut_total - mut_positive
    odds_ratio = (mut_positive * wt_negative) / (mut_negative * wt_positive)
    z = NormalDist().inv_cdf(0.975)
    rr_se = math.sqrt(
        1 / mut_positive
        - 1 / mut_total
        + 1 / wt_positive
        - 1 / wt_total
    )
    or_se = math.sqrt(
        1 / mut_positive + 1 / mut_negative + 1 / wt_positive + 1 / wt_negative
    )
    wt_ci_low, wt_ci_high = wilson_interval(wt_positive, wt_total)
    mut_ci_low, mut_ci_high = wilson_interval(mut_positive, mut_total)
    # Newcombe's interval for the difference between independent proportions,
    # constructed from the two Wilson score intervals (without continuity correction).
    difference_ci_low = difference - math.sqrt(
        (mut_fraction - mut_ci_low) ** 2 + (wt_ci_high - wt_fraction) ** 2
    )
    difference_ci_high = difference + math.sqrt(
        (mut_ci_high - mut_fraction) ** 2 + (wt_fraction - wt_ci_low) ** 2
    )
    return pd.DataFrame(
        [
            {
                "wt_positive": wt_positive,
                "wt_callable": wt_total,
                "wt_fraction": wt_fraction,
                "mut_positive": mut_positive,
                "mut_callable": mut_total,
                "mut_fraction": mut_fraction,
                "mut_minus_wt_fraction": difference,
                "mut_minus_wt_percentage_points": difference * 100,
                "mut_minus_wt_percentage_points_95ci_low_read_level": difference_ci_low
                * 100,
                "mut_minus_wt_percentage_points_95ci_high_read_level": difference_ci_high
                * 100,
                "relative_increase_fraction": risk_ratio - 1,
                "relative_increase_percent": (risk_ratio - 1) * 100,
                "read_pooled_risk_ratio": risk_ratio,
                "read_pooled_risk_ratio_95ci_low": math.exp(math.log(risk_ratio) - z * rr_se),
                "read_pooled_risk_ratio_95ci_high": math.exp(math.log(risk_ratio) + z * rr_se),
                "read_pooled_odds_ratio": odds_ratio,
                "read_pooled_odds_ratio_95ci_low": math.exp(math.log(odds_ratio) - z * or_se),
                "read_pooled_odds_ratio_95ci_high": math.exp(math.log(odds_ratio) + z * or_se),
            }
        ]
    )


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "axes.linewidth": 0.9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8.2,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def draw_figure(group_summary: pd.DataFrame, output_dir: Path, dpi: int) -> None:
    summary = group_summary.set_index("sample_group").reindex(GROUP_ORDER)
    x = np.arange(len(GROUP_ORDER), dtype=float)
    fractions = {
        signal_class: np.array(
            [
                summary.loc[group, f"{signal_class}_fraction"] * 100
                for group in GROUP_ORDER
            ]
        )
        for signal_class in CLASS_ORDER
    }
    totals = fractions["terminal_junction_nonA_like"]
    denominators = summary["n_callable"].astype(int).to_numpy()

    configure_style()
    fig = plt.figure(figsize=(8.9, 4.55))
    grid = fig.add_gridspec(1, 3, width_ratios=[0.95, 1.0, 1.12], wspace=0.48)
    ax_main = fig.add_subplot(grid[0, 0])
    ax_zoom = fig.add_subplot(grid[0, 1])
    ax_key = fig.add_subplot(grid[0, 2])
    width = 0.62

    bottom = np.zeros(len(x))
    for signal_class in CLASS_ORDER:
        values = fractions[signal_class]
        ax_main.bar(
            x,
            values,
            width=width,
            bottom=bottom,
            color=CLASS_COLOR[signal_class],
            edgecolor="white",
            linewidth=0.45,
        )
        bottom += values

    ax_main.set_ylim(0, 102)
    ax_main.set_xlim(-0.6, 1.6)
    ax_main.set_yticks(np.arange(0, 101, 20))
    ax_main.set_xticks(
        x,
        [f"WT\nn={denominators[0]:,}", f"Mut\nn={denominators[1]:,}"],
    )
    ax_main.set_ylabel("Fraction of callable MHV reads (%)")
    ax_main.spines[["top", "right"]].set_visible(False)

    crop_bottom = 80.0
    crop = Rectangle(
        (-0.48, crop_bottom),
        1.96,
        20.0,
        fill=False,
        edgecolor="#F04B43",
        linewidth=1.0,
        linestyle=(0, (2.2, 2.2)),
        clip_on=False,
    )
    ax_main.add_patch(crop)

    zoom_bottom = np.zeros(len(x))
    for signal_class in ["terminal_junction_nonA_like"]:
        values = fractions[signal_class]
        ax_zoom.bar(
            x,
            values,
            width=width,
            bottom=zoom_bottom,
            color=CLASS_COLOR[signal_class],
            edgecolor="white",
            linewidth=0.45,
        )
        zoom_bottom += values

    ax_zoom.set_ylim(0, 20.5)
    ax_zoom.set_xlim(-0.6, 1.6)
    ax_zoom.set_yticks([0, 4, 8, 12, 16, 20])
    ax_zoom.set_xticks(
        x,
        [f"WT\nn={denominators[0]:,}", f"Mut\nn={denominators[1]:,}"],
    )
    ax_zoom.set_ylabel("Non-A-like reads (%)")
    ax_zoom.set_title("Expanded Non-A-like fraction", pad=8, fontweight="bold")
    ax_zoom.spines[["top", "right"]].set_visible(False)

    for xpos, total, group in zip(x, totals, GROUP_ORDER):
        positive = int(summary.loc[group, "n_terminal_junction_nonA_like"])
        ax_zoom.text(
            xpos,
            total + 0.35,
            f"{total:.2f}%\n({positive:,})",
            ha="center",
            va="bottom",
            fontsize=8.2,
            fontweight="bold",
        )

    difference = totals[1] - totals[0]
    ax_zoom.text(
        0.5,
        19.25,
        f"Mut - WT = +{difference:.2f} pp",
        ha="center",
        va="center",
        fontsize=8.3,
    )

    for y_main, y_zoom in [(crop_bottom, 0), (100, 20)]:
        connector = ConnectionPatch(
            xyA=(1.48, y_main),
            coordsA=ax_main.transData,
            xyB=(-0.48, y_zoom),
            coordsB=ax_zoom.transData,
            color="#F04B43",
            linewidth=0.9,
            linestyle=(0, (2.2, 2.2)),
            clip_on=False,
        )
        fig.add_artist(connector)

    handles = [
        Patch(facecolor=CLASS_COLOR[signal_class], edgecolor="none", label=CLASS_LABEL[signal_class])
        for signal_class in CLASS_ORDER
    ]
    ax_key.axis("off")
    ax_key.legend(
        handles=handles,
        title="Terminal-junction signal class",
        title_fontsize=9.5,
        loc="upper left",
        bbox_to_anchor=(0.0, 0.98),
        handlelength=1.2,
        handleheight=1.2,
        labelspacing=0.75,
    )
    ax_key.text(
        0.02,
        0.56,
        (
            "Overall Non-A-like fraction\n"
            f"WT:  {totals[0]:.2f}%\n"
            f"Mut: {totals[1]:.2f}%\n"
            f"Difference: +{difference:.2f} percentage points"
        ),
        transform=ax_key.transAxes,
        ha="left",
        va="top",
        fontsize=9.0,
        linespacing=1.45,
    )
    ax_key.text(
        0.02,
        0.28,
        (
            "Non-A-like = absolute current deviation\n"
            "from the pooled host reference >= 5 pA."
        ),
        transform=ax_key.transAxes,
        ha="left",
        va="top",
        fontsize=8.2,
        color="#4A4A4A",
        linespacing=1.35,
    )

    fig.text(0.018, 0.985, "A", fontsize=22, fontweight="bold", va="top", ha="left")
    fig.subplots_adjust(left=0.095, right=0.985, bottom=0.18, top=0.90)

    for extension in ["pdf", "svg", "png"]:
        kwargs = {"dpi": dpi} if extension == "png" else {}
        fig.savefig(
            output_dir / f"figure6A_overall_mut_vs_wt_virus.{extension}",
            **kwargs,
        )
    plt.close(fig)


def draw_barcode_figure(
    barcode_summary: pd.DataFrame, output_dir: Path, dpi: int
) -> None:
    """Draw the same composition view with each infected barcode separated."""
    barcode_order = ["bc05", "bc06", "bc07", "bc11"]
    summary = barcode_summary.set_index("mapping_barcode").reindex(barcode_order)
    x = np.array([0.0, 1.05, 2.65, 3.70])
    fractions = {
        signal_class: summary[f"{signal_class}_fraction"].to_numpy(dtype=float)
        * 100
        for signal_class in CLASS_ORDER
    }
    totals = fractions["terminal_junction_nonA_like"]
    tick_labels = barcode_order

    configure_style()
    fig = plt.figure(figsize=(10.2, 4.65))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.20, 1.22, 1.08], wspace=0.44)
    ax_main = fig.add_subplot(grid[0, 0])
    ax_zoom = fig.add_subplot(grid[0, 1])
    ax_key = fig.add_subplot(grid[0, 2])
    width = 0.62

    bottom = np.zeros(len(x))
    for signal_class in CLASS_ORDER:
        values = fractions[signal_class]
        ax_main.bar(
            x,
            values,
            width=width,
            bottom=bottom,
            color=CLASS_COLOR[signal_class],
            edgecolor="white",
            linewidth=0.45,
        )
        bottom += values

    ax_main.set_ylim(0, 102)
    ax_main.set_xlim(-0.60, 4.30)
    ax_main.set_yticks(np.arange(0, 101, 20))
    ax_main.set_xticks(x, tick_labels)
    ax_main.tick_params(axis="x", labelsize=8.5)
    ax_main.set_ylabel("Fraction of callable MHV reads (%)")
    ax_main.axvline(1.85, color="#B8B8B8", linewidth=0.8, linestyle=(0, (2, 2)))
    ax_main.spines[["top", "right"]].set_visible(False)
    ax_main.text(
        0.525,
        -0.115,
        "WT-associated",
        transform=ax_main.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=8.2,
        fontweight="bold",
    )
    ax_main.text(
        3.175,
        -0.115,
        "Mut-associated",
        transform=ax_main.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=8.2,
        fontweight="bold",
    )

    crop_bottom = 80.0
    crop_left = -0.46
    crop_right = 4.16
    crop = Rectangle(
        (crop_left, crop_bottom),
        crop_right - crop_left,
        20.0,
        fill=False,
        edgecolor="#F04B43",
        linewidth=1.0,
        linestyle=(0, (2.2, 2.2)),
        clip_on=False,
    )
    ax_main.add_patch(crop)

    ax_zoom.bar(
        x,
        totals,
        width=width,
        color=CLASS_COLOR["terminal_junction_nonA_like"],
        edgecolor="white",
        linewidth=0.45,
    )
    ax_zoom.set_ylim(0, 20.5)
    ax_zoom.set_xlim(-0.60, 4.30)
    ax_zoom.set_yticks([0, 4, 8, 12, 16, 20])
    ax_zoom.set_xticks(x, tick_labels)
    ax_zoom.tick_params(axis="x", labelsize=8.5)
    ax_zoom.set_ylabel("Non-A-like reads (%)")
    ax_zoom.set_title("Barcode-level Non-A-like fraction", pad=8, fontweight="bold")
    ax_zoom.axvline(1.85, color="#B8B8B8", linewidth=0.8, linestyle=(0, (2, 2)))
    ax_zoom.spines[["top", "right"]].set_visible(False)
    ax_zoom.text(
        0.525,
        -0.115,
        "WT-associated",
        transform=ax_zoom.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=8.2,
        fontweight="bold",
    )
    ax_zoom.text(
        3.175,
        -0.115,
        "Mut-associated",
        transform=ax_zoom.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=8.2,
        fontweight="bold",
    )

    for xpos, total in zip(x, totals):
        ax_zoom.text(
            xpos,
            total + 0.30,
            f"{total:.2f}%",
            ha="center",
            va="bottom",
            fontsize=7.6,
            fontweight="bold",
        )

    for y_main, y_zoom in [(crop_bottom, 0), (100, 20)]:
        connector = ConnectionPatch(
            xyA=(crop_right, y_main),
            coordsA=ax_main.transData,
            xyB=(crop_left, y_zoom),
            coordsB=ax_zoom.transData,
            color="#F04B43",
            linewidth=0.9,
            linestyle=(0, (2.2, 2.2)),
            clip_on=False,
        )
        fig.add_artist(connector)

    handles = [
        Patch(
            facecolor=CLASS_COLOR[signal_class],
            edgecolor="none",
            label=CLASS_LABEL[signal_class],
        )
        for signal_class in CLASS_ORDER
    ]
    ax_key.axis("off")
    ax_key.legend(
        handles=handles,
        title="Terminal-junction signal class",
        title_fontsize=9.5,
        loc="upper left",
        bbox_to_anchor=(0.0, 0.98),
        handlelength=1.2,
        handleheight=1.2,
        labelspacing=0.75,
    )
    ax_key.text(
        0.02,
        0.61,
        (
            "Barcode-level Non-A-like fraction\n"
            + "\n".join(
                f"{DISPLAY_GROUP[row['sample_group']]} / {barcode}: "
                f"{int(row['n_terminal_junction_nonA_like']):,} / "
                f"{int(row['n_callable']):,} "
                f"({100 * row['terminal_junction_nonA_like_fraction']:.2f}%)"
                for barcode, row in summary.iterrows()
            )
        ),
        transform=ax_key.transAxes,
        ha="left",
        va="top",
        fontsize=8.8,
        linespacing=1.38,
    )
    ax_key.text(
        0.02,
        0.25,
        (
            "Non-A-like = absolute current deviation\n"
            "from the pooled host reference >= 5 pA.\n"
            "Each bar uses its own callable-read denominator."
        ),
        transform=ax_key.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        color="#4A4A4A",
        linespacing=1.35,
    )

    fig.text(0.018, 0.985, "A", fontsize=22, fontweight="bold", va="top", ha="left")
    fig.subplots_adjust(left=0.083, right=0.988, bottom=0.225, top=0.90)

    for extension in ["pdf", "svg", "png"]:
        kwargs = {"dpi": dpi} if extension == "png" else {}
        fig.savefig(
            output_dir / f"figure6A_barcode_level_mut_vs_wt_virus.{extension}",
            **kwargs,
        )
    plt.close(fig)


def main() -> None:
    args = parse_args()
    source_parquet = args.source_parquet.resolve()
    output_root = args.output_root.resolve()
    data_dir = output_root / "data"
    figure_dir = output_root / "output"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    reads = load_and_filter(source_parquet)
    group_summary = summarize(reads, ["sample_group", "display_group"])
    barcode_summary = summarize(
        reads, ["sample_group", "display_group", "mapping_barcode"]
    )

    expected = {
        "MHV-wt": {
            "n_total_mhv": 23587,
            "n_callable": 19169,
            "n_terminal_junction_nonA_like": 2616,
        },
        "MHV-mut": {
            "n_total_mhv": 61313,
            "n_callable": 50630,
            "n_terminal_junction_nonA_like": 8163,
        },
    }
    for group, values in expected.items():
        row = group_summary[group_summary["sample_group"] == group].iloc[0]
        for column, expected_value in values.items():
            if int(row[column]) != expected_value:
                raise ValueError(
                    f"Unexpected {group} {column}: {row[column]} != {expected_value}"
                )

    raw_columns = [
        "read_id",
        "mapping_barcode",
        "sample_group",
        "display_group",
        "mapping_group",
        "barcode_mapping_consistent",
        "boundary_method",
        "eligible_terminal_nona_direct",
        "terminal_nona_direct_qc_status",
        "terminal_window_mean_centered_pa",
        "stableA_window_mean_centered_pa",
        "terminal_minus_stableA_pa",
        "pooled_long_host_A_DNA_reference_pa",
        "terminal_nona_fixed_cutoff_pa",
        "terminal_nona_signed_difference_pa",
        "terminal_nona_abs_difference_pa",
        "terminal_junction_nonA_like",
        "terminal_nona_direct_class",
        "figure_signal_class",
    ]
    reads[raw_columns].sort_values(
        ["sample_group", "mapping_barcode", "read_id"]
    ).to_csv(
        data_dir / "figure6A_overall_virus_read_level.tsv.gz",
        sep="\t",
        index=False,
        float_format="%.8f",
        compression="gzip",
    )
    barcode_summary.to_csv(
        data_dir / "figure6A_overall_virus_barcode_summary.tsv",
        sep="\t",
        index=False,
        float_format="%.12g",
    )
    if args.figure_mode in {"all", "pooled"}:
        group_summary.to_csv(
            data_dir / "figure6A_overall_virus_group_summary.tsv",
            sep="\t", index=False, float_format="%.12g",
        )
        effect_summary(group_summary).to_csv(
            data_dir / "figure6A_overall_virus_effect_summary.tsv",
            sep="\t", index=False, float_format="%.12g",
        )
        draw_figure(group_summary, figure_dir, args.dpi)
    if args.figure_mode in {"all", "barcode"}:
        draw_barcode_figure(barcode_summary, figure_dir, args.dpi)

    n_callable = int(reads["eligible_terminal_nona_direct"].eq(1).sum())
    print(f"Infected-virus rows: {len(reads):,}")
    print(f"Callable infected-virus rows: {n_callable:,}")
    print(group_summary.to_string(index=False))
    print(f"Wrote source data to {data_dir}")
    print(f"Wrote figure files to {figure_dir}")


if __name__ == "__main__":
    main()
