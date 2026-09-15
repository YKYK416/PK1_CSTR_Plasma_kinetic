"""Create a source-locked pathway-concentration atlas for direct NH formation.

The atlas is a descriptive information-theoretic re-expression of the six
integrated direct-NH family shares already used in Figures 3 and 13.  It does
not introduce a new kinetic model, uncertainty distribution, or independent
validation layer.  Its purpose is to quantify pathway concentration and to
identify the second-ranked direct-NH family across the accepted CSTR map.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "analysis" / "figures_20260829" / "Figure_13_direct_NH_dominance_margin_source_data.csv"
OUT = ROOT / "analysis" / "p17_direct_nh_pathway_concentration_20260901_r1"

FAMILY_COLUMNS = [
    "named_h2star_pct",
    "rydberg_h2star_pct",
    "n2d_h2_pct",
    "n2p_h2_pct",
    "vib_h2_pct",
    "association_pct",
]
SECONDARY_CODE = {
    "n2p_h2_pct": 0,
    "rydberg_h2star_pct": 1,
    "association_pct": 2,
}
SECONDARY_LABEL = {
    "n2p_h2_pct": r"N(2P) + H$_2$",
    "rydberg_h2star_pct": r"Rydberg-H$_2^*$",
    "association_pct": "association",
}
SECONDARY_COLOR = ["#0072B2", "#D55E00", "#009E73"]  # Okabe-Ito derived.


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def grid(frame: pd.DataFrame, column: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    en = np.sort(frame.en_td.unique())
    x_n2 = np.sort(frame.n2_fraction.unique())
    values = frame.pivot(index="en_td", columns="n2_fraction", values=column).loc[en, x_n2]
    return en, x_n2, values.to_numpy(dtype=float)


def heatmap(ax, values, en, x_n2, title, cmap, vmin, vmax, cbar_label, fmt=".2f"):
    image = ax.imshow(values, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(len(x_n2)), [f"{value:.1f}" for value in x_n2])
    ax.set_yticks(np.arange(len(en)), [str(int(value)) for value in en])
    ax.set_xlabel(r"$x_{\mathrm{N_2}}$ in feed")
    ax.set_ylabel(r"$E/N$ (Td)")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=9)
    cb = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(cbar_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    for spine in ax.spines.values():
        spine.set_linewidth(0.65)
    return image


def main() -> None:
    require(SOURCE.exists(), f"Missing source data: {SOURCE}")
    OUT.mkdir(parents=True, exist_ok=True)
    source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    frame = pd.read_csv(SOURCE)
    require(len(frame) == 108, f"Expected 108 accepted states, found {len(frame)}")
    require(frame.run_status.isin(["success", "retry_success"]).all(), "Only accepted map states are permitted")

    percentages = frame[FAMILY_COLUMNS].copy()
    total = percentages.sum(axis=1)
    require(np.allclose(total, 100.0, atol=1e-7), "Direct-NH family shares do not close to 100%")
    probabilities = percentages.div(total, axis=0)
    entropy = -(probabilities.where(probabilities > 0, 1.0) * np.log(probabilities.where(probabilities > 0, 1.0))).sum(axis=1)
    entropy_norm = entropy / np.log(len(FAMILY_COLUMNS))

    frame["pathway_entropy_nats"] = entropy
    frame["pathway_entropy_normalized"] = entropy_norm
    frame["effective_pathway_count"] = np.exp(entropy)
    frame["runner_up_code"] = frame.largest_secondary_column.map(SECONDARY_CODE)
    require(frame.runner_up_code.notna().all(), "Unexpected runner-up family")
    frame["runner_up_label"] = frame.largest_secondary_column.map(SECONDARY_LABEL)
    frame.to_csv(OUT / "Figure_26_direct_nh_pathway_concentration_source_data.csv", index=False)

    en, x_n2, h_grid = grid(frame, "pathway_entropy_normalized")
    _, _, n_eff_grid = grid(frame, "effective_pathway_count")
    _, _, runner_grid = grid(frame, "runner_up_code")

    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(7.15, 6.45), constrained_layout=False)
    grid_spec = fig.add_gridspec(
        2, 2, height_ratios=[1.0, 1.04],
        left=0.085, right=0.975, top=0.855, bottom=0.155, hspace=0.52, wspace=0.52,
    )
    ax_a = fig.add_subplot(grid_spec[0, 0])
    ax_b = fig.add_subplot(grid_spec[0, 1])
    ax_c = fig.add_subplot(grid_spec[1, 0])
    ax_d = fig.add_subplot(grid_spec[1, 1])

    heatmap(
        ax_a, h_grid, en, x_n2,
        "A  Normalized direct-NH pathway entropy",
        "cividis", 0.0, max(0.35, float(h_grid.max())),
        r"$H/\ln(6)$",
    )
    heatmap(
        ax_b, n_eff_grid, en, x_n2,
        "B  Effective number of direct-NH families",
        "viridis", 1.0, 2.0,
        r"$N_{\mathrm{eff}}=\exp(H)$",
    )

    cmap = ListedColormap(SECONDARY_COLOR)
    runner = ax_c.imshow(runner_grid, origin="lower", aspect="auto", cmap=cmap, vmin=-0.5, vmax=2.5)
    ax_c.set_xticks(np.arange(len(x_n2)), [f"{value:.1f}" for value in x_n2])
    ax_c.set_yticks(np.arange(len(en)), [str(int(value)) for value in en])
    ax_c.set_xlabel(r"$x_{\mathrm{N_2}}$ in feed")
    ax_c.set_ylabel(r"$E/N$ (Td)")
    ax_c.set_title("C  Identity of the second-ranked direct-NH family", loc="left", fontweight="bold", fontsize=9)
    code_label = {0: "P", 1: "R", 2: "A"}
    for row in range(runner_grid.shape[0]):
        for col in range(runner_grid.shape[1]):
            ax_c.text(col, row, code_label[int(runner_grid[row, col])], ha="center", va="center",
                      color="white", fontsize=6.3, fontweight="bold")
    ax_c.legend(
        handles=[Line2D([0], [0], marker="s", color="none", markerfacecolor=color, markeredgecolor="none",
                        markersize=7, label=f"{code} = {SECONDARY_LABEL[name]}")
                 for name, code, color in zip(SECONDARY_CODE, ["P", "R", "A"], SECONDARY_COLOR)],
        loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=1, frameon=False, fontsize=6.7,
        handletextpad=0.4,
    )
    for spine in ax_c.spines.values():
        spine.set_linewidth(0.65)

    point_colors = np.array([SECONDARY_COLOR[int(code)] for code in frame.runner_up_code])
    ax_d.scatter(frame.pathway_entropy_normalized, frame.log10_named_to_secondary,
                 c=point_colors, s=25, edgecolors="white", linewidths=0.35, alpha=0.95)
    ax_d.axhline(0.0, color="0.3", ls="--", lw=0.8, label="crossover threshold")
    ax_d.set_xlabel(r"Normalized pathway entropy, $H/\ln(6)$")
    ax_d.set_ylabel(r"$\log_{10}$(named-H$_2^*$/runner-up)")
    ax_d.set_title("D  Concentration and direct-NH dominance", loc="left", fontweight="bold", fontsize=9)
    ax_d.set_xlim(-0.01, max(0.35, float(frame.pathway_entropy_normalized.max()) + 0.02))
    ax_d.set_ylim(-0.15, max(2.5, float(frame.log10_named_to_secondary.max()) + 0.1))
    ax_d.grid(axis="y", color="0.9", lw=0.6)
    ax_d.spines[["top", "right"]].set_visible(False)
    ax_d.text(
        0.98, 0.06,
        f"108 accepted states\nmax H/ln(6) = {frame.pathway_entropy_normalized.max():.3f}\n"
        f"max N_eff = {frame.effective_pathway_count.max():.3f}\n"
        f"min ratio = {frame.named_to_secondary_ratio.min():.3f}",
        transform=ax_d.transAxes, ha="right", va="bottom", fontsize=7,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "0.7", "alpha": 0.95},
    )

    fig.suptitle("Direct-NH pathway concentration", x=0.5, y=0.972,
                 fontsize=10.4, fontweight="bold")
    fig.text(0.5, 0.936,
             "Descriptive transform of the Figure 3/13 terminal integrated-ROP shares; not independent kinetic evidence.",
             ha="center", va="center", fontsize=7.1, color="0.28")
    fig.savefig(OUT / "Figure_26_direct_nh_pathway_concentration.png", dpi=400, bbox_inches="tight")
    fig.savefig(OUT / "Figure_26_direct_nh_pathway_concentration.pdf", bbox_inches="tight")
    plt.close(fig)

    max_row = frame.loc[frame.pathway_entropy_normalized.idxmax()]
    report = f"""# P17 direct-NH pathway-concentration audit

**Scope:** An information-theoretic re-expression of the six direct-NH family shares already reported in the accepted continuous-wave CSTR map. It is descriptive and does not introduce a new mechanism, probability distribution, or independent validation dataset.

## Source lock

- Source CSV: `{SOURCE.relative_to(ROOT)}`
- SHA-256: `{source_hash}`
- Accepted states: {len(frame)}
- Family-share closure: {total.min():.12f}–{total.max():.12f}%

## Definitions

For direct-NH family shares \(p_i\), \(H=-\\sum_i p_i\\ln p_i\), normalized entropy is \(H/\\ln(6)\), and the effective number of active families is \(N_{{\\mathrm{{eff}}}}=\\exp(H)\). A larger value indicates a less concentrated distribution of the already integrated family fluxes; it is not a kinetic rate, route fraction beyond the supplied ROP grouping, or uncertainty metric.

## Results

- Normalized entropy range: {frame.pathway_entropy_normalized.min():.6f}–{frame.pathway_entropy_normalized.max():.6f}.
- Effective pathway count range: {frame.effective_pathway_count.min():.6f}–{frame.effective_pathway_count.max():.6f}.
- Lowest named-H2*/runner-up ratio: {frame.named_to_secondary_ratio.min():.6f}.
- Maximum entropy occurs at {int(max_row.en_td)} Td and \\(x_{{\\mathrm{{N_2}}}}={max_row.n2_fraction:.1f}\\): named-H2* share = {max_row.named_h2star_pct:.6f}%, runner-up = {max_row.runner_up_label}, ratio = {max_row.named_to_secondary_ratio:.6f}.
- Runner-up identity counts: {frame.runner_up_label.value_counts().to_dict()}.

## Interpretation boundary

The atlas quantifies that the direct-NH distribution remains concentrated even at its most diverse accepted state. It reinforces the no-crossover observation by a complementary descriptive metric, but it must not be described as independent evidence, a formal uncertainty analysis, a surface-pathway result, or a universal plasma-ammonia conclusion.
"""
    (OUT / "P17_direct_nh_pathway_concentration_audit.md").write_text(report, encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
