# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Full-plane kernel timing maps and paired summaries at equal accuracy."""
import csv
import html
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as effects
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

sys.dont_write_bytecode = True
from experiments.section51.run import GRID, LANGUAGES, METHODS, write_json
from marcum._paths import ROOT
OUT = ROOT / "data/section51"
FIGURES = ROOT / "results/figures/section51"
FIGURES.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(OUT.parent / "hybrid_v1_v2_quick"))
from bound import upper
def kappa_bound(a,b,nu):
    return upper(a,b,nu)[0]


LANGUAGE_NAMES = {"python": "Python", "cpp": "C++", "fortran": "Fortran"}
METHOD_NAMES = {"se": "SE", "af": "AF", "ni": "NI", "se_rr": "SE-RR"}
ORDERS = (1, 10, 20)
GRAY = "#d4d4d4"
ROW_MAX_US = {"python": 100., "cpp": 2., "fortran": 2.}
OVER_COLOR = GRAY


def table_html(rows):
    labels = ("nu", "language", "points", "SE", "AF", "NI", "SE-RR", "SE / SE-RR")
    output = ["<table><thead><tr>" + "".join("<th>%s</th>" % html.escape(x) for x in labels)
              + "</tr></thead><tbody>"]
    for row in rows:
        output.append("<tr>" + "".join("<td>%s</td>" % (
            ("%.3g" % row[x]) if isinstance(row[x], float) else html.escape(str(row[x])))
            for x in labels) + "</tr>")
    output.append("</tbody></table>")
    return "".join(output)


def main():
    status = json.loads((OUT / "status.json").read_text(encoding="utf-8"))
    assert status["complete"] and status["publication_timing_allowed"]
    rows = [json.loads(line) for line in (OUT / "timings.jsonl").read_text(encoding="utf-8").splitlines()]
    assert len(rows) == 57600
    index = {(r["nu"], r["language"], r["method"], r["b"], r["a"]): r for r in rows}
    assert len(index) == len(rows)
    arrays, coverage = {}, []
    for nu in ORDERS:
        for language in LANGUAGES:
            for method in METHODS:
                selected = [[index[(nu, language, method, b, a)] for a in GRID] for b in GRID]
                time = np.array([[r["microseconds"] if r["setting"] is not None else np.nan
                                  for r in line] for line in selected])
                assert np.all(time[np.isfinite(time)] > 0)
                arrays[(nu, language, method)] = time
                coverage.append(dict(nu=nu, language=language, method=method,
                                     qualified=int(np.isfinite(time).sum()), total=time.size))

    common_rows, core_rows = [], []
    for nu in ORDERS:
        all_common = np.logical_and.reduce([
            np.isfinite(arrays[(nu, language, method)])
            for language in LANGUAGES for method in METHODS])
        core_common = np.logical_and.reduce([
            np.isfinite(arrays[(nu, language, method)])
            for language in LANGUAGES for method in ("se", "ni", "se_rr")])
        assert all_common.any() and core_common.any()
        for mask, output in ((all_common, common_rows), (core_common, core_rows)):
            for language in LANGUAGES:
                row = dict(nu=nu, language=LANGUAGE_NAMES[language], points=int(mask.sum()))
                for method in METHODS:
                    time = arrays[(nu, language, method)]
                    row[METHOD_NAMES[method]] = float(time[mask].mean()) if np.isfinite(time[mask]).all() else "—"
                row["SE / SE-RR"] = row["SE"] / row["SE-RR"]
                row["AF / SE-RR"] = row["AF"] / row["SE-RR"] if isinstance(row["AF"],float) else None
                row["NI / SE-RR"] = row["NI"] / row["SE-RR"]
                output.append(row)
    for filename, data in (("summary_all_four.csv", common_rows), ("summary_core.csv", core_rows),
                           ("coverage.csv", coverage)):
        with (FIGURES / filename).open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    with (FIGURES / "point_results.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    latex = [r"\begin{table}[htbp]", r"\centering", r"\begin{tabular}{c l r r r r r r r}",
             r"\toprule", r"$\nu$ & Language & SE & AF & NI & SE-RR & vs. SE & vs. AF & vs. NI \\", r"\midrule"]
    for row in common_rows:
        times = [row[name] for name in ("SE", "AF", "NI", "SE-RR")]
        fastest = min(times)
        formatted = [r"\textbf{%.3g}" % value if value == fastest else "%.3g" % value for value in times]
        ratios = [row[key] for key in ("SE / SE-RR", "AF / SE-RR", "NI / SE-RR")]
        latex.append("%d & %s & %s & %s \\\\" % (row["nu"], row["language"], " & ".join(formatted), " & ".join("%.2f" % value for value in ratios)))
    latex.extend([r"\bottomrule", r"\end{tabular}",
                  r"\caption{Mean computation times in microseconds and SE-RR speedups on the common grid points.}",
                  r"\label{tab:individual-method-timing}", r"\end{table}"])
    (FIGURES / "summary_table.tex").write_text("\n".join(latex)+"\n", encoding="utf-8")

    # Each language has one linear range shared by its methods and all orders.
    maximum = {language: max(float(np.nanmax(arrays[(nu, language, method)]))
                            for nu in ORDERS for method in METHODS)
               for language in LANGUAGES}
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.linewidth": .65,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })
    width_mm = 183
    for nu in ORDERS:
        fig, axes = plt.subplots(3, 4, figsize=(width_mm / 25.4, 151 / 25.4),
                                 sharex=True, sharey=True)
        fig.subplots_adjust(left=.105, right=.86, top=.91, bottom=.115, wspace=.12, hspace=.18)
        for row, language in enumerate(LANGUAGES):
            cmap = plt.get_cmap("viridis_r")
            cmap.set_bad(GRAY)
            cmap.set_over(OVER_COLOR)
            norm = Normalize(0, ROW_MAX_US[language], clip=False)
            for col, method in enumerate(METHODS):
                ax = axes[row, col]
                data = arrays[(nu, language, method)]
                mesh = ax.imshow(np.ma.masked_invalid(data), origin="lower",
                                 extent=(.25, 20.25, .25, 20.25), interpolation="nearest",
                                 cmap=cmap, norm=norm, aspect="equal")
                if row == 0:
                    ax.set_title(METHOD_NAMES[method], pad=7)
                if col == 0:
                    ax.set_ylabel(r"$b$", labelpad=2)
                if row == 2:
                    ax.set_xlabel(r"$a$", labelpad=2)
                ax.set_xticks([.5, 10, 20])
                ax.set_yticks([.5, 10, 20])
                ax.set_xticklabels(["0.5", "10", "20"])
                ax.set_yticklabels(["0.5", "10", "20"])
                ax.tick_params(length=2.5, pad=2)
                if method == "se_rr":
                    fine = np.linspace(.25, 20.25, 401)
                    a, b = np.meshgrid(fine, fine)
                    region = (b <= 7) & (kappa_bound(a, b, np.full_like(a, nu)) <= 1000)
                    lines = ax.contour(fine, fine, region, levels=[.5], colors=["#151515"],
                                       linewidths=[.8], linestyles=["--"])
                    for line in lines.collections:
                        line.set_path_effects([effects.withStroke(linewidth=1.8, foreground="white")])
                ax.set_xlim(.25, 20.25)
                ax.set_ylim(.25, 20.25)
            position = axes[row, -1].get_position()
            cax = fig.add_axes([.89, position.y0, .018, position.height])
            ticks = np.linspace(0, ROW_MAX_US[language], 5)
            bar = fig.colorbar(mesh, cax=cax, extend="max", ticks=ticks)
            bar.ax.set_yticklabels(["%g" % value for value in ticks])
            bar.set_label("Time (µs)", fontsize=8, labelpad=4)
            bar.ax.tick_params(labelsize=7.5, length=2.5)
            ypos = (position.y0 + position.y1) / 2
            fig.text(.021, ypos, LANGUAGE_NAMES[language], ha="center", va="center",
                     rotation=90, rotation_mode="anchor", fontsize=9)
        fig.legend(handles=[
            Patch(facecolor=GRAY, edgecolor="none", label="Accuracy not attained or time above row scale"),
            Line2D([0], [0], color="#151515", linestyle="--", linewidth=.8,
                   label=r"Boundary of $b\leq7,\ B_{a,b,\nu}\leq10^3$"),
        ], loc="lower center", bbox_to_anchor=(.5, .016), ncol=2,
            frameon=False, fontsize=8, columnspacing=1.5)
        fig.savefig(str(FIGURES / ("heatmap_language_scales_nu%d.pdf" % nu)), dpi=450)
        fig.savefig(str(FIGURES / ("heatmap_language_scales_nu%d.svg" % nu)), dpi=450)
        fig.savefig(str(FIGURES / ("heatmap_language_scales_nu%d.png" % nu)), dpi=600)
        plt.close(fig)

    np.savez_compressed(FIGURES / "heatmap_data.npz", a=GRID, b=GRID,
                        **{"nu%d_%s_%s" % key: value for key, value in arrays.items()})
    write_json(FIGURES / "summary.json", dict(all_four_common=common_rows, core_common=core_rows,
               coverage=coverage, linear_maximum_microseconds=maximum,
               display_ranges_microseconds={lang: [0, top] for lang, top in ROW_MAX_US.items()},
               scale_type="Linear within each language; each language retains its range across all orders.",
               over_range_color=OVER_COLOR,
               common_sets="Identical point set for all three languages at each fixed order.",
               gray="Either no allowed setting attains relative error 1e-10 or qualified time exceeds that language's display range. Slow completed cells remain qualified in all raw data and summaries.",
               search_limits={"SE": [0, 1000], "AF": [1, 100], "NI": [0, 10], "SE-RR": [1, 1000]}))
    figures = "".join('<h2>ν = %d</h2><img src="heatmap_language_scales_nu%d.png" alt="Four-method timing maps">' %
                      (nu, nu) for nu in ORDERS)
    report = """<!doctype html><html lang="en"><meta charset="utf-8">
<title>Full-plane timing at equal accuracy</title><style>
body{font:16px/1.6 system-ui,sans-serif;max-width:1180px;margin:40px auto;color:#20242b;padding:0 24px}
img{width:100%;height:auto}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;margin:24px 0}
th,td{padding:10px 12px;border-bottom:1px solid #ddd;text-align:right}th{background:#eef2f4}
h1{font-size:26px}h2{margin-top:36px;font-size:21px}p{max-width:1000px}
</style><h1>Full-plane timing of four methods at equal accuracy</h1>
<p>a,b = 0.5,1.0,...,20, with 1,600 points per order. The target relative error is 1e-10.
Independent Bessel-integral references and accuracy calibration are computed outside the timers.
Each cell reports the median time over five batches of repeated evaluations.
The linear color scales are 0-100 microseconds for Python and 0-2 microseconds for C++ and Fortran.
All four methods in each row share one scale, and each language uses the same range across the three orders.
Gray indicates failure to attain the target within the tested search range, or a qualified computation time above the row's display range.
Completed computations above the display limit remain in the raw data and summary statistics.
Dashed lines in the SE-RR panels bound the region where b &lt;= 7 and the condition-number bound is at most 1000.</p>
<p>The external search ranges are 0-1000 for SE, 1-100 for AF, refinement levels 0-10 for NI,
and 1-1000 for SE-RR. Integer settings are tested in ascending order, and each method is rerun with the first qualifying setting for timing.
Cells that fail to qualify within this search range are gray.</p>
<h2>Mean time on common qualifying points (microseconds)</h2>
<p>All three languages and all four methods use the same common qualifying point set at each order.</p>""" + table_html(common_rows)
    report += "<h2>Mean time on common qualifying points for SE, NI and SE-RR (microseconds)</h2>"
    report += "<p>This table uses the larger point set where these three methods qualify in all languages. A dash is shown when AF does not qualify over the entire set.</p>"
    report += table_html(core_rows) + figures + "</html>"
    (FIGURES / "report.html").write_text(report, encoding="utf-8")
    print(json.dumps(dict(all_four=common_rows, core=core_rows, maxima=maximum), indent=2))


if __name__ == "__main__":
    main()
