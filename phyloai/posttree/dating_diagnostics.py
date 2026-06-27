"""MCMCtree output parsing and diagnostic plot generation."""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr


def parse_mcmctree_out(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        m = re.match(
            r"^\s*(t_n\d+)\s+"
            r"([\d.]+)\s+\(\s*([\d.]+),\s*([\d.]+)\)",
            line,
        )
        if m:
            node, mean, lower, upper = m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))
            rows.append({
                "node": node,
                "mean": mean,
                "lower": lower,
                "upper": upper,
                "ci_width": upper - lower,
            })
    return rows


def extract_node_tree(text: str) -> str | None:
    marker = "Species tree for FigTree."
    idx = text.find(marker)
    if idx == -1:
        return None
    after = text[idx:]
    first = re.search(r"(\([\s\S]+?\)[\s\d]*;)", after)
    return first.group(1).strip() if first else None


def build_time_table(run1: list[dict], run2: list[dict],
                     label1: str = "run1", label2: str = "run2") -> list[dict]:
    by_node = {r["node"]: r for r in run2}
    result = []
    for r in run1:
        node = r["node"]
        r2 = by_node.get(node, {})
        result.append({
            "node": node,
            f"mean_{label1}": r["mean"],
            f"lower_{label1}": r["lower"],
            f"upper_{label1}": r["upper"],
            f"ci_width_{label1}": r["ci_width"],
            f"mean_{label2}": r2.get("mean", float("nan")),
            f"lower_{label2}": r2.get("lower", float("nan")),
            f"upper_{label2}": r2.get("upper", float("nan")),
            f"ci_width_{label2}": r2.get("ci_width", float("nan")),
        })
    return result


def write_time_table_csv(table: list[dict], path: Path) -> None:
    if not table:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(table[0].keys()), delimiter=",")
        writer.writeheader()
        writer.writerows(table)


def _spearman(x: list[float], y: list[float]) -> tuple[float, float]:
    if len(x) < 3:
        return float("nan"), float("nan")
    rho, pval = spearmanr(x, y)
    return float(rho), float(pval)


def _linear_fit(x: list[float], y: list[float]) -> tuple[float, float, float]:
    if len(x) < 2:
        return float("nan"), float("nan"), float("nan")
    m, b = np.polyfit(x, y, 1)
    residuals = [yi - (m * xi + b) for xi, yi in zip(x, y)]
    rmse = float(np.sqrt(np.mean([r * r for r in residuals])))
    return float(m), float(b), rmse


def plot_convergence(
    table: list[dict],
    x_col: str,
    y_col: str,
    out_path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    node_field: str = "node",
) -> dict[str, float]:
    pairs = [
        (r[x_col], r[y_col], r.get(node_field, ""))
        for r in table
        if not np.isnan(r.get(x_col, float("nan")))
        and not np.isnan(r.get(y_col, float("nan")))
    ]
    x = [p[0] for p in pairs]
    y = [p[1] for p in pairs]
    labels = [p[2].replace("t_n", "n") if p[2].startswith("t_n") else p[2] for p in pairs]
    if not x:
        rho, pval, slope, intercept, rmse = float("nan"), float("nan"), float("nan"), float("nan"), float("nan")
    else:
        rho, pval = _spearman(x, y)
        slope, intercept, rmse = _linear_fit(x, y)

    all_vals = x + y
    lo, hi = (min(all_vals), max(all_vals)) if all_vals else (0, 1)
    pad = (hi - lo) * 0.06 if hi > lo else 0.5

    fig, ax = plt.subplots(figsize=(5, 5))
    if x:
        ax.scatter(x, y, s=40, alpha=0.8, color="steelblue", zorder=3)
        for xi, yi, lbl in zip(x, y, labels):
            ax.annotate(lbl, (xi, yi), textcoords="offset points",
                        xytext=(3, 3), fontsize=6, alpha=0.7)
    if len(x) >= 2:
        xs_line = np.linspace(lo - pad, hi + pad, 100)
        y_fit = slope * xs_line + intercept
        ax.plot(xs_line, y_fit, color="firebrick", lw=1.5, ls="--",
                label=_fit_label(slope, intercept))
        ax.legend(fontsize=8)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return {"rho": rho, "pvalue": pval, "slope": slope,
            "intercept": intercept, "rmse": rmse}


def _fit_label(slope: float, intercept: float) -> str:
    sign = "+" if intercept >= 0 else "-"
    return f"fit: y = {slope:.3f}x {sign} {abs(intercept):.3f}"


def plot_line(
    x: list[float],
    y: list[float],
    out_path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    labels: list[str] | None = None,
    square: bool = False,
) -> dict[str, float]:
    pairs = sorted(zip(x, y, labels or [""] * len(x)))
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    lbs = [p[2].replace("t_n", "n") if p[2].startswith("t_n") else p[2] for p in pairs]
    rho, pval = _spearman(xs, ys)
    slope, intercept, rmse = _linear_fit(xs, ys)

    all_vals = xs + ys
    lo, hi = min(all_vals), max(all_vals)
    pad = (hi - lo) * 0.06 if hi > lo else 0.5

    figsize = (5, 5) if square else (6, 4)
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(xs, ys, marker="o", markersize=4, color="steelblue", lw=1.2)
    for xi, yi, lb in zip(xs, ys, lbs):
        ax.annotate(lb, (xi, yi), textcoords="offset points",
                    xytext=(3, 3), fontsize=6, alpha=0.7)

    if len(xs) >= 2:
        x_line = np.linspace(lo - pad, hi + pad, 100) if square else np.linspace(min(xs), max(xs), 2)
        y_fit = slope * x_line + intercept
        ax.plot(x_line, y_fit, color="firebrick", lw=1.5, ls="--",
                label=_fit_label(slope, intercept))
        ax.legend(fontsize=8)

    if square:
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return {"rho": rho, "pvalue": pval, "slope": slope,
            "intercept": intercept, "rmse": rmse}


def plot_trace(mcmc_txt: Path, out_path: Path, title: str) -> None:
    if not mcmc_txt.exists():
        return
    lines = mcmc_txt.read_text(errors="ignore").splitlines()
    if len(lines) < 2:
        return
    header = lines[0].split()
    rows = []
    for line in lines[1:]:
        parts = line.split()
        try:
            rows.append([float(v) for v in parts])
        except ValueError:
            continue
    if not rows:
        return
    arr = np.array(rows)

    col_indices = [
        i for i, h in enumerate(header)
        if h.lower() not in ("gen", "iter", "time") and i < arr.shape[1]
    ]
    if not col_indices:
        col_indices = list(range(min(len(header), arr.shape[1], 12)))

    from matplotlib.backends.backend_pdf import PdfPages
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out_path) as pdf:
        for col_i in col_indices:
            col_name = header[col_i] if col_i < len(header) else f"col{col_i}"
            fig, ax = plt.subplots(figsize=(8, 3))
            ax.plot(arr[:, col_i], lw=0.8, color="steelblue")
            ax.set_title(col_name)
            ax.set_xlabel("iteration")
            ax.set_ylabel(col_name)
            fig.suptitle(title, fontsize=10)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def generate_all_diagnostics(
    *,
    run_dirs: list[Path],
    diag_dir: Path,
    n_runs: int = 2,
) -> dict[str, Any]:
    import itertools
    warnings: list[str] = []
    corr_rows: list[dict] = []
    output_files: dict[str, dict[str, str]] = {}

    def _add_corr(comparison: str, stats: dict[str, float]) -> None:
        row = {"comparison": comparison}
        for k in ("rho", "pvalue", "slope", "intercept", "rmse"):
            v = stats.get(k, float("nan"))
            row[k] = "" if (isinstance(v, float) and np.isnan(v)) else v
        corr_rows.append(row)
        if row["rho"] == "" or row["slope"] == "":
            warnings.append(
                f"{comparison}: too few valid points (need >= 3 for "
                f"Spearman, >= 2 for slope/intercept/RMSE)."
            )

    post_times: list[list[dict]] = []
    prior_times: list[list[dict]] = []
    generated: list[str] = []
    skipped: list[dict] = []

    for i, run_dir in enumerate(run_dirs):
        run_label = f"run{i+1}"

        out_file = run_dir / "mcmctree.out"
        if out_file.exists():
            rows = parse_mcmctree_out(out_file.read_text(errors="ignore"))
            post_times.append(rows)
            node_tree = extract_node_tree(out_file.read_text(errors="ignore"))
            if node_tree:
                (run_dir / "FigTree.node.tre").write_text(node_tree + "\n")
            else:
                skipped.append({"reason": "No FigTree tree found in mcmctree.out",
                                "file": str(out_file)})
        else:
            post_times.append([])

        prior_out = run_dir / "prior" / "mcmctree.out"
        if prior_out.exists():
            rows = parse_mcmctree_out(prior_out.read_text(errors="ignore"))
            prior_times.append(rows)
            node_tree = extract_node_tree(prior_out.read_text(errors="ignore"))
            if node_tree:
                (run_dir / "prior" / "FigTree.node.tre").write_text(node_tree + "\n")
            else:
                skipped.append({"reason": "No FigTree tree found in prior mcmctree.out",
                                "file": str(prior_out)})
        else:
            prior_times.append([])

        for kind, mcmc_file in [
            ("posterior", run_dir / "mcmc.txt"),
            ("prior", run_dir / "prior" / "mcmc.txt"),
        ]:
            trace_pdf = diag_dir / "traces" / f"mcmc_trace_{run_label}_{kind}.pdf"
            plot_trace(
                mcmc_file,
                trace_pdf,
                title=f"MCMC trace — {run_label} {kind}",
            )
            output_files[f"trace_{run_label}_{kind}"] = {"path": str(trace_pdf), "description": f"MCMC trace plot: {run_label} {kind}, parameter sampling over iterations"}

    if n_runs >= 2:
        for kind_label, times_list, out_name in [
            ("posterior", post_times, "posterior_times"),
            ("prior", prior_times, "prior_times"),
        ]:
            valid_runs = [(i, rows) for i, rows in enumerate(times_list) if rows]
            if len(valid_runs) >= 2:
                all_nodes: dict[str, dict] = {}
                for i, rows in valid_runs:
                    label = f"run{i+1}"
                    for r in rows:
                        node = r["node"]
                        if node not in all_nodes:
                            all_nodes[node] = {"node": node}
                        all_nodes[node][f"mean_{label}"] = r["mean"]
                        all_nodes[node][f"lower_{label}"] = r["lower"]
                        all_nodes[node][f"upper_{label}"] = r["upper"]
                        all_nodes[node][f"ci_width_{label}"] = r["ci_width"]
                combined_table = list(all_nodes.values())
                conv_csv = diag_dir / "convergence" / f"{out_name}.csv"
                write_time_table_csv(
                    combined_table,
                    conv_csv,
                )
                output_files[f"convergence_{out_name}"] = {"path": str(conv_csv), "description": f"Combined {kind_label} node age estimates from all valid runs"}

            pairs = list(itertools.combinations(range(len(times_list)), 2))
            for a, b in pairs:
                label_a = f"run{a+1}"
                label_b = f"run{b+1}"
                if (a < len(times_list) and b < len(times_list)
                        and times_list[a] and times_list[b]):
                    pair_table = build_time_table(times_list[a], times_list[b],
                                                  label1=label_a, label2=label_b)
                    conv_pdf = diag_dir / "convergence" / f"convergence_{kind_label}_{label_a}_vs_{label_b}.pdf"
                    stats = plot_convergence(
                        pair_table,
                        f"mean_{label_a}", f"mean_{label_b}",
                        conv_pdf,
                        title=f"Convergence — {kind_label} means ({label_a} vs {label_b})",
                        xlabel=f"Mean age {label_a} (100 Mya)",
                        ylabel=f"Mean age {label_b} (100 Mya)",
                    )
                    _add_corr(f"convergence_{kind_label}_{label_a}_vs_{label_b}", stats)
                    generated.append(f"convergence_{kind_label}_{label_a}_vs_{label_b}")
                    output_files[f"convergence_{kind_label}_{label_a}_vs_{label_b}"] = {"path": str(conv_pdf), "description": f"{kind_label.capitalize()} convergence diagnostic: scatter plot with regression for {label_a} vs {label_b}"}
                else:
                    skipped.append({"reason": f"{kind_label} mcmctree.out empty/missing",
                                    "run": f"{label_a}_vs_{label_b}"})
    else:
        skipped.append({"reason": "n_runs < 2; need >= 2 runs for convergence plots",
                        "n_runs": n_runs})

    for i, run_dir in enumerate(run_dirs):
        run_label = f"run{i+1}"
        post = post_times[i] if i < len(post_times) else []
        prior = prior_times[i] if i < len(prior_times) else []

        for kind, rows in [("posterior", post), ("prior", prior)]:
            if not rows:
                skipped.append({"reason": f"{kind} mcmctree.out empty for {run_label}",
                                "run": run_label})
                continue
            x = [r["mean"] for r in rows]
            y = [r["ci_width"] for r in rows]
            lbs = [r.get("node", "") for r in rows]
            inf_pdf = diag_dir / "infinite_sites" / f"infinite_sites_{run_label}_{kind}.pdf"
            stats = plot_line(
                x, y,
                inf_pdf,
                title=f"Infinite-sites — {run_label} {kind}",
                xlabel="Mean age (100 Mya)",
                ylabel="95% CI width (100 Mya)",
                labels=lbs,
            )
            _add_corr(f"infinite_sites_{run_label}_{kind}", stats)
            generated.append(f"infinite_sites_{run_label}_{kind}")
            output_files[f"infinite_sites_{run_label}_{kind}"] = {"path": str(inf_pdf), "description": f"Infinite-sites diagnostic: {run_label} {kind}, mean age vs 95% credible interval width"}

        if post and prior:
            post_by_node = {r["node"]: r["mean"] for r in post}
            prior_by_node = {r["node"]: r["mean"] for r in prior}
            shared = [n for n in post_by_node if n in prior_by_node]
            if shared:
                xp = [post_by_node[n] for n in shared]
                yp = [prior_by_node[n] for n in shared]
                pvp_pdf = diag_dir / "posterior_vs_prior" / f"posterior_vs_prior_{run_label}.pdf"
                stats = plot_line(
                    xp, yp,
                    pvp_pdf,
                    title=f"Posterior vs prior — {run_label}",
                    xlabel="Posterior mean age (100 Mya)",
                    ylabel="Prior mean age (100 Mya)",
                    labels=shared,
                    square=True,
                )
                _add_corr(f"posterior_vs_prior_{run_label}", stats)
                generated.append(f"posterior_vs_prior_{run_label}")
                output_files[f"posterior_vs_prior_{run_label}"] = {"path": str(pvp_pdf), "description": f"Posterior vs prior mean node age comparison for {run_label}"}
            else:
                skipped.append({"reason": "No shared nodes between posterior and prior",
                                "run": run_label})
        else:
            skipped.append({"reason": "Need both posterior and prior mcmctree.out",
                            "run": run_label})

    if corr_rows:
        corr_path = diag_dir / "spearman_correlations.csv"
        corr_path.parent.mkdir(parents=True, exist_ok=True)
        with open(corr_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["comparison", "rho", "pvalue",
                            "slope", "intercept", "rmse"],
            )
            writer.writeheader()
            writer.writerows(corr_rows)
        output_files["spearman_correlations"] = {"path": str(corr_path), "description": "Spearman rank correlations and regression statistics for convergence assessments"}

    return {"spearman": corr_rows, "warnings": warnings,
            "generated": generated, "skipped": skipped,
            "output_files": output_files}
