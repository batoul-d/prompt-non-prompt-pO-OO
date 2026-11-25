
from pathlib import Path
import base64
import html
import pandas as pd
import sys

# ============================================================
# CONFIGURATION
# ============================================================

if len(sys.argv) < 3:
    print("Usage: python makeReport.py PO_CASE_NAME OO_CASE_NAME [DRAW_DIAGNOSTICS]")
    sys.exit(1)

PO_CASE_NAME = sys.argv[1]
OO_CASE_NAME = sys.argv[2]

# Optional third argument: true/false
DRAW_DIAGNOSTICS = (
    len(sys.argv) > 3
    and sys.argv[3].lower() in ("true", "1", "yes", "y")
)

PO_PLOTS = Path(f"output/outputFits_pO_{PO_CASE_NAME}")
OO_PLOTS = Path(f"output/outputFits_OO_{OO_CASE_NAME}")

SUMMARY_PLOTS = Path("plots")

PO_TABLE = Path(f"output/fitSummary_pO_{PO_CASE_NAME}.tsv")
OO_TABLE = Path(f"output/fitSummary_OO_{OO_CASE_NAME}.tsv")

OUTPUT = Path("plots/comparison_report.html")

PLOT_WIDTH = 420

# Prefixes for NLL and sWeight plots
DIAGNOSTIC_PREFIXES = ("c0NLL", "cNLL", "cWeights")

print("PO directory:", PO_PLOTS.resolve())
print("PO exists:", PO_PLOTS.exists())
print("PO is directory:", PO_PLOTS.is_dir())

print("OO directory:", OO_PLOTS.resolve())
print("OO exists:", OO_PLOTS.exists())
print("OO is directory:", OO_PLOTS.is_dir())

print("Draw NLL and sWeight plots:", DRAW_DIAGNOSTICS)

if PO_PLOTS.exists():
    print("PO contents:", list(PO_PLOTS.iterdir())[:10])

if OO_PLOTS.exists():
    print("OO contents:", list(OO_PLOTS.iterdir())[:10])


# ============================================================
# HELPERS
# ============================================================

def load_plots(directory, diagnostics=False):
    """
    Load and sort PDF files.

    diagnostics=False:
        Normal fit plots, excluding NLL and sWeights.

    diagnostics=True:
        Only NLL and sWeight plots.

    Both exclude the integrated [0,20] bin.
    """

    if not directory.is_dir():
        print(f"[WARNING] Plot directory does not exist: {directory}")
        return []

    plots = []

    for p in directory.rglob("*"):
        if not p.is_file() or p.suffix.lower() != ".pdf":
            continue

        if "0_200" in p.name:
            continue

        is_diagnostic = p.stem.startswith(DIAGNOSTIC_PREFIXES)

        if diagnostics and is_diagnostic:
            plots.append(p)
        elif not diagnostics and not is_diagnostic:
            plots.append(p)

    return sorted(plots)


def image_to_html(path):
    """Embed PDF files directly into the HTML."""

    filename = html.escape(path.stem)
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")

    content = f"""
        <iframe
            src="data:application/pdf;base64,{encoded}"
            width="100%"
            height="500px">
        </iframe>
    """

    return f"""
    <div class="plot">
        {content}
        <div class="caption">{filename}</div>
    </div>
    """


def plots_to_html(directory, diagnostics=False):
    plots = load_plots(directory, diagnostics)

    if not plots:
        return "<p>No PDF files found.</p>"

    return "\n".join(image_to_html(p) for p in plots)


def table_to_html(path):
    if not path.is_file():
        print(f"[WARNING] Table does not exist: {path}")
        return f"<p>Table not found: {html.escape(str(path))}</p>"

    df = pd.read_csv(path, sep="\t")
    return df.to_html(index=False, border=0, classes="data-table")


def summary_plots_to_html(directory, prefix):
    if not directory.is_dir():
        return "<p>Summary plot directory not found.</p>"

    plots = sorted(
        p for p in directory.rglob("*")
        if p.is_file()
        and p.suffix.lower() == ".pdf"
        and p.stem.startswith(prefix + "_")
        and "0_200" not in p.name
    )

    if not plots:
        return "<p>No summary plots found.</p>"

    return "\n".join(image_to_html(p) for p in plots)


# ============================================================
# BUILD REPORT
# ============================================================

# Normal fit plots
po_plots_html = plots_to_html(PO_PLOTS, diagnostics=False)
oo_plots_html = plots_to_html(OO_PLOTS, diagnostics=False)

# NLL and sWeight plots (only if enabled)
if DRAW_DIAGNOSTICS:
    po_diagnostics_html = plots_to_html(PO_PLOTS, diagnostics=True)
    oo_diagnostics_html = plots_to_html(OO_PLOTS, diagnostics=True)

# Summary plots
po_summary_html = summary_plots_to_html(SUMMARY_PLOTS, "pO")
oo_summary_html = summary_plots_to_html(SUMMARY_PLOTS, "OO")

# Numerical tables
po_table_html = table_to_html(PO_TABLE)
oo_table_html = table_to_html(OO_TABLE)


# ============================================================
# OPTIONAL DIAGNOSTICS SECTION
# ============================================================

diagnostics_section = ""

if DRAW_DIAGNOSTICS:
    diagnostics_section = f"""
    <h2>NLL and sWeight plots</h2>

    <h3>pO NLL and sWeights — {html.escape(PO_CASE_NAME)}</h3>
    <div class="plot-row">
    {po_diagnostics_html}
    </div>

    <h3>OO NLL and sWeights — {html.escape(OO_CASE_NAME)}</h3>
    <div class="plot-row">
    {oo_diagnostics_html}
    </div>
    """


# ============================================================
# HTML REPORT
# ============================================================

report = f"""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>pO and OO comparison</title>

<style>

body {{
    font-family: Arial, sans-serif;
    margin: 30px;
    background: #ffffff;
}}

h1 {{
    margin-bottom: 35px;
}}

h2 {{
    margin-top: 35px;
    border-bottom: 2px solid #444;
    padding-bottom: 8px;
}}

h3 {{
    margin-top: 25px;
}}

.plot-row {{
    display: flex;
    flex-wrap: nowrap;
    gap: 15px;
    overflow-x: auto;
    padding: 15px 0;
}}

.plot {{
    flex: 0 0 {PLOT_WIDTH}px;
    text-align: center;
}}

.plot iframe {{
    border: none;
    width: 100%;
    height: 500px;
}}

.caption {{
    font-size: 12px;
    margin-top: 5px;
    color: #555;
}}

.tables {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 30px;
    margin-top: 20px;
}}

.table-container {{
    overflow-x: auto;
}}

.data-table {{
    border-collapse: collapse;
    font-size: 12px;
    white-space: nowrap;
}}

.data-table th,
.data-table td {{
    border: 1px solid #ccc;
    padding: 6px 9px;
    text-align: center;
}}

.data-table th {{
    background: #eeeeee;
}}

@media print {{
    .plot-row {{
        overflow: visible;
        flex-wrap: wrap;
    }}

    .plot {{
        flex: 0 0 23%;
    }}
}}

</style>
</head>

<body>

<h1>pO and OO comparison</h1>

<h2>Fit plots</h2>

<h3>pO fits — {html.escape(PO_CASE_NAME)}</h3>
<div class="plot-row">
{po_plots_html}
</div>

<h3>OO fits — {html.escape(OO_CASE_NAME)}</h3>
<div class="plot-row">
{oo_plots_html}
</div>

{diagnostics_section}

<h2>Summary plots</h2>

<h3>pO summary</h3>
<div class="plot-row">
{po_summary_html}
</div>

<h3>OO summary</h3>
<div class="plot-row">
{oo_summary_html}
</div>

<h2>Numerical results</h2>

<div class="tables">

    <div class="table-container">
        <h3>pO table</h3>
        {po_table_html}
    </div>

    <div class="table-container">
        <h3>OO table</h3>
        {oo_table_html}
    </div>

</div>

</body>
</html>
"""

# ============================================================
# SAVE REPORT
# ============================================================

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(report, encoding="utf-8")

print(f"\nReport created: {OUTPUT.resolve()}")

print(f"pO normal fits: {len(load_plots(PO_PLOTS, diagnostics=False))}")
print(f"OO normal fits: {len(load_plots(OO_PLOTS, diagnostics=False))}")

if DRAW_DIAGNOSTICS:
    print(f"pO NLL and sWeights: {len(load_plots(PO_PLOTS, diagnostics=True))}")
    print(f"OO NLL and sWeights: {len(load_plots(OO_PLOTS, diagnostics=True))}")

print(f"pO summary: {len(list(SUMMARY_PLOTS.glob('pO_*.pdf')))}")
print(f"OO summary: {len(list(SUMMARY_PLOTS.glob('OO_*.pdf')))}")
