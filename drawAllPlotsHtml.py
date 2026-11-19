
from pathlib import Path
import base64
import html
import pandas as pd

# ============================================================
# CONFIGURATION: change these paths
# ============================================================

PO_PLOTS = Path("output/outputFits_pO_fineBins_copy")
OO_PLOTS = Path("output/outputFits_OO_fineBins_VWG_4Gaussians")

SUMMARY_PLOTS = Path("plots")

PO_TABLE = Path("output/fitSummary_pO_fineBins_copy.tsv")
OO_TABLE = Path("output/fitSummary_OO_fineBins_VWG_4Gaussians.tsv")

OUTPUT = Path("plots/comparison_report.html")

PLOT_WIDTH = 420

print("PO directory:", PO_PLOTS.resolve())
print("PO exists:", PO_PLOTS.exists())
print("PO is directory:", PO_PLOTS.is_dir())

print("OO directory:", OO_PLOTS.resolve())
print("OO exists:", OO_PLOTS.exists())
print("OO is directory:", OO_PLOTS.is_dir())

print("PO contents:", list(PO_PLOTS.iterdir())[:10])
print("OO contents:", list(OO_PLOTS.iterdir())[:10])

# ============================================================
# HELPERS
# ============================================================

def load_plots(directory):
    """Load and sort PDF files, excluding [0,20] GeV/c^2 bin."""
    return sorted(
        p for p in directory.rglob("*")
        if p.is_file()
        and p.suffix.lower() in (".pdf")
        and "0_200" not in p.name
    )


def image_to_html(path):
    """Embed PDF files directly into the HTML."""

    filename = html.escape(path.stem)

    if path.suffix.lower() == ".pdf":
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


def plots_to_html(directory):
    plots = load_plots(directory)

    if not plots:
        return "<p>No PNG files found.</p>"

    return "\n".join(image_to_html(p) for p in plots)


def table_to_html(path):
    df = pd.read_csv(path, sep="\t")
    return df.to_html(index=False, border=0, classes="data-table")


def summary_plots_to_html(directory, prefix):
    plots = sorted(
        p for p in directory.rglob("*")
        if p.is_file()
        and p.suffix.lower() in (".pdf")
        and p.stem.startswith(prefix + "_")
        and "0_200" not in p.name
    )

    if not plots:
        return "<p>No plots found.</p>"

    return "\n".join(image_to_html(p) for p in plots)

# ============================================================
# BUILD REPORT
# ============================================================

# Fit plots
po_plots_html = plots_to_html(PO_PLOTS)
oo_plots_html = plots_to_html(OO_PLOTS)

# Summary plots
po_summary_html = summary_plots_to_html(SUMMARY_PLOTS, "pO")
oo_summary_html = summary_plots_to_html(SUMMARY_PLOTS, "OO")

# Numerical tables
po_table_html = table_to_html(PO_TABLE)
oo_table_html = table_to_html(OO_TABLE)


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

.plot img {{
    width: 100%;
    height: auto;
    object-fit: contain;
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

<h3>pO fits</h3>
<div class="plot-row">
{po_plots_html}
</div>

<h3>OO fits</h3>
<div class="plot-row">
{oo_plots_html}
</div>


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

OUTPUT.write_text(report, encoding="utf-8")

print(f"Report created: {OUTPUT.resolve()}")
print(f"pO fits: {len(load_plots(PO_PLOTS))}")
print(f"OO fits: {len(load_plots(OO_PLOTS))}")
print(f"pO summary: {len(list(SUMMARY_PLOTS.glob('pO_*.pdf')))}")
print(f"OO summary: {len(list(SUMMARY_PLOTS.glob('OO_*.pdf')))}")