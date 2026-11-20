"""Common figure style (PRE, REVTeX two-column). Figures are produced at final width."""
import matplotlib as mpl
mpl.use("Agg")
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "lines.linewidth": 1.0, "lines.markersize": 4,
    "mathtext.fontset": "stix", "font.family": "STIXGeneral",
    "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "legend.frameon": False,
})
COLW, TEXTW = 3.375, 7.0
# Okabe-Ito
BLUE, VERM, GREEN, ORANGE, SKY, PURPLE, YELLOW, GREY = "#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#999999"
C_POS, C_NEG = VERM, BLUE          # vacancy / charge on the two sublattices (same meaning in every figure)
C_FROZEN = "#DDDDDD"
def panel(ax, s, x=-0.02, y=1.02, **kw):
    ax.text(x, y, s, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, **kw)
def save(fig, name):
    fig.savefig(name + ".pdf"); fig.savefig(name + ".png", dpi=300)
def xcolor(x):
    """Colour for the scaling variable x = N*/N (same in every figure)."""
    return mpl.colormaps["viridis"](0.92 * (1 - x))
MARK = {0.1: "o", 0.2: "s", 0.3: "^", 0.5: "D", 0.7: "v", 0.9: "<"}
