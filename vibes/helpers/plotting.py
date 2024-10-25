"""helpers for plotting"""

# Tableau colors:
tableau_colors_int = [
    (0, 107, 164),
    (255, 128, 14),
    (171, 171, 171),
    (89, 89, 89),
    (95, 158, 209),
    (200, 82, 0),
    (137, 137, 137),
]

tableau_colors = [tuple(c / 255 for c in color) for color in tableau_colors_int]


fontsize = 14

rc_params = {
    "font.family": "sans-serif",
    "font.sans-serif": "DejaVu Sans",
    "font.size": fontsize,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "savefig.dpi": 150,
    "figure.titlesize": fontsize + 2,
    "axes.facecolor": "white",  # axes background color
    "axes.edgecolor": "1b1b1b",  # axes edge color
    "axes.linewidth": 0.8,  # edge linewidth
    "axes.titlesize": fontsize,
    "xtick.major.size": 6,  # major tick size in points
    "xtick.minor.size": 4,  # minor tick size in points
    "xtick.major.width": 0.8,  # major tick width in points
    "xtick.minor.width": 0.6,  # minor tick width in points
    "xtick.color": "1b1b1b",  # color of the tick labels
    "ytick.major.size": 6,  # major tick size in points
    "ytick.minor.size": 4,  # minor tick size in points
    "ytick.major.width": 0.8,  # major tick width in points
    "ytick.minor.width": 0.6,  # minor tick width in points
    "ytick.color": "1b1b1b",  # color of the tick labels
}
