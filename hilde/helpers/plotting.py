""" helpers for plotting """

# Tableau colors:
tableau_colors_int = [
    (171, 171, 171),
    (0, 107, 164),
    (255, 128, 14),
    (137, 137, 137),
    (89, 89, 89),
    (200, 82, 0),
]

tableau_colors = [
    tuple(c / 255 for c in color) for color in tableau_colors_int
]