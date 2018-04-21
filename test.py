from motif_scanner import MotifScanner

scanner = MotifScanner("GAACGTCCCTGCCGTATTATTGTTTTTATAGG")
# background distribution
background = {
    "A": 0.25,
    "C": 0.25,
    "G": 0.25,
    "T": 0.25
}
# Position-specific Distribution Matrix
psdm = [
    {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.7},
    {'A': 0.7, 'C': 0.1, 'G': 0.1, 'T': 0.1},
    {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.7},
    {'A': 0.7, 'C': 0.1, 'G': 0.1, 'T': 0.1}
]
for p_value, position, piece in scanner.scan(p_value=0.05, background=background, psdm=psdm):
    # do something with p_value, position, piece, for example:
    print(p_value, position, piece)
