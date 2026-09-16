"""Gradient interlayer screening pipeline.

The script is self-contained. It reads ``parameters.json``, runs the
thermodynamic screening, applies the CTE filter, and writes a single
integrated ``results.csv`` table.
"""

import csv
import json
import math
import os
import sys


REPO_DIR = os.path.dirname(os.path.abspath(__file__))
PARAMETERS_PATH = os.path.join(REPO_DIR, "parameters.json")
RESULTS_PATH = os.path.join(REPO_DIR, "results.csv")
COMPOSITION_STEP = 0.005


ELEMENTS = {
    "H":  {"phi": 5.20, "n_ws": 1.50, "v23": 1.42,"K": None, "G": None, "Tm": 14, "etype": "NTM", "Rp": 0.0, "dH_trans": 100, "valence": 1, "z": 1},
    "Li": {"phi": 2.85, "n_ws": 0.96, "v23": 5.5, "K": 11.6, "G": 4.3, "Tm": 454, "etype": "NTM", "Rp": 0.0, "dH_trans": 0, "valence": 1, "z": 3},
    "Be": {"phi": 4.50, "n_ws": 1.67, "v23": 2.9, "K": 130, "G": 150, "Tm": 1560, "etype": "NTM", "Rp": 0.4, "dH_trans": 0, "valence": 2, "z": 4},
    "B":  {"phi": 5.30, "n_ws": 1.55, "v23": 2.8, "K": 200, "G": 170, "Tm": 2349, "etype": "NTM", "Rp": 0.0, "dH_trans": 30, "valence": 3, "z": 5},
    "C":  {"phi": 6.20, "n_ws": 1.80, "v23": 2.0, "K": 442, "G": 534, "Tm": 3823, "etype": "NTM", "Rp": 0.0, "dH_trans": 180, "valence": 4, "z": 6},
    "N":  {"phi": 7.30, "n_ws": 1.90, "v23": 1.6, "K": None, "G": None, "Tm": 63, "etype": "NTM", "Rp": 0.0, "dH_trans": 310, "valence": 5, "z": 7},
    "Na": {"phi": 2.70, "n_ws": 0.82, "v23": 8.3, "K": 6.7, "G": 2.5, "Tm": 371, "etype": "NTM", "Rp": 0.0, "dH_trans": 0, "valence": 1, "z": 11},
    "Mg": {"phi": 3.45, "n_ws": 1.17, "v23": 5.8, "K": 35.6, "G": 17.3, "Tm": 923, "etype": "NTM", "Rp": 0.4, "dH_trans": 0, "valence": 2, "z": 12},
    "Al": {"phi": 4.20, "n_ws": 1.39, "v23": 4.6, "K": 73.5, "G": 26.9, "Tm": 933, "etype": "NTM", "Rp": 1.9, "dH_trans": 0, "valence": 3, "z": 13},
    "Si": {"phi": 4.70, "n_ws": 1.50, "v23": 4.2, "K": 99.0, "G": 67.0, "Tm": 1687, "etype": "NTM", "Rp": 2.1, "dH_trans": 34, "valence": 4, "z": 14},
    "P":  {"phi": 5.55, "n_ws": 1.60, "v23": 4.0, "K": 35.0, "G": 17.0, "Tm": 317, "etype": "NTM", "Rp": 0.0, "dH_trans": 17, "valence": 5, "z": 15},
    "K":  {"phi": 2.25, "n_ws": 0.65, "v23": 12.8, "K": 3.3, "G": 1.3, "Tm": 337, "etype": "NTM", "Rp": 0.0, "dH_trans": 0, "valence": 1, "z": 19},
    "Ca": {"phi": 2.55, "n_ws": 0.91, "v23": 8.8, "K": 15.2, "G": 7.4, "Tm": 1115, "etype": "NTM", "Rp": 0.4, "dH_trans": 0, "valence": 2, "z": 20},
    "Sc": {"phi": 3.25, "n_ws": 1.27, "v23": 6.1, "K": 55.0, "G": 30.0, "Tm": 1814, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 21},
    "Ti": {"phi": 3.80, "n_ws": 1.52, "v23": 4.8, "K": 108.4, "G": 43.4, "Tm": 1941, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 4, "z": 22},
    "V":  {"phi": 4.25, "n_ws": 1.64, "v23": 4.1, "K": 161.9, "G": 47.3, "Tm": 2183, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 5, "z": 23},
    "Cr": {"phi": 4.65, "n_ws": 1.73, "v23": 3.7, "K": 190.3, "G": 115.3, "Tm": 2180, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 6, "z": 24},
    "Mn": {"phi": 4.45, "n_ws": 1.61, "v23": 3.8, "K": 131, "G": 82, "Tm": 1519, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 7, "z": 25},
    "Fe": {"phi": 4.93, "n_ws": 1.77, "v23": 3.7, "K": 168.3, "G": 81.8, "Tm": 1811, "etype": "TM", "Rp": 0.85, "dH_trans": 0, "valence": 8, "z": 26},
    "Co": {"phi": 5.40, "n_ws": 1.79, "v23": 3.5, "K": 191.4, "G": 81.8, "Tm": 1768, "etype": "TM", "Rp": 0.85, "dH_trans": 0, "valence": 9, "z": 27},
    "Ni": {"phi": 5.30, "n_ws": 1.75, "v23": 3.5, "K": 187.6, "G": 94.7, "Tm": 1728, "etype": "TM", "Rp": 0.85, "dH_trans": 0, "valence": 10, "z": 28},
    "Cu": {"phi": 4.45, "n_ws": 1.47, "v23": 3.7, "K": 137.8, "G": 48.4, "Tm": 1358, "etype": "TM", "Rp": 0.3, "dH_trans": 0, "valence": 11, "z": 29},
    "Zn": {"phi": 3.95, "n_ws": 1.32, "v23": 4.4, "K": 64.8, "G": 37.9, "Tm": 693, "etype": "NTM", "Rp": 1.4, "dH_trans": 0, "valence": 2, "z": 30},
    "Ga": {"phi": 4.10, "n_ws": 1.31, "v23": 5.2, "K": 56.9, "G": 23.0, "Tm": 303, "etype": "NTM", "Rp": 1.9, "dH_trans": 0, "valence": 3, "z": 31},
    "Ge": {"phi": 4.60, "n_ws": 1.37, "v23": 4.6, "K": 74.9, "G": 53.4, "Tm": 1211, "etype": "NTM", "Rp": 2.1, "dH_trans": 25, "valence": 4, "z": 32},
    "As": {"phi": 4.80, "n_ws": 1.44, "v23": 5.2, "K": 35.0, "G": 17.0, "Tm": 1090, "etype": "NTM", "Rp": 0.0, "dH_trans": 0, "valence": 5, "z": 33},
    "Se": {"phi": 4.80, "n_ws": 1.40, "v23": 5.0, "K": 9.0, "G": 4.0, "Tm": 494, "etype": "NTM", "Rp": 0.0, "dH_trans": 0, "valence": 6, "z": 34},
    "Rb": {"phi": 2.10, "n_ws": 0.60, "v23": 14.6, "K": 2.9, "G": 1.1, "Tm": 312, "etype": "NTM", "Rp": 0.0, "dH_trans": 0, "valence": 1, "z": 37},
    "Sr": {"phi": 2.40, "n_ws": 0.84, "v23": 10.2, "K": 11.0, "G": 5.5, "Tm": 1050, "etype": "NTM", "Rp": 0.4, "dH_trans": 0, "valence": 2, "z": 38},
    "Y":  {"phi": 3.20, "n_ws": 1.21, "v23": 7.3, "K": 40.8, "G": 26.0, "Tm": 1799, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 39},
    "Zr": {"phi": 3.45, "n_ws": 1.39, "v23": 5.8, "K": 95.3, "G": 35.1, "Tm": 2128, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 4, "z": 40},
    "Nb": {"phi": 4.05, "n_ws": 1.64, "v23": 4.9, "K": 168.4, "G": 37.4, "Tm": 2750, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 5, "z": 41},
    "Mo": {"phi": 4.65, "n_ws": 1.77, "v23": 4.4, "K": 261.2, "G": 122.6, "Tm": 2896, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 6, "z": 42},
    "Ru": {"phi": 5.42, "n_ws": 1.83, "v23": 3.8, "K": 320, "G": 173, "Tm": 2607, "etype": "TM", "Rp": 0.0, "dH_trans": 0, "valence": 8, "z": 44},
    "Rh": {"phi": 5.55, "n_ws": 1.82, "v23": 3.8, "K": 275, "G": 150, "Tm": 2237, "etype": "TM", "Rp": 0.0, "dH_trans": 0, "valence": 9, "z": 45},
    "Pd": {"phi": 5.60, "n_ws": 1.67, "v23": 4.3, "K": 187.1, "G": 45.0, "Tm": 1828, "etype": "TM", "Rp": 0.0, "dH_trans": 0, "valence": 10, "z": 46},
    "Ag": {"phi": 4.35, "n_ws": 1.36, "v23": 4.7, "K": 103.6, "G": 30.3, "Tm": 1235, "etype": "TM", "Rp": 0.15, "dH_trans": 0, "valence": 11, "z": 47},
    "Cd": {"phi": 4.05, "n_ws": 1.24, "v23": 5.5, "K": 48.7, "G": 22.2, "Tm": 594, "etype": "NTM", "Rp": 1.4, "dH_trans": 0, "valence": 2, "z": 48},
    "In": {"phi": 3.90, "n_ws": 1.17, "v23": 6.3, "K": 41.5, "G": 8.7, "Tm": 430, "etype": "NTM", "Rp": 1.9, "dH_trans": 0, "valence": 3, "z": 49},
    "Sn": {"phi": 4.15, "n_ws": 1.24, "v23": 6.4, "K": 54.1, "G": 18.5, "Tm": 505, "etype": "NTM", "Rp": 2.1, "dH_trans": 0, "valence": 4, "z": 50},
    "Sb": {"phi": 4.40, "n_ws": 1.26, "v23": 6.6, "K": 39.0, "G": 13.0, "Tm": 904, "etype": "NTM", "Rp": 0.0, "dH_trans": 0, "valence": 5, "z": 51},
    "Te": {"phi": 4.65, "n_ws": 1.32, "v23": 5.6, "K": 23.0, "G": 11.0, "Tm": 723, "etype": "NTM", "Rp": 0.0, "dH_trans": 0, "valence": 6, "z": 52},
    "Cs": {"phi": 1.95, "n_ws": 0.55, "v23": 16.8, "K": 2.1, "G": 0.8, "Tm": 302, "etype": "NTM", "Rp": 0.0, "dH_trans": 0, "valence": 1, "z": 55},
    "Ba": {"phi": 2.32, "n_ws": 0.81, "v23": 11.0, "K": 9.6, "G": 4.8, "Tm": 1000, "etype": "NTM", "Rp": 0.4, "dH_trans": 0, "valence": 2, "z": 56},
    "La": {"phi": 3.20, "n_ws": 1.18, "v23": 8.0, "K": 27.9, "G": 14.5, "Tm": 1193, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 57},
    "Ce": {"phi": 3.18, "n_ws": 1.19, "v23": 7.5, "K": 21.5, "G": 12.0, "Tm": 1071, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 58},
    "Pr": {"phi": 3.19, "n_ws": 1.20, "v23": 7.5, "K": 28.8, "G": 14.8, "Tm": 1208, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 59},
    "Nd": {"phi": 3.19, "n_ws": 1.20, "v23": 7.4, "K": 31.8, "G": 15.8, "Tm": 1297, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 60},
    "Sm": {"phi": 3.20, "n_ws": 1.21, "v23": 7.3, "K": 37.8, "G": 19.6, "Tm": 1345, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 62},
    "Eu": {"phi": 2.50, "n_ws": 0.88, "v23": 9.5, "K": 11.5, "G": 5.6, "Tm": 1095, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 2, "z": 63},
    "Gd": {"phi": 3.22, "n_ws": 1.21, "v23": 7.1, "K": 38.8, "G": 22.0, "Tm": 1585, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 64},
    "Tb": {"phi": 3.22, "n_ws": 1.22, "v23": 7.0, "K": 40.0, "G": 23.0, "Tm": 1629, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 65},
    "Dy": {"phi": 3.22, "n_ws": 1.22, "v23": 6.9, "K": 40.5, "G": 25.8, "Tm": 1685, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 66},
    "Ho": {"phi": 3.22, "n_ws": 1.23, "v23": 6.9, "K": 40.2, "G": 27.0, "Tm": 1747, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 67},
    "Er": {"phi": 3.22, "n_ws": 1.23, "v23": 6.8, "K": 44.4, "G": 28.3, "Tm": 1802, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 68},
    "Tm": {"phi": 3.22, "n_ws": 1.24, "v23": 6.8, "K": 44.7, "G": 30.5, "Tm": 1818, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 69},
    "Yb": {"phi": 2.55, "n_ws": 0.92, "v23": 8.5, "K": 13.3, "G": 6.9, "Tm": 1097, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 2, "z": 70},
    "Lu": {"phi": 3.22, "n_ws": 1.24, "v23": 6.7, "K": 47.6, "G": 27.2, "Tm": 1936, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 3, "z": 71},
    "Hf": {"phi": 3.60, "n_ws": 1.45, "v23": 5.6, "K": 109, "G": 55.6, "Tm": 2506, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 4, "z": 72},
    "Ta": {"phi": 4.05, "n_ws": 1.63, "v23": 4.9, "K": 196.3, "G": 69.2, "Tm": 3290, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 5, "z": 73},
    "W":  {"phi": 4.80, "n_ws": 1.81, "v23": 4.5, "K": 310, "G": 160, "Tm": 3695, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 6, "z": 74},
    "Re": {"phi": 5.20, "n_ws": 1.85, "v23": 4.2, "K": 370, "G": 178, "Tm": 3459, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 7, "z": 75},
    "Os": {"phi": 6.00, "n_ws": 1.90, "v23": 4.0, "K": 418, "G": 222, "Tm": 3306, "etype": "TM", "Rp": 0.0, "dH_trans": 0, "valence": 8, "z": 76},
    "Ir": {"phi": 5.55, "n_ws": 1.90, "v23": 4.0, "K": 371, "G": 210, "Tm": 2719, "etype": "TM", "Rp": 0.0, "dH_trans": 0, "valence": 9, "z": 77},
    "Pt": {"phi": 5.65, "n_ws": 1.78, "v23": 4.4, "K": 278, "G": 61, "Tm": 2041, "etype": "TM", "Rp": 0.0, "dH_trans": 0, "valence": 10, "z": 78},
    "Au": {"phi": 5.15, "n_ws": 1.57, "v23": 4.7, "K": 173.2, "G": 27.8, "Tm": 1337, "etype": "TM", "Rp": 0.3, "dH_trans": 0, "valence": 11, "z": 79},
    "Hg": {"phi": 4.20, "n_ws": 1.24, "v23": 5.8, "K": 27.8, "G": 13.0, "Tm": 234, "etype": "NTM", "Rp": 1.4, "dH_trans": 0, "valence": 2, "z": 80},
    "Tl": {"phi": 3.90, "n_ws": 1.12, "v23": 6.6, "K": 35.7, "G": 5.7, "Tm": 577, "etype": "NTM", "Rp": 1.9, "dH_trans": 0, "valence": 3, "z": 81},
    "Pb": {"phi": 4.10, "n_ws": 1.15, "v23": 6.9, "K": 43.0, "G": 8.6, "Tm": 601, "etype": "NTM", "Rp": 2.1, "dH_trans": 0, "valence": 4, "z": 82},
    "Bi": {"phi": 4.15, "n_ws": 1.16, "v23": 7.2, "K": 33.2, "G": 8.3, "Tm": 545, "etype": "NTM", "Rp": 0.0, "dH_trans": 0, "valence": 5, "z": 83},
    "Th": {"phi": 3.30, "n_ws": 1.28, "v23": 7.3, "K": 54.5, "G": 31.0, "Tm": 2023, "etype": "TM", "Rp": 0.7, "dH_trans": 0, "valence": 4, "z": 90},
    "U":  {"phi": 4.05, "n_ws": 1.56, "v23": 5.6, "K": 113, "G": 74, "Tm": 1408, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 6, "z": 92},
    "Pu": {"phi": 3.60, "n_ws": 1.44, "v23": 5.2, "K": 97, "G": 43, "Tm": 913, "etype": "TM", "Rp": 1.0, "dH_trans": 0, "valence": 6, "z": 94},
}


ATOMIC_RADIUS = {
    "H":0.53,"Li":1.52,"Be":1.12,"B":0.87,"C":0.77,"N":0.71,"Na":1.86,"Mg":1.60,
    "Al":1.43,"Si":1.17,"P":1.10,"K":2.27,"Ca":1.97,"Sc":1.64,"Ti":1.46,"V":1.34,
    "Cr":1.28,"Mn":1.32,"Fe":1.27,"Co":1.25,"Ni":1.25,"Cu":1.28,"Zn":1.37,
    "Ga":1.35,"Ge":1.22,"As":1.21,"Se":1.22,"Rb":2.48,"Sr":2.15,"Y":1.78,
    "Zr":1.60,"Nb":1.46,"Mo":1.39,"Tc":1.36,"Ru":1.34,"Rh":1.35,"Pd":1.38,
    "Ag":1.44,"Cd":1.57,"In":1.66,"Sn":1.58,"Sb":1.45,"Te":1.35,"Cs":2.65,
    "Ba":2.22,"La":1.88,"Hf":1.59,"Ta":1.46,"W":1.39,"Re":1.37,"Os":1.35,
    "Ir":1.36,"Pt":1.39,"Au":1.44,"Hg":1.60,"Tl":1.72,"Pb":1.75,"Bi":1.63,
    "Th":1.80,"U":1.55,"Pu":1.60,
    "Ce":1.83,"Pr":1.82,"Nd":1.82,"Pm":1.80,"Sm":1.80,"Eu":2.04,"Gd":1.80,
    "Tb":1.78,"Dy":1.77,"Ho":1.76,"Er":1.75,"Tm":1.75,"Yb":1.94,"Lu":1.73,
    "Pa":1.61,"Np":1.55,"Am":1.73,
}


R_GAS = 0.008314
QP_RATIO = 9.4
P_TM_TM = 14.2
P_NTM_NTM = 10.7
P_TM_NTM = 12.35
MAX_VOL_ITER = 10
VOL_TOL = 1e-6
LATTICE_STABILITY = {
    3: [0.0, 0.0, 0.05],
    4: [0.0, 0.01, 0.10],
    5: [0.08, 0.0, 0.0],
    6: [0.10, 0.05, 0.0],
    7: [0.15, 0.10, 0.02],
    8: [0.02, 0.10, 0.0],
    9: [0.02, 0.0, 0.03],
    10: [0.0, 0.03, 0.05],
    11: [0.0, 0.05, 0.10],
}


def has_elastic_data(symbol):
    element = ELEMENTS.get(symbol)
    return bool(element and element.get("K") and element.get("G"))


def _get_pq(ea, eb):
    da, db = ELEMENTS[ea], ELEMENTS[eb]
    if da["etype"] == "TM" and db["etype"] == "TM":
        return P_TM_TM, P_TM_TM * QP_RATIO, False
    if da["etype"] == "NTM" and db["etype"] == "NTM":
        return P_NTM_NTM, P_NTM_NTM * QP_RATIO, False
    return P_TM_NTM, P_TM_NTM * QP_RATIO, True


def _calc_R(ea, eb, p_value, transition_metal):
    if not transition_metal:
        return 0.0
    return p_value * ELEMENTS[ea]["Rp"] * ELEMENTS[eb]["Rp"]


def calc_gamma(ea, eb):
    da, db = ELEMENTS[ea], ELEMENTS[eb]
    dp = da["phi"] - db["phi"]
    dn = da["n_ws"] - db["n_ws"]
    ni = 0.5 * (1.0 / da["n_ws"] + 1.0 / db["n_ws"])
    p_value, q_value, tm = _get_pq(ea, eb)
    return (-p_value * dp * dp + q_value * dn * dn - _calc_R(ea, eb, p_value, tm)) / ni


def volume_correction(ea, eb, ca, cb):
    da, db = ELEMENTS[ea], ELEMENTS[eb]
    vap, vbp = da["v23"], db["v23"]
    na, nb = da["n_ws"], db["n_ws"]
    dp = da["phi"] - db["phi"]
    va, vb = vap, vbp
    for _ in range(MAX_VOL_ITER):
        total = ca * va + cb * vb
        cas = 0.5 if total < 1e-12 else ca * va / total
        cbs = 0.5 if total < 1e-12 else cb * vb / total
        van = vap + (1.5 / na) * dp * cbs
        vbn = vbp + (1.5 / nb) * (-dp) * cas
        if abs(van - va) < VOL_TOL and abs(vbn - vb) < VOL_TOL:
            va, vb = van, vbn
            break
        va, vb = van, vbn
    total = ca * va + cb * vb
    if total < 1e-12:
        return va, vb, 0.5, 0.5
    return va, vb, ca * va / total, cb * vb / total


def _calc_elastic(ea, eb, ca, cb, va, vb):
    ka, ga = ELEMENTS[ea]["K"], ELEMENTS[ea]["G"]
    kb, gb = ELEMENTS[eb]["K"], ELEMENTS[eb]["G"]
    wa, wb = va ** 1.5, vb ** 1.5

    def _energy(ks, gm, ws, wm):
        numerator = 2.0 * ks * gm * (ws - wm) ** 2
        denominator = 3.0 * ks * wm + 4.0 * gm * ws
        return 0.0 if denominator < 1e-12 else numerator / denominator

    return ca * cb * (
        cb * _energy(ka, gb, wa, wb) + ca * _energy(kb, ga, wb, wa)
    )


def _calc_structural(ea, eb, ca, cb):
    za = ELEMENTS[ea]["valence"]
    zb = ELEMENTS[eb]["valence"]
    zk = min(max(round(ca * za + cb * zb), 3), 11)
    stability = LATTICE_STABILITY.get(zk)
    if stability is None:
        return 0.0
    energy_mix = min(stability)
    sa = LATTICE_STABILITY.get(min(max(za, 3), 11))
    sb = LATTICE_STABILITY.get(min(max(zb, 3), 11))
    energy_ref = ca * min(sa) + cb * min(sb) if sa and sb else 0.0
    return (energy_mix - energy_ref) * 96.485


def calc_binary_solid_solution(ea, eb, ca):
    cb = 1.0 - ca
    gamma = calc_gamma(ea, eb)
    va, vb, cas, cbs = volume_correction(ea, eb, ca, cb)
    dc = ca * cb * (cbs * va + cas * vb) * gamma
    de = (
        _calc_elastic(ea, eb, ca, cb, va, vb)
        if has_elastic_data(ea) and has_elastic_data(eb)
        else 0.0
    )
    return dc + de + _calc_structural(ea, eb, ca, cb)


def calc_delta(elements, fractions):
    r_avg = sum(f * ATOMIC_RADIUS.get(el, 1.4) for el, f in zip(elements, fractions))
    if r_avg < 1e-12:
        return 0.0
    delta_sq = sum(
        f * (1 - ATOMIC_RADIUS.get(el, 1.4) / r_avg) ** 2
        for el, f in zip(elements, fractions)
    )
    return math.sqrt(delta_sq) * 100.0


def composition_fractions(step=COMPOSITION_STEP):
    step_count = int(round(1.0 / step))
    return [index / step_count for index in range(1, step_count)]


def mixture_enthalpy_delta(end, candidate, fraction):
    composition = {}
    for element, share in zip(end["elements"], end["fractions"]):
        composition[element] = composition.get(element, 0.0) + share * (1.0 - fraction)
    composition[candidate] = composition.get(candidate, 0.0) + fraction

    elements = list(composition.keys())
    fractions = [composition[element] for element in elements]
    if len(elements) < 2:
        return 0.0, 0.0

    enthalpy = 0.0
    for i in range(len(elements)):
        for j in range(i + 1, len(elements)):
            pair_sum = fractions[i] + fractions[j]
            x_i = fractions[i] / max(pair_sum, 1e-12)
            enthalpy += pair_sum * pair_sum * calc_binary_solid_solution(
                elements[i], elements[j], x_i
            )
    delta = calc_delta(elements, fractions)
    return enthalpy, delta


def end_extremes(end, candidate, step=COMPOSITION_STEP):
    worst_enthalpy = None
    worst_delta = 0.0
    for fraction in composition_fractions(step):
        enthalpy, delta = mixture_enthalpy_delta(end, candidate, fraction)
        if worst_enthalpy is None or abs(enthalpy) > abs(worst_enthalpy):
            worst_enthalpy = enthalpy
        worst_delta = max(worst_delta, delta)
    return worst_enthalpy, worst_delta


def binary_extremes(element_x, element_y, step=COMPOSITION_STEP):
    worst_enthalpy = None
    worst_delta = 0.0
    for fraction in composition_fractions(step):
        enthalpy = calc_binary_solid_solution(element_x, element_y, fraction)
        delta = calc_delta([element_x, element_y], [fraction, 1.0 - fraction])
        if worst_enthalpy is None or abs(enthalpy) > abs(worst_enthalpy):
            worst_enthalpy = enthalpy
        worst_delta = max(worst_delta, delta)
    return worst_enthalpy, worst_delta


def passes(enthalpy, delta, thresholds):
    return (
        thresholds["dh_min"] <= enthalpy <= thresholds["dh_max"]
        and delta <= thresholds["delta_max"]
    )


def scan_end(end, candidates, thresholds, step=COMPOSITION_STEP):
    rows = []
    for candidate in candidates:
        enthalpy, delta = end_extremes(end, candidate, step)
        rows.append(
            {
                "candidate": candidate,
                "worst_dh": enthalpy,
                "abs_worst_dh": abs(enthalpy),
                "worst_delta": delta,
                "passed": passes(enthalpy, delta, thresholds),
            }
        )
    rows.sort(key=lambda row: (row["abs_worst_dh"], row["worst_delta"]))
    return rows


def cross_paths(end_a_candidates, end_b_candidates, thresholds, step=COMPOSITION_STEP):
    paths = []
    for a_row in end_a_candidates:
        for b_row in end_b_candidates:
            x = a_row["candidate"]
            y = b_row["candidate"]
            enthalpy, delta = binary_extremes(x, y, step)
            if passes(enthalpy, delta, thresholds):
                paths.append(
                    {
                        "x": x,
                        "y": y,
                        "worst_dh": enthalpy,
                        "abs_worst_dh": abs(enthalpy),
                        "worst_delta": delta,
                    }
                )
    paths.sort(key=lambda row: (row["x"], row["y"]))
    return paths


def _path_cte(x, y, end_a_name, end_b_name, alpha_a, alpha_b, candidate_cte):
    if x == y:
        nodes = [end_a_name, x, end_b_name]
        alphas = [alpha_a, candidate_cte[x], alpha_b]
    else:
        nodes = [end_a_name, x, y, end_b_name]
        alphas = [alpha_a, candidate_cte[x], candidate_cte[y], alpha_b]
    return nodes, alphas


def _elin(alphas):
    node_count = len(alphas)
    interval_count = node_count - 1
    ideal = [
        alphas[0] + (index / interval_count) * (alphas[-1] - alphas[0])
        for index in range(node_count)
    ]
    squared_error = sum(
        (actual - target) ** 2 for actual, target in zip(alphas, ideal)
    )
    return math.sqrt(squared_error / node_count)


def cte_screen(paths, end_a_name, end_b_name, alpha_a, alpha_b, candidate_cte):
    ranked = []
    for path in paths:
        nodes, alphas = _path_cte(
            path["x"],
            path["y"],
            end_a_name,
            end_b_name,
            alpha_a,
            alpha_b,
            candidate_cte,
        )
        if any(alphas[i] >= alphas[i + 1] for i in range(len(alphas) - 1)):
            continue
        ranked.append(
            {
                "x": path["x"],
                "y": path["y"],
                "path": "-".join(nodes),
                "n_nodes": len(nodes),
                "alpha_values": alphas,
                "elin": _elin(alphas),
            }
        )
    ranked.sort(key=lambda row: (row["elin"], row["path"]))
    return ranked


def load_parameters(path=None):
    with open(path or PARAMETERS_PATH, "r", encoding="utf-8") as handle:
        parameters = json.load(handle)
    _validate_parameters(parameters)
    return parameters


def _validate_parameters(parameters):
    for key in ("end_a", "end_b"):
        end = parameters.get(key)
        if not isinstance(end, dict):
            raise ValueError(f"`{key}` must be an object")
        if not isinstance(end.get("name"), str) or not end["name"]:
            raise ValueError(f"`{key}.name` must be a non-empty string")
        if not isinstance(end.get("elements"), list) or not end["elements"]:
            raise ValueError(f"`{key}.elements` must be a non-empty list")
        if not isinstance(end.get("fractions"), list) or len(end["fractions"]) != len(
            end["elements"]
        ):
            raise ValueError(f"`{key}.fractions` must have the same length as elements")
        for element in end["elements"]:
            if element not in ELEMENTS:
                raise ValueError(f"unknown element `{element}` in `{key}`")
        if abs(sum(end["fractions"]) - 1.0) > 1e-6:
            raise ValueError(f"`{key}.fractions` must sum to 1")
        if "cte" in end and end["cte"] is not None:
            float(end["cte"])

    thresholds = parameters.get("thresholds", {})
    if thresholds.get("dh_min", 0.0) >= thresholds.get("dh_max", 0.0):
        raise ValueError("`thresholds.dh_min` must be smaller than `dh_max`")
    if thresholds.get("delta_max", 0.0) < 0.0:
        raise ValueError("`thresholds.delta_max` must be non-negative")

    candidates = parameters.get("candidates", {})
    include = candidates.get("include")
    exclude = candidates.get("exclude", [])
    if include is not None and not isinstance(include, list):
        raise ValueError("`candidates.include` must be a list or null")
    if not isinstance(exclude, list):
        raise ValueError("`candidates.exclude` must be a list")


def candidate_symbols(parameters):
    candidates = parameters.get("candidates", {})
    include = candidates.get("include")
    exclude = set(candidates.get("exclude", []))
    symbols = list(ELEMENTS.keys()) if include is None else list(include)
    for symbol in list(symbols) + list(exclude):
        if symbol not in ELEMENTS:
            raise ValueError(f"unknown candidate element `{symbol}`")
    return [symbol for symbol in symbols if symbol not in exclude]


def clean_value(value, digits=4):
    if value is None:
        return None
    rounded_value = round(value, digits)
    return 0.0 if abs(rounded_value) < 10 ** (-digits) else rounded_value


def _empty_row(section):
    return {
        "section": section,
        "id": None,
        "path": None,
        "candidate": None,
        "x": None,
        "y": None,
        "worst_dh": None,
        "abs_worst_dh": None,
        "worst_delta": None,
        "rank": None,
        "n_nodes": None,
        "alpha_0": None,
        "alpha_1": None,
        "alpha_2": None,
        "alpha_3": None,
        "elin": None,
    }


def build_results(end_a_candidates, end_b_candidates, paths, ranked_paths):
    rows = []
    for row in end_a_candidates:
        item = _empty_row("end_a")
        item["id"] = row["candidate"]
        item["candidate"] = row["candidate"]
        item["worst_dh"] = clean_value(row["worst_dh"], 4)
        item["abs_worst_dh"] = clean_value(row["abs_worst_dh"], 4)
        item["worst_delta"] = clean_value(row["worst_delta"], 4)
        rows.append(item)

    for row in end_b_candidates:
        item = _empty_row("end_b")
        item["id"] = row["candidate"]
        item["candidate"] = row["candidate"]
        item["worst_dh"] = clean_value(row["worst_dh"], 4)
        item["abs_worst_dh"] = clean_value(row["abs_worst_dh"], 4)
        item["worst_delta"] = clean_value(row["worst_delta"], 4)
        rows.append(item)

    for path in paths:
        item = _empty_row("cross_paths")
        item["id"] = f"{path['x']}-{path['y']}"
        item["x"] = path["x"]
        item["y"] = path["y"]
        item["worst_dh"] = clean_value(path["worst_dh"], 4)
        item["abs_worst_dh"] = clean_value(path["abs_worst_dh"], 4)
        item["worst_delta"] = clean_value(path["worst_delta"], 4)
        rows.append(item)

    for rank, row in enumerate(ranked_paths, start=1):
        item = _empty_row("cte_ranked")
        item["id"] = rank
        item["path"] = row["path"]
        item["x"] = row["x"]
        item["y"] = row["y"]
        item["rank"] = rank
        item["n_nodes"] = row["n_nodes"]
        alphas = row["alpha_values"] + [None] * (4 - len(row["alpha_values"]))
        item["alpha_0"] = clean_value(alphas[0], 3)
        item["alpha_1"] = clean_value(alphas[1], 3)
        item["alpha_2"] = clean_value(alphas[2], 3)
        item["alpha_3"] = clean_value(alphas[3], 3)
        item["elin"] = clean_value(row["elin"], 6)
        rows.append(item)

    return rows


FIELD_NAMES = [
    "section",
    "id",
    "path",
    "candidate",
    "x",
    "y",
    "worst_dh",
    "abs_worst_dh",
    "worst_delta",
    "rank",
    "n_nodes",
    "alpha_0",
    "alpha_1",
    "alpha_2",
    "alpha_3",
    "elin",
]


def write_results(rows):
    with open(RESULTS_PATH, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELD_NAMES)
        writer.writeheader()
        writer.writerows(rows)


def run(parameters_path=None):
    parameters = load_parameters(parameters_path)
    thresholds = parameters["thresholds"]
    end_a = parameters["end_a"]
    end_b = parameters["end_b"]
    candidates = candidate_symbols(parameters)

    end_a_scan = scan_end(end_a, candidates, thresholds)
    end_b_scan = scan_end(end_b, candidates, thresholds)
    end_a_candidates = [row for row in end_a_scan if row["passed"]]
    end_b_candidates = [row for row in end_b_scan if row["passed"]]
    paths = cross_paths(end_a_candidates, end_b_candidates, thresholds)

    cte = dict(parameters["cte"])
    cte[end_a["name"]] = end_a["cte"]
    cte[end_b["name"]] = end_b["cte"]
    used_candidates = sorted({path["x"] for path in paths} | {path["y"] for path in paths})
    missing = [element for element in used_candidates if cte.get(element) is None]
    if missing:
        raise SystemExit(
            "Missing CTE values for: "
            + ", ".join(missing)
            + "\nSupplement the `cte` object in `parameters.json`."
        )

    ranked_paths = cte_screen(
        paths,
        end_a["name"],
        end_b["name"],
        cte[end_a["name"]],
        cte[end_b["name"]],
        cte,
    )
    rows = build_results(end_a_candidates, end_b_candidates, paths, ranked_paths)
    write_results(rows)

    print(f"{end_a['name']} candidates: {len(end_a_candidates)}")
    print(f"{end_b['name']} candidates: {len(end_b_candidates)}")
    print(f"Cross paths: {len(paths)}")
    print(f"CTE-ranked paths: {len(ranked_paths)}")
    for row in ranked_paths:
        print(f"  {row['path']:<28} Elin={row['elin']:.6f}")
    print(f"Results written to: {RESULTS_PATH}")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else None)
