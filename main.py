# Build a Python program that:
# - Predicts molecular geometry using VSEPR
# - Computes bond dipole vectors from electronegativity differences
# - Determines overall molecular polarity
# - Visualizes dipoles and geometry


import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


main_group_elements = {
    # Group 1 (Alkali Metals + H) - 1 Valence Electron
    "H":  {"name": "Hydrogen",    "valence": 1, "electronegativity": 2.20},
    "Li": {"name": "Lithium",     "valence": 1, "electronegativity": 0.98},
    "Na": {"name": "Sodium",      "valence": 1, "electronegativity": 0.93},
    "K":  {"name": "Potassium",   "valence": 1, "electronegativity": 0.82},
    "Rb": {"name": "Rubidium",    "valence": 1, "electronegativity": 0.82},
    "Cs": {"name": "Cesium",      "valence": 1, "electronegativity": 0.79},
    "Fr": {"name": "Francium",    "valence": 1, "electronegativity": 0.70},

    # Group 2 (Alkaline Earth Metals) - 2 Valence Electrons
    "Be": {"name": "Beryllium",   "valence": 2, "electronegativity": 1.57},
    "Mg": {"name": "Magnesium",   "valence": 2, "electronegativity": 1.31},
    "Ca": {"name": "Calcium",     "valence": 2, "electronegativity": 1.00},
    "Sr": {"name": "Strontium",   "valence": 2, "electronegativity": 0.95},
    "Ba": {"name": "Barium",      "valence": 2, "electronegativity": 0.89},
    "Ra": {"name": "Radium",      "valence": 2, "electronegativity": 0.90},

    # Group 13 (Boron Group) - 3 Valence Electrons
    "B":  {"name": "Boron",       "valence": 3, "electronegativity": 2.04},
    "Al": {"name": "Aluminum",    "valence": 3, "electronegativity": 1.61},
    "Ga": {"name": "Gallium",     "valence": 3, "electronegativity": 1.81},
    "In": {"name": "Indium",      "valence": 3, "electronegativity": 1.78},
    "Tl": {"name": "Thallium",    "valence": 3, "electronegativity": 1.62},
    "Nh": {"name": "Nihonium",    "valence": 3, "electronegativity": None},

    # Group 14 (Carbon Group) - 4 Valence Electrons
    "C":  {"name": "Carbon",      "valence": 4, "electronegativity": 2.55},
    "Si": {"name": "Silicon",     "valence": 4, "electronegativity": 1.90},
    "Ge": {"name": "Germanium",   "valence": 4, "electronegativity": 2.01},
    "Sn": {"name": "Tin",         "valence": 4, "electronegativity": 1.96},
    "Pb": {"name": "Lead",        "valence": 4, "electronegativity": 2.33},
    "Fl": {"name": "Flerovium",   "valence": 4, "electronegativity": None},

    # Group 15 (Pnictogens) - 5 Valence Electrons
    "N":  {"name": "Nitrogen",    "valence": 5, "electronegativity": 3.04},
    "P":  {"name": "Phosphorus",  "valence": 5, "electronegativity": 2.19},
    "As": {"name": "Arsenic",     "valence": 5, "electronegativity": 2.18},
    "Sb": {"name": "Antimony",    "valence": 5, "electronegativity": 2.05},
    "Bi": {"name": "Bismuth",     "valence": 5, "electronegativity": 2.02},
    "Mc": {"name": "Moscovium",   "valence": 5, "electronegativity": None},

    # Group 16 (Chalcogens) - 6 Valence Electrons
    "O":  {"name": "Oxygen",      "valence": 6, "electronegativity": 3.44},
    "S":  {"name": "Sulfur",      "valence": 6, "electronegativity": 2.58},
    "Se": {"name": "Selenium",    "valence": 6, "electronegativity": 2.55},
    "Te": {"name": "Tellurium",   "valence": 6, "electronegativity": 2.10},
    "Po": {"name": "Polonium",    "valence": 6, "electronegativity": 2.00},
    "Lv": {"name": "Livermorium", "valence": 6, "electronegativity": None},

    # Group 17 (Halogens) - 7 Valence Electrons
    "F":  {"name": "Fluorine",    "valence": 7, "electronegativity": 3.98},
    "Cl": {"name": "Chlorine",    "valence": 7, "electronegativity": 3.16},
    "Br": {"name": "Bromine",     "valence": 7, "electronegativity": 2.96},
    "I":  {"name": "Iodine",      "valence": 7, "electronegativity": 2.66},
    "At": {"name": "Astatine",    "valence": 7, "electronegativity": 2.20},
    "Ts": {"name": "Tennessine",  "valence": 7, "electronegativity": None},

    # Group 18 (Noble Gases) - 8 Valence Electrons
    "He": {"name": "Helium",      "valence": 2, "electronegativity": None},
    "Ne": {"name": "Neon",        "valence": 8, "electronegativity": None},
    "Ar": {"name": "Argon",       "valence": 8, "electronegativity": None},
    "Kr": {"name": "Krypton",     "valence": 8, "electronegativity": 3.00},
    "Xe": {"name": "Xenon",       "valence": 8, "electronegativity": 2.60},
    "Rn": {"name": "Radon",       "valence": 8, "electronegativity": 2.20},
    "Og": {"name": "Oganesson",   "valence": 8, "electronegativity": None},
}






# AXE notation A reperesents central atom, X represents atoms bonded to A, and E represents lone pairs
vsepr_table = {
    # linear
    "AX2E0" : {"Electron Pair Geometry" : "Linear", "Bond Angle" : 180},

    # Trigonal planar
    "AX3E0" : {"Electron Pair Geometry" : "Trigonal Planar", "Bond Angle" : 120},
    "AX2E1" : {"Electron Pair Geometry" : "Bent", "Bond Angle" : 119},

    # Tetrahedral 
    "AX4E0" : {"Electron Pair Geometry" : "Tetrahedral", "Bond Angle" : 109.5},

    # based off of tetrahedral not so not accurate enough
    "AX3E1" : {"Electron Pair Geometry" : "Trigonal pyramidial", "Bond Angle" : 107.3},
    "AX2E2" : {"Electron Pair Geometry" : "Bent", "Bond Angle" : 104.5},
}








def get_electrons_and_elements(molecule):

    electrons = 0
    elements = {}

    i = 0
    # gets each element and adds up the number of electrons
    while i < len(molecule):
        if i+1 < len(molecule) and ((molecule[i] + molecule[i+1]) in main_group_elements):
            elements[molecule[i] + molecule[i+1]] = {"name": main_group_elements[molecule[i] + molecule[i+1]]["name"], "valence" : main_group_elements[molecule[i] + molecule[i+1]]["valence"], "electronegativity" : main_group_elements[molecule[i] + molecule[i+1]]["electronegativity"], "number" : 1}
            if i+2 < len(molecule):
                if molecule[i+2].isdigit():
                    electrons += int(molecule[i+2])*main_group_elements[molecule[i]+molecule[i+1]]["valence"]
                    elements[molecule[i] + molecule[i+1]] = {"name": main_group_elements[molecule[i] + molecule[i+1]]["name"], "valence" : main_group_elements[molecule[i] + molecule[i+1]]["valence"], "electronegativity" : main_group_elements[molecule[i] + molecule[i+1]]["electronegativity"], "number" : int(molecule[i+2])}
                    i += 3
                else:
                    electrons += main_group_elements[molecule[i]+molecule[i+1]]["valence"]
                    i += 2

            else:
                electrons += main_group_elements[molecule[i]+molecule[i+1]]["valence"]
                i += 2
            

        elif molecule[i] in main_group_elements:
            elements[molecule[i]] = {"name": main_group_elements[molecule[i]]["name"], "valence" : main_group_elements[molecule[i]]["valence"], "electronegativity" : main_group_elements[molecule[i]]["electronegativity"], "number" : 1}
            if i+1 < len(molecule):
                if molecule[i+1].isdigit():
                    electrons += int(molecule[i+1])*main_group_elements[molecule[i]]["valence"]
                    elements[molecule[i]] = {"name": main_group_elements[molecule[i]]["name"], "valence" : main_group_elements[molecule[i]]["valence"], "electronegativity" : main_group_elements[molecule[i]]["electronegativity"], "number" : int(molecule[i+1])}
                    i += 2
                else:
                    electrons += main_group_elements[molecule[i]]["valence"]
                    i += 1

            else:
                electrons += main_group_elements[molecule[i]]["valence"]
                i += 1

        else:
            break



    #checks the charge for the molecule and adjusts electron number accordingly
    if molecule[len(molecule) - 1] == "+":
        electrons -= 1

    elif molecule[len(molecule) - 1] == "-":
        electrons += 1

    if molecule[len(molecule) - 2] == "+":
        electrons -= int(molecule[len(molecule)-1])

    elif molecule[len(molecule) - 2] == "-":
        electrons += int(molecule[len(molecule)-1])


    return electrons, elements




def get_central_atom_and_surrounding_atoms(molecule):

    # gets the central atom and surrounding atoms by checking which atom has lowest electronegativity and if it is H. Set the central atom to most electronegative element in PT and compares using for loop except for the exception of HF for which this algorithm doesn't work.
    _, elements = get_electrons_and_elements(molecule)
    central_atom = "F"
    surrounding_atoms = {}

    if molecule[0] + molecule[1] == "HF":
        central_atom = "F"
        surrounding_atoms["H"] = {"name" : "Hydrogen", "valence" : 1, "electronegativity" : 2.20, "number" : 1}
    else:
        for i in elements:
            if elements[i]["electronegativity"] < main_group_elements[central_atom]["electronegativity"] and i != "H":
                central_atom = i
            else:
                surrounding_atoms[i] = {"name": elements[i]["name"], "valence": elements[i]["valence"], "electronegativity": elements[i]["electronegativity"], "number": elements[i]["number"]}
    
    return central_atom, surrounding_atoms



# gets number of bonded atoms by getting the number of each atom in surrounding_atoms and gets lone_pairs by taking electrons - 2*bonded atoms - 6*bonded atoms(if its not H) all divided by 2
def get_bonded_atoms_and_lone_pairs(molecule):

    electrons, _ = get_electrons_and_elements(molecule)
    _, surrounding_atoms = get_central_atom_and_surrounding_atoms(molecule)
    bonded_atoms = 0

    for i in surrounding_atoms:
        bonded_atoms += surrounding_atoms[i]["number"]


    lone_pairs = electrons-2*bonded_atoms
    for i in surrounding_atoms:
        if i != "H":
            lone_pairs -= 6*surrounding_atoms[i]["number"]


    lone_pairs /= 2
    lone_pairs = int(lone_pairs)

    return bonded_atoms, lone_pairs


# predicts the geometry using vsepr theory

def predict_geometry(molecule):
    bonded_atoms, lone_pairs = get_bonded_atoms_and_lone_pairs(molecule)
    key = f"AX{bonded_atoms}E{lone_pairs}"
    return vsepr_table[key]["Electron Pair Geometry"], vsepr_table[key]["Bond Angle"], key

def normalize(v):
    return v / np.linalg.norm(v)


def get_azimuthals(molecule):
    bonds, _ = get_bonded_atoms_and_lone_pairs(molecule)
    azimuthals = []

    i = 0
    if bonds == 1 or bonds == 2 or bonds == 3:
        while i < bonds:
            azimuthals.append((2*np.pi*i)/bonds)
            i += 1




    return azimuthals


def get_thetas(molecule):
    bonds, _ = get_bonded_atoms_and_lone_pairs(molecule)
    _, _, key = predict_geometry(molecule)
    bond_angle = (vsepr_table[key]["Bond Angle"]/180) * np.pi
    theta = []

    if bonds == 2:
        theta.append(bond_angle/2)

    elif bonds == 3:
        theta.append(np.arcsin(np.sqrt(-2/3*(np.cos(bond_angle)-1))))





    return theta


def get_vectors(molecule):
    azimuthals = get_azimuthals(molecule)
    theta = get_thetas(molecule)
    bonds, _ = get_bonded_atoms_and_lone_pairs(molecule)
    vectors = []
    if bonds == 1 or bonds == 2 or bonds == 3:
        for i in azimuthals:
            vectors.append(np.array([np.sin(theta[0])*np.cos(i), np.sin(theta[0])*np.sin(i), np.cos(theta[0])]))
    elif bonds == 4:
        vectors.append(np.array([0, 0, 1]))
        vectors.append(np.array([np.sin((109.5/180)*np.pi)*np.cos(0), np.sin((109.5/180)*np.pi)*np.sin(0), np.cos((109.5/180)*np.pi)]))
        vectors.append(np.array([np.sin((109.5/180)*np.pi)*np.cos((120/180)*np.pi), np.sin((109.5/180)*np.pi)*np.sin((120/180)*np.pi), np.cos((109.5/180)*np.pi)]))
        vectors.append(np.array([np.sin((109.5/180)*np.pi)*np.cos((240/180)*np.pi), np.sin((109.5/180)*np.pi)*np.sin((240/180)*np.pi), np.cos((109.5/180)*np.pi)]))




    

    return vectors

def get_bond_polarity_vectors(molecule):
    central_atom, surrounding_atoms = get_central_atom_and_surrounding_atoms(molecule)
    _, _, key = predict_geometry(molecule)
    directional_vectors = get_vectors(molecule)
    normalized_vectors = []
    bond_polarity_vectors = []

    for i in directional_vectors:
        normalized_vectors.append(normalize(i))



    i = 0
    for atom in surrounding_atoms:
        j = 0
        while j < surrounding_atoms[atom]["number"]:
            bond_polarity_vectors.append(normalized_vectors[j+i] * abs(surrounding_atoms[atom]["electronegativity"]-main_group_elements[central_atom]["electronegativity"]))
            j += 1
        i += j
    return bond_polarity_vectors



def get_molecular_polarity_vector(molecule):
    bond_polarity_vectors = get_bond_polarity_vectors(molecule)
    molecular_polarity_vector = np.array([0.0, 0.0, 0.0])

    for i in bond_polarity_vectors:
        molecular_polarity_vector += i

    return molecular_polarity_vector




polarity = ""

molecule = input("Molecule: ")

central_atom, surrounding_atoms = get_central_atom_and_surrounding_atoms(molecule)


molecular_geometry, bond_angles, _ = predict_geometry(molecule)


vectors = get_vectors(molecule)

molecular_polarity_vector = get_molecular_polarity_vector(molecule)


threshold = 1e-3

if np.linalg.norm(molecular_polarity_vector) < threshold:
    polarity = "Nonpolar"
else: 
    polarity = "Polar"




fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

for i in vectors:
    ax.quiver(0, 0, 0, i[0], i[1], i[2], color='b', arrow_length_ratio=0.1)



if polarity == "Polar":

    bond_dipole_pointing_out = False
    for atom in surrounding_atoms:
        if central_atom == "C" and main_group_elements[central_atom]["electronegativity"] < surrounding_atoms[atom]["electronegativity"]:
            ax.quiver(-1.5, 0, molecular_polarity_vector[2], molecular_polarity_vector[0], molecular_polarity_vector[1], molecular_polarity_vector[2])
            ax.text(-2, 0, 0.1, "net dipole", zdir=molecular_polarity_vector)
            bond_dipole_pointing_out = True
            break
    if bond_dipole_pointing_out == False:
        ax.quiver(-1.5, 0, molecular_polarity_vector[2], -molecular_polarity_vector[0], -molecular_polarity_vector[1], -molecular_polarity_vector[2])
        ax.text(-2, 0, 0.1, "net dipole", zdir=molecular_polarity_vector)




ax.text(-0.09,0,-0.25, central_atom)



i = 0

for atom in surrounding_atoms:
    j = 0
    while j < surrounding_atoms[atom]["number"]:
        ax.text(vectors[i][0], vectors[i][1], vectors[i][2]+0.1, atom)
        j += 1
        i += 1
    

# Set the axis limits for a clear view
ax.set_xlim([-2, 2])
ax.set_ylim([-2, 2])
ax.set_zlim([-2, 2])
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

# Display the plot
plt.title(f'{molecule}: {molecular_geometry}, {polarity}')
plt.show()



    




