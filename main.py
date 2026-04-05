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





# normalizes the vector making the magnitude 1
def normalize(v):
    return v / np.linalg.norm(v)


def get_repulsion_energy(vectors):
    repulsion_energy = 0
    for i in vectors:
        for j in vectors:
            eps = 1e-8
            r = vectors[i]["vector"]-vectors[j]["vector"]
            dist = np.linalg.norm(r)
            if dist < eps:
                continue

            if vectors[i]["type"] == "bp" and vectors[j]["type"] == "bp":
                weight = 1
            elif vectors[i]["type"] == "lp" and vectors[j]["type"] == "lp":
                weight = 1.5376
            else:
                weight = 1.24

            
            repulsion_energy += weight/dist

    return repulsion_energy


def get_force_direction(vectors, index):
    force_direction = np.array([0, 0, 0])
    vector = vectors[index]["vector"]

    for j in vectors:
        if j == index:
            continue
        eps = 1e-8
        r = vector - vectors[j]["vector"]
        dist = np.linalg.norm(r)
        if dist < eps:
            continue
        
        if vectors[index]["type"] == "bp" and vectors[j]["type"] == "bp":
            weight = 1
        elif vectors[index]["type"] == "lp" and vectors[j]["type"] == "lp":
            weight = 1.5376
        else:
            weight = 1.24

        force_direction = force_direction + weight*(r / (dist ** 3))

    return force_direction


def minimize_repulsion(molecule):
    bonded_atoms, lone_pairs = get_bonded_atoms_and_lone_pairs(molecule)
    electron_domains = bonded_atoms + lone_pairs
    vectors = {}
    rotation_size = 1e-4

    iterations = 0

    i=0

    while i < electron_domains:
        v = np.random.normal(0, 1, 3)
        v = normalize(v)



        if i < bonded_atoms:
            vectors[i] = {"type": "bp", "vector" : v}
        else:
            vectors[i] = {"type": "lp", "vector" : v}
        i += 1


    E = get_repulsion_energy(vectors)
    while True:
        for index in vectors:
            v = vectors[index]["vector"]
            force_direction = get_force_direction(vectors, index)
            F = force_direction - ((np.dot(force_direction, v)) * v)
            v_temp = v + (rotation_size * F)
            vectors[index]["vector"] = normalize(v_temp)
        
        E_prev = E
        E = get_repulsion_energy(vectors)


        
        forces_zeroed = 0
        for index in vectors:
            v = vectors[index]["vector"]
            force_direction = get_force_direction(vectors, index)
            f_tan = force_direction - ((np.dot(force_direction, v)) * v)
            if np.linalg.norm(f_tan) < 1e-6:
                forces_zeroed += 1

        if abs(E-E_prev) < 1e-8 and forces_zeroed == len(vectors):
            break

    return vectors



# Finds the magnitude of the bond vectors by subtracting the electronegativities and multiplies the normalized calculated vectors

def get_bond_dipoles(vectors):
    central_atom, surrounding_atoms = get_central_atom_and_surrounding_atoms(molecule)
    bond_dipoles = []





    i = 0
    for atom in surrounding_atoms:
        j = 0
        while j < surrounding_atoms[atom]["number"]:
            if vectors[i+j]["type"] == "bp":
                vector = normalize(vectors[i+j]["vector"])
                bond_dipoles.append(vector * abs(surrounding_atoms[atom]["electronegativity"]-main_group_elements[central_atom]["electronegativity"]))
                j += 1
        i += j
    return bond_dipoles



# adds up the bond vectors to determine the molecule vector

def get_molecular_dipole(bond_dipoles: list):
    molecular_dipole = np.array([0.0, 0.0, 0.0])
    for i in bond_dipoles:
        molecular_dipole += i
        

    return molecular_dipole

def draw_molecule(vectors):
    central_atom, surrounding_atoms = get_central_atom_and_surrounding_atoms(molecule)
    for i in vectors:
        if vectors[i]["type"] == "bp":
            ax.quiver(0, 0, 0, vectors[i]["vector"][0], vectors[i]["vector"][1], vectors[i]["vector"][2], color='b', arrow_length_ratio=0.1)

    ax.text(-0.09,0,-0.25, central_atom)

    i = 0

    for atom in surrounding_atoms:
        j = 0
        while j < surrounding_atoms[atom]["number"]:
            if vectors[i]["type"] == "bp":
                ax.text(vectors[i]["vector"][0], vectors[i]["vector"][1], vectors[i]["vector"][2]+0.1, atom)
                j += 1
            i += 1



def draw_dipole_moment(molecular_dipole):
    central_atom, surrounding_atoms = get_central_atom_and_surrounding_atoms(molecule)

    if polarity == "Polar":

        bond_dipole_pointing_out = False
        for atom in surrounding_atoms:
            if central_atom == "C" and main_group_elements[central_atom]["electronegativity"] < surrounding_atoms[atom]["electronegativity"]:
                ax.quiver(-1.5, 0, molecular_dipole[2], molecular_dipole[0], molecular_dipole[1], molecular_dipole[2])
                ax.text(-2, 0, 0.1, "net dipole", zdir=molecular_dipole)
                bond_dipole_pointing_out = True
                break
        if bond_dipole_pointing_out == False:
            ax.quiver(-1.5, 0, molecular_dipole[2], -molecular_dipole[0], -molecular_dipole[1], -molecular_dipole[2])
            ax.text(-2, 0, 0.1, "net dipole", zdir=molecular_dipole)




def get_bond_angles(vectors):
    for i in range(len(vectors)):
        for j in range(i+1, len(vectors)):
            if vectors[i]["type"] == "bp" and vectors[j]["type"] == "bp":
                angle = np.degrees(np.arccos(np.dot(vectors[i]["vector"], vectors[j]["vector"])))
                print(angle)
    


polarity = ""

molecule = input("Molecule: ")




#molecular_geometry, bond_angles, _ = predict_geometry(molecule)

#molecular_polarity_vector = get_molecular_polarity_vector(molecule)

vectors = minimize_repulsion(molecule)


get_bond_angles(vectors)

dipoles = get_bond_dipoles(vectors)


molecular_dipole = get_molecular_dipole(dipoles)

threshold = 1e-3

if np.linalg.norm(molecular_dipole) < threshold:
    polarity = "Nonpolar"
else: 
    polarity = "Polar"




fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')






draw_molecule(vectors)
draw_dipole_moment(molecular_dipole)





    

# Set the axis limits for a clear view
ax.set_xlim([-2, 2])
ax.set_ylim([-2, 2])
ax.set_zlim([-2, 2])
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

# Display the plot
plt.title(f'{molecule}: {polarity}')
plt.show()



    




