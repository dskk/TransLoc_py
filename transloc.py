import os,shutil,sys
import toml
import numpy as np
from scipy.integrate import simps
from scipy.interpolate import splrep, splev
if os.path.isfile("RTA_L2.pyd"):  # You can run this script without RTA_L2.pyd (but very slow)
    try:
        import RTA_L2 as cpp_RTA_L2
    except:
        print("An error occurred while importing RTA_L2.pyd")
        sys.exit(1)
else:
    cpp_RTA_L2 = None


class Connection:  # This donetes a connection between two molecules
    def __init__(self, delta_x, delta_y, energy, index1, index2):
        self.delta_x = delta_x  # Difference in X-coordinates
        self.delta_y = delta_y  # Difference in Y-coordinates
        self.energy = energy  # Coupling energy
        self.index1 = index1  # Index of molecule 1
        self.index2 = index2  # Index of molecule 2


def constract_H(connections, nmol):  # This constracts the Hamiltonian using connection information
    H = np.zeros((nmol, nmol))
    for con in connections:
        i, j = con.index1, con.index2
        H[i][j] = con.energy
        H[j][i] = con.energy
    return H


# A reimplementation of Ciuchi's original RTA_L2 function. It is highly recommended to use faster code in C++ by importing RTA_L2.pyd.
def RTA_L2(eig_vals, eig_vecs, connections, inv_tau):
    ncon = len(connections)
    nmol = len(eig_vals)
    x = np.zeros(nmol)
    y = np.zeros(nmol)
    for i in range(nmol):
        for j in range(nmol):
            temp_x = 0.0
            temp_y = 0.0
            for con in connections:
                v1 = eig_vecs[con.index1]
                v2 = eig_vecs[con.index2]
                temp_x += (v1[i]*v2[j]-v1[j]*v2[i]) * con.delta_x * con.energy
                temp_y += (v1[i]*v2[j]-v1[j]*v2[i]) * con.delta_y * con.energy
            x[i] += (temp_x**2) * 2 / (inv_tau**2+(eig_vals[i]-eig_vals[j])**2)
            y[i] += (temp_y**2) * 2 / (inv_tau**2+(eig_vals[i]-eig_vals[j])**2)
    return x, y


# Another implementation of the RTA_L2 function (this returns the same result). Provided as an aid to understanding.
def RTA_L2_modified(eig_vals, eig_vecs, connections, inv_tau):
    ncon = len(connections)
    nmol = len(eig_vals)
    x = np.zeros(nmol)
    y = np.zeros(nmol)
    HX_x = np.zeros((nmol, nmol))
    HX_y = np.zeros((nmol, nmol))
    l = []
    for con in connections:
        HX_x[con.index1][con.index2] = con.delta_x * con.energy
        HX_y[con.index1][con.index2] = con.delta_y * con.energy
        l.append((con.index1, con.index2))
    for i in range(nmol):
        for j in range(nmol):
            su_x = 0
            su_y = 0
            v1 = eig_vecs.T[i]
            v2 = eig_vecs.T[j]
            mat = np.outer(v1, v2)
            mat -= mat.T
            for tup in l:
                i1, i2 = tup
                su_x += mat[i1][i2] * HX_x[i1][i2]
                su_y += mat[i1][i2] * HX_y[i1][i2]
            x[i] += (su_x**2) * 2 / (inv_tau**2+(eig_vals[i]-eig_vals[j])**2)
            y[i] += (su_y**2) * 2 / (inv_tau**2+(eig_vals[i]-eig_vals[j])**2)
    return x, y


# Initial Y-data is y_in, X-data is x_in. This function resamples (x_resample will be used as new X-data) and smoothens Y-data.
def gauss_broad(x_resample, x_in, y_in, fwhm):
    integ_in = simps(y_in, x_in)
    iota = [i for i in range(len(x_in))]
    f = splrep(iota, x_in, k=5, s=3)
    waight_arr = np.gradient(splev(iota, f))
    # Broaden y_in
    yy = np.zeros(len(x_resample))
    for i in range(len(x_in)):
        x = x_in[i]
        y = y_in[i] * waight_arr[i]
        for j in range(len(x_resample)):
            d = x - x_resample[j]
            s = fwhm**2 / 4 / np.log(2)  # 2*sigma**2
            yy[j] += y * np.exp(-d**2/s)
    # Keep integrated value the same
    scale_factor = integ_in/simps(yy, x_resample)
    return yy * scale_factor


assert os.path.isfile("system.toml"), "system.toml not found"
if os.path.isfile("settings.toml"):
    settings_toml = toml.load("settings.toml")
else:
    settings_toml = dict()

system_toml = toml.load("system.toml")
inv_tau = system_toml["inv_tau"]
gauss_broadening_width = system_toml["gauss_broadening_width"]
nrepeat = system_toml["nrepeat"]
energy_range_min = system_toml["energy_range_min"]
energy_range_max = system_toml["energy_range_max"]
energy_sep = system_toml["energy_sep"]
cell_vec_a = np.array(system_toml["cell_vec_a"])
cell_vec_b = np.array(system_toml["cell_vec_b"])
na = system_toml["supercell_a"]
nb = system_toml["supercell_b"]
molecules = system_toml["molecules"]
nmol_cell = len(molecules.keys())
interactions = system_toml["interactions"]
transfer_integrals = system_toml["transfer_integrals"]

# Load and parse settings.toml
output_folder_name = None  # Output folder name
if "output_folder_name" in settings_toml.keys():
    output_folder_name = settings_toml["output_folder_name"]
if not output_folder_name:
    print("Enter output folder name:", end="")
    output_folder_name = input()
dump_flag = None  # Whether to output a dump file or not
if "dump" in settings_toml.keys():
    dump_flag = settings_toml["dump"].lower() in ["y", "yes"]
if dump_flag is None:
    print("Dump Hamiltonian, eigenvalues, eigenvectors etc?(y/N):", end="")
    dump_flag = input().lower() in ["y", "yes"]
if os.path.isdir(output_folder_name):
    if "overwrite" in settings_toml.keys() and settings_toml["overwrite"].lower() in ["y", "yes"]:
        shutil.rmtree(output_folder_name)
    elif "overwrite" in settings_toml.keys() and settings_toml["overwrite"].lower() in ["n", "no"]:
        print(f"Folder '{output_folder_name}' already exists. Aborting.")
        sys.exit()
    else:
        print("Overwrite?(y/N):", end="")
        if input().lower() in ["y", "yes"]:
            shutil.rmtree(output_folder_name)
        else:
            print(f"Aborting.")
            sys.exit()

os.mkdir(output_folder_name)
os.mkdir(f"{output_folder_name}/input_files")
shutil.copy2("system.toml", f"{output_folder_name}/input_files/system.toml")
if os.path.isfile("settings.toml"):
    shutil.copy2("settings.toml",
                 f"{output_folder_name}/input_files/settings.toml")

resample_energy = np.linspace(energy_range_min, energy_range_max, energy_sep)
np.savetxt(f"{output_folder_name}/energy.txt", resample_energy)

molecule_position_dic = dict()
molecule_name_to_index_dic = dict()
molecule_index_to_position_arr = []
i = 0
for k in molecules.keys():
    pos = np.array([molecules[k]["x"], molecules[k]["y"]])
    molecule_position_dic[k] = pos
    molecule_index_to_position_arr.append(pos)
    molecule_name_to_index_dic[k] = i
    i += 1

# For dump output
molecule_positions_arr = [None for i in range(na*nb*nmol_cell)]
for i in range(nb):
    for j in range(na):
        for k in range(nmol_cell):
            molecule_positions_arr[i*na*nmol_cell+j*nmol_cell+k] = \
                cell_vec_a * j + cell_vec_b * i + \
                molecule_index_to_position_arr[k]


for rep_cnt in range(nrepeat):
    connections = []
    for i in range(nb):  # list up all connections
        for j in range(na):
            for k in interactions.keys():
                interaction = interactions[k]
                molname1 = interaction["mol1"]
                molname2 = interaction["mol2"]
                ind1 = (i*na+j) * nmol_cell+molecule_name_to_index_dic[molname1]  # Index of molecule1
                ii = (i+interaction["cell_translation_b"]) % nb  # i for molecule2. Periodic boundary is considered.
                jj = (j+interaction["cell_translation_a"]) % na  # j for molecule2. Periodic boundary is considered.
                ind2 = (ii*na+jj) * nmol_cell+molecule_name_to_index_dic[molname2]  # Index of molecule2
                transfer_integral = np.random.choice(transfer_integrals[interaction["type"]])  # Coupling energy
                cell_origin = cell_vec_a * j + cell_vec_b * i
                r1 = cell_origin + molecule_position_dic[molname1]  # Position of molecule1
                r2 = cell_origin + molecule_position_dic[molname2] + \
                     cell_vec_a * interaction["cell_translation_a"] + \
                     cell_vec_b * interaction["cell_translation_b"]  # Position of molecule2
                delta_x, delta_y = r2 - r1
                connections.append(Connection(delta_x, delta_y, transfer_integral, ind1, ind2))
    H = constract_H(connections, na*nb*nmol_cell)
    # this function calculates eigenvalues of symmetric real matrix
    eig_vals, eig_vecs = np.linalg.eigh(H)
    nmol = len(eig_vals)
    if cpp_RTA_L2 is not None:  # If RTA_L2.pyd is provided
        # The list "connections" cannot be passed to a C++ function.
        # Following lines unpack it and do some pre-calculations to make C++ vectors (HX_x and HX_y).
        HX_x = np.zeros((nmol, nmol))
        HX_y = np.zeros((nmol, nmol))
        inds = []
        for con in connections:
            HX_x[con.index1][con.index2] = con.delta_x * con.energy
            HX_y[con.index1][con.index2] = con.delta_y * con.energy
            inds.append([con.index1, con.index2])
        x, y = cpp_RTA_L2.RTA_L2(eig_vals, eig_vecs, HX_x, HX_y, inds, inv_tau)
    else:  # You can run this script without RTA_L2.pyd (but very slow)
        x, y = RTA_L2(eig_vals, eig_vecs, connections, inv_tau)
    py_broad_x = gauss_broad(resample_energy, eig_vals, x,
                             gauss_broadening_width)
    py_broad_y = gauss_broad(resample_energy, eig_vals, y,
                             gauss_broadening_width)
    os.mkdir(f"{output_folder_name}/{rep_cnt+1:0>4d}")
    np.savetxt(f"{output_folder_name}/{rep_cnt+1:0>4d}/x.txt", py_broad_x)
    np.savetxt(f"{output_folder_name}/{rep_cnt+1:0>4d}/y.txt", py_broad_y)
    if dump_flag:
        with open(f"{output_folder_name}/{rep_cnt+1:0>4d}/dump.txt", "w") as fo:
            fo.write(f"{nmol}\n")
            for i in range(nmol):
                fo.write(
                    f"{molecule_positions_arr[i][0]} {molecule_positions_arr[i][1]}\n")
            for i in range(nmol):
                for j in range(nmol-1):
                    fo.write(f"{H[i][j]} ")
                fo.write(f"{H[i][-1]}\n")
            for i in range(nmol):
                fo.write(f"{eig_vals[i]}\n")
            for i in range(nmol):
                for j in range(nmol-1):
                    fo.write(f"{eig_vecs[i][j]} ")
                fo.write(f"{eig_vecs[i][-1]}\n")
