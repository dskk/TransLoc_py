# TransLoc_py

This is a reimplementation of Ciuchi's [TransLoc](https://github.com/CiuK1469/TransLoc) (accompanying the paper "Practical computation of the charge mobility in molecular semiconductors using transient localization theory") in Python with additional features.

This code outputs square of the transient localization length (L2) along $x$ and $y$ axes at each eigenenergies of the system.


Requirements
---

- Python 3.6+ (The code works fine in the author's Python 3.9.7 environment)
- `numpy`, `scipy` and `toml` libraries


Usage
---

- Make your (empty) working directory and copy `transloc.py` into it.
- Prepare input files listed below and locate them in the same directory.
- Run `transloc.py`.
Output files will be generated in a subdirectory.

Sample
---

You can copy `transloc.py` and `RTA_L2.pyd` to a sample folder and run `transloc.py`. This will reproduce our results. If there is an error in importing `RTA_L2.pyd`, you can delete it (`transloc.py` automatically switches to the *slow* mode without it).

When executing a trial run, you can set smaller values to the entries `nrepeat`, `supercell_a` and `supercell_b` in `system.toml` to shorten execution time.

Input files
---

|Filename|Mandatory / Optional|Description|
|:---|:---:|:---|
|system.toml|Mandatory|Describes the system (2D lattice) to calculate|
|settings.toml|Optional|Controls output options|
|RTA_L2.pyd|Optional<br>(Highly recommended to use)|Faster implementation of the `RTA_L2` function in C++|

Detailed information of each file is written in the following sections.


system.toml
---

This file describes the system (2D lattice) to calculate. Entries are listed below. Refer to the sample file for how to write.

**You can use any unit of length and energy, but they must all be the same throughout the file.** The units in the output file will be the same.

<img src="https://github.com/dskk/TransLoc_py/blob/71ff555332ad5ac886e56998b47ee3e38a5c4469/readme_image.png" width="50%" />

This image is a visualization of `pentacene_sample/system.toml`. The horizontal axis is the $x$ axis and the vertical axis is the $y$ axis. $a$ axis is parallel to the $x$ axis and $b$ axis is slightly tilted from $y$ axis. The orange parallelogram is the unit cell containing the molecules m1 and m2. Interactions are indicated by arrows. Replicated cells are numbered in ascending order from left to right and from bottom to top, respectively. 

### system_name
Name of the system. This value will not be used for the calculation. Use as a memo.

### energy_range_min, energy_range_max, energy_sep
The lower and upper limits and increments of the energy for outputting the calculation results. All eigenenergies must fall within this energy range, so it is recommended to set a range bit wider than the bandwidth (That is, setting $\pm \textrm{bandwidth}/2$ as the upper and lower limit, respectively. Note that the median of eigenenergies is expected to be zero). See also the sample file and "Output files" section.

### inv_tau
$\hbar/\tau$ where $\tau$ is the representative length of time that each molecule stays in place (=typical intermolecular oscillation period).

### gauss_broadening_width
FWHM of a Gaussian function for smoothing to convert the descrete result into a continuous function. Specify a minute value like 5 meV. See the code for detail.

### nrepeat
Number of times the calculation will be repeated.

### cell_vec_a, cell_vec_b
Cell vectors along $a$ and $b$ axes. Use Cartesian coordinates.

### supercell_a, supercell_b
Number of times the cell will be replicated along the $a$ and $b$ axes.

### molecules
Definition of molecules in the unit cell. Label uniquely like m1, m2 etc. and define Cartesian coordinates (see the sample file).

### interactions
Definition of interactions between molecules. Specify labels of interacting molecules and interaction type. If the second molecule (mol2) is not in the same cell as the first molecule (mol1), specify the (non-zero) difference in the cell number along $a$ and $b$ axes.

### transfer_integrals
Labels for interaction types (I1, I2... in the sample file) and corresponding binding energies. If a binding energy has some randomness, you can make a list with length greater than 1 and have the code randomly select it upon Hamiltonian construction. This random selection will be made independently for each of the Hamiltonian matrix elements corresponding to the interaction type.

settings.toml
---

This file controls the following output options. Refer to the sample file for how to write.
|Key name|Allowed values|Notes|
|:---|:---|:---|
|output_folder_name|Any strings|Name of a subdirectory to be created to store output files|
|dump|"yes" / "no"|"yes" to enable dump output|
|overwrite|"yes" / "no"|Effective only when a directory with the same name as output directory already exists.<br>"yes" to force overwrite. "no" to force abort. Otherwise, a confirmation message will be prompted.|

If `settings.toml` does not contain any of the keys above (or you do not have `settings.toml` in your working directory), you will be asked to set options.


RTA_L2.pyd
---

This is a binary file containing Python-callable function generated by [pybind11](https://github.com/pybind/pybind11).

Follow the tutorial in the [documentation](https://pybind11.readthedocs.io/en/latest/) and compile `RTA_L2.cpp` to generate `RTA_L2.pyd`.

The `RTA_L2.pyd` file included in this project (generated on a ordinary laptop with an 11th gen Intel Core i7 Processor) may also work on your computer. You do not need to compile `RTA_L2.cpp` in this case.

Without this file, `transloc.py` automatically call `RTA_L2` funcion implemented in Python. It works, but takes dozens of times longer to process.


Output files
------------

### {iteration_number}/x(y).txt
Contains square of the transient localization length (L2).
The L2 value in each line is for the energy level written in the corresponding line of `energy.txt`.
The dimension of L2 value is the square of length and its unit is the same as the unit of length used in `system.toml`.
Multiply by $e/2\tau k_B T$ to get the mobility.

This is the Y-data of L2 vs energy plot.

### energy.txt
Contains an arithmetic sequence defined by `energy_range_min`, `energy_range_max` and `energy_sep` in `system.toml`.

This is the X-data of L2 vs energy plot.

### {iteration_number}/dump.txt
Contains values in the following format.

```math
\begin{array}{ll}
N &\textrm{Number of molecules}\\
x_1,~y_1 &\textrm{Position of molecule }1\\
\vdots & \vdots \\
x_N,~y_N &\textrm{Position of molecule }N\\
H_{1,1} \ldots H_{1,N} &\textrm{Hamiltonian first row}\\
\vdots & \vdots \\
H_{N,1} \ldots H_{N,N} &\textrm{Hamiltonian last row}\\
E_1 &\textrm{Lowest eigenenergy}\\
\vdots & \vdots \\
E_N &\textrm{Highest eigenenergy}\\
V_{1,1} \ldots V_{1,N} &\textrm{Eigenvector corresponding to } E_1\\
\vdots & \vdots \\
V_{N,1} \ldots V_{N,N} &\textrm{Eigenvector corresponding to } E_N
\end{array}
```
