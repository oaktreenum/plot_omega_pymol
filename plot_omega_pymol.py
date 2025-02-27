import pymol
from pymol import cmd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from scipy.special import i0
import seaborn as sns
import numpy as np
from typing import Tuple, Optional

"""
Scripts for plotting torsion angles for proteins and peptoids, either selections or multiple objects

Allon Goldberg
Research Assistant, Flatiron Institute (NYC), Center for Computational Biology, Biomolecular Design Group
2/25/2025

USAGE:
1. If not installed your PyMOL environment, (pip) install matplotlib, seaborn, scipy, numpy         (in PyMOL terminal)
2. run <path_to_script>/plot_omega_pymol.py             (in PyMOL terminal)
3. Plot! 

COMMANDS:
    help(<function_name>)                             ——— print detailed specific function info
    plot_omega_info                                   ——— print info on functions
    get_omega_angles SELECTION_STRING, PRINT_BOOL     ——— retrieve selection's omega angles and optionally print values to terminal
    get_phipsi_angles SELECTION_STRING, PRINT_BOOL    ——— retrieve selection's phi/psi angles and optionally print values to terminal
    plot_omega SELECTION_STRING, PRINT_BOOL           ——— plot selection's omega angles and optionally print values to terminal
    plot_rama SELECTION_STRING, PRINT_BOOL            ——— plot selection's Ramachandran plot (psi v phi angles) and optionally print values to terminal
    plot_all_angles SELECTION_STRING, PRINT_BOOL      ——— plot all selection's torsion angles (omega, phi, psi) and optionally print values to terminal
"""



def get_omega_angles(selection: str = 'All', print_true: Optional[bool] = 0) -> Tuple[list,list,list,list]:
    """
    Grab omega (ω) angles for PyMOL selection and OPTIONAL print
    Omega is the dihedral angle between CA(i-1), C(i-1), N(i), CA(i)

    Parameters:
        selection (str): The selection string (e.g., 'chain A', 'resi 6-13')
        print_bool (int): A boolean-like value to control printing. Use 1/True to print or 0/False (default).

    Example: 
        get_omega_angles chain A or resi 6-13, 1
        get_omega_angles("chain A or resi 6-13", 1)

    For More Info:
        plot_omega_info
    """
    object_list = cmd.get_object_list(selection)
    omega_angles = []
    residue_names = []
    residue_index_list = []  

    for obj in object_list:
        model = cmd.get_model(f"{obj} and {selection}")
        residue_names_temp = []
        residue_index_temp = []
        cmd.iterate(f"{obj} and name CA", "residue_names_temp.append(resn)", space=locals())
        cmd.iterate(f"{obj} and {selection} and name CA", "residue_index_temp.append(int(resi))", space=locals())
        # residue_names_temp = [atom.resn for atom in model.atom if atom.name == "CA"]
        # residue_index_temp = [int(atom.resi) for atom in model.atom if atom.name == "CA"]
        residue_names.append(residue_names_temp)
        omega_temp = []
        residue_index_used = []
        for i in residue_index_temp:
            if i>1:
                omega = cmd.get_dihedral(
                    f"{obj} and resi {i-1} and name CA",
                    f"{obj} and resi {i-1} and name C",
                    f"{obj} and resi {i} and name N",
                    f"{obj} and resi {i} and name CA",
                )
                omega_temp.append(omega)
                residue_index_used.append(i)
        residue_index_list.append(residue_index_used)
        omega_angles.append(omega_temp)
    
    # Print if print_true is 1 or True
    if print_true:
        for i in range(len(object_list)):
            for j in range(len(omega_angles[i])):
                print(f'{object_list[i]} [{residue_names[i][j+1]}]{residue_index_list[i][j]}: Pre-OMG = {omega_angles[i][j]:.3f}')
    return object_list, residue_names, omega_angles, residue_index_list

def get_phipsi_angles(selection: str = 'All', print_true: Optional[bool] = 0) -> Tuple[list,list,list,list]:
    """
    Grab Phi (φ) and Psi (ψ) angles for PyMOL selection and optionally print them to terminal in readable fashion
    Phi is the dihedral angle between C(i-1), N(i), CA(i), and C(i)
    Psi is the dihedral angle between N(i), CA(i), and C(i), N(i+1)
    
    This utilizes the built-in cmd.get_phipsi() command

    Parameters:
        selection (str): The selection string (e.g., 'chain A', 'resi 6-13')
        print_bool (int): A boolean-like value to control printing. Use 1/True to print or 0/False (default).

    Example: 
        get_phipsi_angles chain A or resi 6-13, 1
        get_phipsi_angles("chain A or resi 6-13", 1)

    For More Info:
        plot_omega_info
    """
    object_list = cmd.get_object_list(selection)
    phi_angles = []
    psi_angles = []
    colors = []
    resnames = []
    resis = []

    angles_obj = cmd.get_phipsi(selection)
    for (object_name, atom_num), (phi, psi) in angles_obj.items():
        phi_angles.append(phi)
        psi_angles.append(psi)
        cmd.iterate(f"{object_name} and index {atom_num}","colors.append(color); resis.append(resi); resnames.append(resn)", space=locals())
        # resnames.append({var_resn.get('r', 'N/A')})
        # Print if print_true is 1 or True
        if print_true:
            print(f"{object_name} [{resnames[-1]}]{resis[-1]}: Phi = {phi:.3f}, Psi = {psi:.3f}")
    
    for i in range(0,len(colors)):
        colors[i] = cmd.get_color_tuple(colors[i])

    return phi_angles, psi_angles, colors, resnames

def plot_omega(selection: str ='All', print_true: Optional[bool] = 0):
    """
    Plot omega (ω) angles for a given selection with matplotlib

    Parameters:
        selection (str): The selection string (e.g., 'chain A', 'resi 6-13')
        print_bool (int): A boolean-like value to control printing. Use 1/True to print or 0/False (default).

    Example: 
        plot_omega chain A or resi 6-13, 1
        plot_omega("chain A or resi 6-13", 1)

    For More Info:
        plot_omega_info
    """

    if print_true:
            object_name_list, residue_names_array, omega_angles_array, angles_indeces= get_omega_angles(selection, 1)
    else:
        object_name_list, residue_names_array, omega_angles_array, angles_indeces= get_omega_angles(selection)

    np_res_names = np.array(residue_names_array).flatten()
    np_omega_angles = np.array(omega_angles_array).flatten()
    np_omega_indeces = np.array(angles_indeces).flatten()

    # Custom cmap
    colors_for_cmap = ['royalblue', 'orange',  'tomato',  'orange', 'springgreen',  'orange', 'tomato', 'orange', 'royalblue'] 
    cyclic_cmap = mcolors.LinearSegmentedColormap.from_list('my_cyclic', colors_for_cmap, N=180)
    # Setup figure
    plt.figure(figsize=(4, 4))
    plt.scatter(np_omega_indeces, np_omega_angles, marker='o', linestyle='-', c=np_omega_angles, cmap=cyclic_cmap)
    
    ### Change y to highlight angle, i.e. ±180, 0 for cis/trans
    plt.axhline(y=180, color='g', linestyle='--')
    plt.axhline(y=-180, color='g', linestyle='--')
    plt.axhline(y=0, color='g', linestyle='--')

    ### Change x range to highlight residue ranges, i.e. loops
    # plt.axvspan(10.5, 16.5, color='grey', alpha=0.13, label='Loops')
    # plt.axvspan(23.5, 29.5, color='grey', alpha=0.13)
    # plt.axvspan(36.5, 42.5, color='grey', alpha=0.13)      

    plt.xlabel('Residue Number')
    plt.ylabel('Omega Angle (°)')
    plt.xticks(np.arange(0, np_omega_indeces[-1] + 1, 5))   
    plt.yticks(range(-180, 185, 30))     
    plt.title(f'Omega (ω) Angles for {selection}')
    plt.ylim(-180, 180)
    # plt.legend()

    manager = plt.get_current_fig_manager()
    manager.set_window_title(f"omega_plot_{selection}")
    plt.grid()
    plt.show(block=False)
    
    return

def vonmises_kde(data, kappa, n_bins=180):
    """
    Adjust data for polar hist plot 
    * May not be working
    """
    theta = np.linspace(-np.pi, np.pi, n_bins)
    kde = np.exp(kappa * np.cos(theta[:, None] - data[None, :])).sum(1) / (2 * np.pi * i0(kappa))
    return theta, kde

def plot_rama(selection: str ='All', print_true: Optional[bool] = 0):
    """
    Plot Ramachandran plots for a given selection with matplotlib
    POINTS COLORED BY THE COLOR OF THE RESIDUE IN THE PyMOL SESSION

    Parameters:
        selection (str): The selection string (e.g., 'chain A', 'resi 6-13')
        print_bool (int): A boolean-like value to control printing. Use 1/True to print or 0/False (default).

    Example: 
        plot_rama chain A or resi 6-13, 1
        plot_rama("chain A or resi 6-13", 1)

    For More Info:
        plot_omega_info
    """
    
    phi, psi, colors, resnames = get_phipsi_angles(selection, print_true)

    fig, ax = plt.subplots(figsize=(4.3, 4))

    ax.scatter(phi, psi, marker='o', linestyle='-', c=colors, alpha=0.888, linewidths=.3, edgecolors='grey')
    ax.set_xlabel('Phi φ (°)')
    ax.set_ylabel('Psi ψ (°)')
    ax.set_xticks(range(-180, 185, 45))   
    ax.set_yticks(range(-180, 185, 45))     
    ax.set_title(f'Rama Plot for {selection}')
    ax.set_xlim(-180, 180) 
    ax.set_ylim(-180, 180) 
    ax.grid()
    ax.axvline(x=0, color='black', linestyle='-')
    ax.axhline(y=0, color='black', linestyle='-')
    
    manager = plt.get_current_fig_manager()
    manager.set_window_title(f"rama_plot_{selection}")
    plt.tight_layout()
    plt.show(block=False)

    return

def plot_all_angles(selection: str ='All', print_true: Optional[bool] = 0):
    """
    Plot all angle plots for a given selection with matplotlib inlcuding histogram plots for each angle

    Parameters:
        selection (str): The selection string (e.g., 'chain A', 'resi 6-13')
        print_bool (int): A boolean-like value to control printing. Use 1/True to print or 0/False (default).

    Example: 
        plot_all_angles chain A or resi 6-13, 1
        plot_all_angles("chain A or resi 6-13", 1)

    For More Info:
        plot_omega_info
    """
    object_name_list, residue_names_array, omega_angles_array, angles_indeces= get_omega_angles(selection, print_true)
    np_res_names = np.array(residue_names_array).flatten()
    np_omega_angles = np.array(omega_angles_array).flatten()
    np_omega_indeces = np.array(angles_indeces).flatten()
    phi, psi, colors, resnames = get_phipsi_angles(selection, print_true)

    # Custom cmap
    colors_for_cmap = ['royalblue', 'orange',  'tomato',  'orange', 'springgreen',  'orange', 'tomato', 'orange', 'royalblue'] 
    cyclic_cmap = mcolors.LinearSegmentedColormap.from_list('my_cyclic', colors_for_cmap, N=180)
    # Setup figure
    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(2, 6, height_ratios=[1, 3])
    ax0 = fig.add_subplot(gs[0, 0:2], projection=None)
    ax1 = fig.add_subplot(gs[0, 2:4], projection=None)
    ax2 = fig.add_subplot(gs[0, 4:6], projection=None)
    ax3 = fig.add_subplot(gs[1, 0:3])
    ax4 = fig.add_subplot(gs[1, 3:6])

    # For polar hist plot *May not be working correctly
    omg_theta, omg_kde =  vonmises_kde(np_omega_angles, 20)
    phi_theta, phi_kde =  vonmises_kde(np.array(phi), 20)
    psi_theta, psi_kde =  vonmises_kde(np.array(psi), 20)

    # Omega kde histogram
    sns.kdeplot(np_omega_angles, fill=True, color="darkviolet", bw_adjust=0.1, ax=ax0)
    # For polar hist plot (below)
    # ax0.plot(omg_theta, omg_kde, color="darkviolet", lw=1) 
    ax0.set_title(f'Omega (ω) Hist', y=1)
    ax0.tick_params(axis='x', labelsize=8)
    ax0.tick_params(axis='y', labelsize=0)
    ax0.set_xticks(range(-180,183,90))
    ax0.set_xlim(-180,180)
    ax0.set_yticks([])
    ax0.set_ylabel('')
    
    # Phi kde histogram
    sns.kdeplot(phi, fill=True, color="deepskyblue", bw_adjust=.23, ax=ax1)
    # For polar hist plot (below)
    # ax1.plot(phi_theta, phi_kde, color="deepskyblue", lw=1)
    ax1.set_title(f'Phi (φ) Hist', y=1)
    ax1.tick_params(axis='x', labelsize=8)
    ax1.tick_params(axis='y', labelsize=0)
    ax1.set_xticks(range(-180,183,90))
    ax1.set_xlim(-180,180)
    ax1.set_yticks([])
    ax1.set_ylabel('')

    # Psi kde histogram
    sns.kdeplot(psi, fill=True, color="cadetblue", bw_adjust=.23, ax=ax2)
    # For polar hist plot (below)
    # ax2.plot(psi_theta, psi_kde, color="cadetblue", lw=1)
    ax2.set_title(f'Psi (ψ) Hist', y=1)
    ax2.tick_params(axis='x', labelsize=8)
    ax2.tick_params(axis='y', labelsize=0)
    ax2.set_xticks(range(-180,180,90))
    ax2.set_xlim(-180,180)
    ax2.set_yticks([])
    ax2.set_ylabel('')

    # Omega plot
    ax3.scatter(np_omega_indeces, np_omega_angles, marker='o', linestyle='-', c=np_omega_angles, cmap=cyclic_cmap)
    ax3.axhline(y=180, color='g', linestyle='--')
    ax3.axhline(y=-180, color='g', linestyle='--')
    ax3.axhline(y=0, color='g', linestyle='--')
    ax3.set_xlabel('Residue Number')
    ax3.set_ylabel('Omega ω (°)')
    ax3.set_xticks(np.arange(0, np_omega_indeces[-1] + 1, 5))   
    ax3.set_yticks(range(-180, 185, 30))     
    ax3.set_title(f'Omega (ω) Angles for {selection}')
    ax3.set_ylim(-180, 180) 
    ax3.grid()
    
    # Ramachandran plot
    ax4.scatter(phi, psi, marker='o', linestyle='-', c=colors, alpha=0.888, linewidths=.3, edgecolors='grey')
    ax4.set_xlabel('Phi φ (°)')
    ax4.set_ylabel('Psi ψ (°)')
    ax4.set_xticks(range(-180, 185, 45))   
    ax4.set_yticks(range(-180, 185, 45))     
    ax4.set_title(f'Rama Plot for {selection}')
    ax4.set_xlim(-180, 180) 
    ax4.set_ylim(-180, 180) 
    ax4.grid()
    ax4.axvline(x=0, color='black', linestyle='-')
    ax4.axhline(y=0, color='black', linestyle='-')

    manager = plt.get_current_fig_manager()
    manager.set_window_title(f"all_torsion_plot_{selection}")
    plt.tight_layout()
    plt.show(block=False)
    
    return

def plot_omega_info():

    print("""
    Scripts for plotting torsion angles for proteins and peptoids, either selections or multiple objects

    Allon Goldberg
    Research Assistant, Flatiron Institute (NYC), Center for Computational Biology, Biomolecular Design Group
    2/25/2025

    USAGE:
    1. If not installed your PyMOL environment, (pip) install matplotlib, seaborn, scipy, numpy         (in PyMOL terminal)
    2. run <path_to_script>/plot_omega_pymol.py             (in PyMOL terminal)
    3. Plot! 

    COMMANDS:
        help(<function_name>)                             ——— print detailed specific function info
        plot_omega_info                                   ——— print info on functions
        get_omega_angles SELECTION_STRING, PRINT_BOOL     ——— retrieve selection's omega angles and optionally print values to terminal
        get_phipsi_angles SELECTION_STRING, PRINT_BOOL    ——— retrieve selection's phi/psi angles and optionally print values to terminal
        plot_omega SELECTION_STRING, PRINT_BOOL           ——— plot selection's omega angles and optionally print values to terminal
        plot_rama SELECTION_STRING, PRINT_BOOL            ——— plot selection's Ramachandran plot (psi v phi angles) and optionally print values to terminal
        plot_all_angles SELECTION_STRING, PRINT_BOOL      ——— plot all selection's torsion angles (omega, phi, psi) and optionally print values to terminal
    """)

    return



# Extend commands
# (essentially creates a shorthand version of the command for the PyMOL terminal)
cmd.extend("plot_omega_info", plot_omega_info)
cmd.extend("get_omega_angles", get_omega_angles)
cmd.extend("get_phipsi_angles", get_phipsi_angles)
cmd.extend("plot_omega", plot_omega)
cmd.extend("plot_rama", plot_rama)
cmd.extend("plot_all_angles", plot_all_angles)

# Enable selection autocompletion for extensions
cmd.auto_arg[0]['get_omega_angles'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['get_phipsi_angles'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['plot_omega'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['plot_rama'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['plot_all_angles'] = cmd.auto_arg[0]['zoom']