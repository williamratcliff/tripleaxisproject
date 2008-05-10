@echo off
echo **-------------------------------------------------------**
echo **----                                               ----**
echo **---- CRYSTALLOGRAPHIC FORTRAN MODULES LIBRARY 3.00 ----**
echo **----                                               ----**
echo **---- Authors: JRC-JGP                   (1999-2005)----**
echo **-------------------------------------------------------**
rem
if x%1 == x goto FIRST
if x%1 == xlibc goto LIBC
rem
echo Special Module for F compiler
g95 -c CFML_mod.f95              -O
rem
:FIRST
echo Level 0: Mathematical Modules
g95 -c CFML_math_gen.f95         -O
g95 -c CFML_random.f95           -O
g95 -c CFML_fft.f95              -O
g95 -c CFML_fft_nF.f95           -O
g95 -c CFML_Profile_TOF.f95      -O
g95 -c CFML_Profile_Finger.f95   -O
g95 -c CFML_Profile_Functs.f95   -O
rem
echo Level 0: Graphical Modules
g95 -c CFML_Isosurface.f95       -O
echo Level 0: Optimization Modules
g95 -c CFML_optimization.f95     -O
g95 -c CFML_optimization_lsq.f95 -O
echo End Compilation Level 0
rem
echo Level 1:  Module Applications
g95 -c CFML_math_3D.f95          -O
g95 -c CFML_spher_harm.f95       -O
g95 -c CFML_string_util.f95      -O
echo End Compilation Level 1
rem
echo Level 2:  Module Applications
g95 -c CFML_sym_table.f95        -O0
g95 -c CFML_chem_scatt.f95       -O0
g95 -c CFML_cryst_types.f95      -O
g95 -c CFML_diffpatt.f95         -O
echo End Compilation Level 2
rem
echo Level 3:  Module Applications
g95 -c CFML_bonds_table.f95      -O0
g95 -c CFML_symmetry.f95         -O
echo End Compilation Level 3
rem
echo Level 4:  Module Applications
g95 -c CFML_reflct_util.f95      -O
g95 -c CFML_atom_mod.f95         -O
echo End Compilation Level 4
rem
echo Level 5:  Module Applications
g95 -c CFML_sfac.f95             -O
g95 -c CFML_propagk.f95          -O
g95 -c CFML_geom_calc.f95        -O
g95 -c CFML_maps.f95             -O
g95 -c CFML_molecules.f95        -O
echo End Compilation Level 5
rem
echo Level 6:  Module Applications
g95 -c CFML_form_cif.f95         -O
g95 -c CFML_conf_calc.f95        -O
g95 -c CFML_refcodes.f95         -O
echo End Compilation Level 6
rem
echo Making CrysFML Library
g95 -c CFML_io_mess.f95          -O
rem
echo Level 7: Optimization Simulated Annealing
g95 -c CFML_optimization_san.f95 -O
echo End Compilation Level 7
rem
ar cr libgcrysfml.a *.o
copy *.mod .\LibG95\.
move *.a   .\LibG95\.
del *.o
del *.mod
