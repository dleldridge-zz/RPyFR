**Author:** Daniel L. Eldridge, Ph.D. (daniel.lee.eldridge@gmail.com)

Spring 2020

#RPyFR Script
**A script for wrangling vibrational frequencies from Gaussian output files and computing RPFRs
over user-specified temperature ranges.**

**Author:** Daniel L. Eldridge, Ph.D. (daniel.lee.eldridge@gmail.com)

For this script to work, you must follow the following file naming format for freqchk output files (.txt):

    MoleculeID_Theory_BasisSet_IsotopicMoleculeID_note.txt

**Examples**:
1. *Methane*_CCSD_augccpVTZ_0nosub_atom2.txt
2. *Methane*_CCSD_augccpVTZ_12CH3D_atom2.txt
3. *Sulfate40H2O*_B3LYP_631Gdp_33SO4_atom1.txt
4. *Thiosulfate34H2O*_B3LYP_631Gdp_32S32S(16O)3_none.txt

The unsubstituted molecule must have the IsotopicMoleculeID of "0nosub" in order for the RPFRs to be
computed correctly.

Tips:
Avoid special characters for basis set (and other) descriptions and use shorthand whenever possible. Also, only place files pertaining to a single optimization/frequency calculation in the data_files folder at one time for processing.

Try it out:
The script will work with the example output files involving methane isotopologues that are included. The specific calculations for methane isotopologues were recently published in the following article:

Lloyd, Eldridge, and Stolper (2021) *Geochim. Cosmochim. Acta* (https://doi.org/10.1016/j.gca.2020.10.008).
