This script fits Kd values to chemical-shift perturbation data.

Usage:
run fit_Kd_paul.py -d WT_shifts_ratio -conc 50e-6 -n 1 -Kd 1e-6  -o output -Dmax option_fixedDmax.in

It can fix the maximum chemical-shift perturbation to a value specified in the inputfile option_fixedDmax.in.
It can also use a global fixed Kd value for everyone.
Right now it cannot do a global fit of a single Kd to all residues.

It also creates a file with the fitted curve (for plotting eg in xmgrace)
