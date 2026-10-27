C_Mathematica_derivations  (Supplementary File 6)
==================================================
Wolfram Mathematica 15 notebooks of the analytical results of the manuscript,
each with a PDF of the evaluated notebook (all outputs and figures).
Open a notebook in Mathematica (desktop or Wolfram Cloud) and evaluate all
cells (Evaluation > Evaluate Notebook). No external files are needed.

Analytical_Derivations_StepByStep_final.nb   (evaluated: ..._final.pdf)
  Section 1  species, complexes and reaction directions (Methods, Reaction network)
  Section 2  conservation laws and deficiency                  main-text Eqs (3)-(6)
  Section 3  reduced mass-action model and analytic Jacobian   Eqs (1), (2)
  Section 4  compact form of the Jacobian (support/jac6.m)
  Section 5  exact positive steady-state parametrisation u=K/P Eqs (7)-(9)
  Section 6  symbolic checks of the flux and total balances    Eqs (7)-(9)
  Section 7  detached-branch and branch-point conditions       Eqs (10), (11)
  Section 8  slope of the crossing branch; Class 7/8 criterion Eqs (13), (14)
  Section 9  fold and branch-point normal-form coefficients    Methods

Excel_Parameter_Examples_Condition_Checks_final.nb   (evaluated: ..._final.pdf)
  One ensemble parameter set of every exact class (sampling iteration and row
  in Ensemble_iter<i>_v3_1e6/Excel_Files/ML_Ready_Dataset_v3_iter<i>.csv):
  Michaelis constants, both inequalities of criterion (10) (Section 3) and the
  crossing-slope formula (13)/(14) for Classes 6-8 (Section 4).

Figure cells
  Every section of both notebooks ends with a text cell "Added figure. ..." and a
  cell that plots the result of that section numerically. The set-up cell for the
  figures is at the end of Section 1 and must be evaluated before the other figure
  cells (Evaluate Notebook does this). All figure code uses names beginning with
  "fig", so it does not change any symbolic result. Each figure is also saved as a
  PDF in a sub-folder "plots" next to the notebook. The figure cells carry the cell
  tag "AddedPlotCell".
  Analytical notebook:
    1 reaction network; 2 stoichiometric matrix and conservation in the full model;
    3 two stable states at A_TOT = 23.21 nM, Jacobian eigenvalues, dose-up and
      dose-down sweeps (hysteresis through the detached branch);
    4 compact Jacobian against finite differences, stability of the three states;
    5 A(u), A_TOT(u) and the complete diagram of the main-text Fig 3 set;
    6 residuals of the parametrisation; 7 criterion (10) and the branch-point
      condition; 8 pitchfork vs transcritical, Eq. (13) slope; 9 normal forms.
  Examples notebook:
    1 parameters of the ten examples; 2 exact diagram of every example;
    3 criterion (10) plane and Michaelis constants; 4 branch points and slopes of
      Classes 6-8; 5 summary of what the parametrisation finds for each example.
