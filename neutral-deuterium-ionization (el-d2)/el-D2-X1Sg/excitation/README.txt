MCCC calculations of electron scattering on molecular hydrogen and its isotopologues
Adiabatic nuclei calculations performed with the spheroidal MCCC-S(210) model.
Reference: Scarlett et al, Atom. Data Nucl. Data Tables (2021).

- Files are organised in directories according to the initial vibrational level
  in the ground electronic (X1Sg) state.

- File names are in the format:
    MCCC-el-D2-[f]_[final vib. state].[i]_[initial vib. state].txt
  with [f]/[i] = final/initial electronic state.
  E.g. MCCC-D2-B1Su_vf=5.X1Sg_vi=0.txt

- In addition to the fully vibrationally resolved transitions, there are files named as follows:
    MCCC-el-D2-[f]_bound.X1Sg_vi=0.txt     (sum over final bound vibrational levels)
    MCCC-el-D2-[f]_DE.X1Sg_vi=0.txt        (dissociative excitation)
    MCCC-el-D2-[f]_total.X1Sg_vi=0.txt     (sum over all final vibrational levels)

- Fitting parameters are provided in similar files in the fits directory

- All electronic states have a number of bound vibrational levels, except for the
  b3Su state which is purely dissociative.

- List of vibrational level energies in each electronic state are provided in the
  vibrational_energies directory

- Electronic states are labelled according to their molecular symmetries:
  
  Symmetry    Description      States
  ------------------------------------------------
  1Sg         singlet (s=0)    X1Sg, EF1Sg, GK1Sg, H1Sg
              Sigma   (M=0)
              gerade

  1Su         singlet (s=0)    B1Su, Bp1Su (B'1Su)
              Sigma   (M=0)
              ungerade

  1Pu         singlet (s=0)    C1Pu, D1Pu
              Pi      (|M|=1)
              ungerade

  1Pg         singlet (s=0)    I1Pg
              Pi      (|M|=1)
              gerade
  
  1Dg         singlet (s=0)    J1Dg
              Delta   (|M|=2)
              gerade

  3Sg         triplet (s=0)    a3Sg, g3Sg, h3Sg
              Sigma   (M=0)
              gerade

  3Su         triplet (s=0)    b3Su, e3Su
              Sigma   (M=0)
              ungerade

  3Pu         triplet (s=0)    c3Pu, d3Pu
              Pi      (|M|=1)
              ungerade

  3Pg         triplet (s=0)    i3Pg
              Pi      (|M|=1)
              gerade

  3Dg         triplet (s=0)    j3Dg
              Delta   (|M|=2)
              gerade
