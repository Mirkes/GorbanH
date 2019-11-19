# Gorban's H family of Lyapunove functions: software and figures
Software and figures for paper  "Universal Gorban’s Entropies: Geometric Case Study" in the Entropy journal
This repository contains four folders for four considered models. Each model contains several submodels with different figure names.
Description of file names is presented below.

Presented code was not specially optimised. Calculation of partial equilibria and level sets can be significantly accelerated on base of subsection "Partial equilibria for several typical reactions" of paper.

## Linear model. Folder name is "Linear Isomerisation"
There are two submodels: with detailed ballance and with complex balance but without detained balance. Six different set of parameters is considered for each submodel.

Main file to start is ModelFigures.m
EqSel.mat contains three selected equilibria. Script "ModelFigures.m" allows you to generate another equilibria if necessary.

All figures are started from "Mod2" prefix.
There are following types of figures:
1. DirectNNN is partial equilibria lines and projections from point c for equilibrium witn number NNN
1. LevelsNNN is partial equilibria lines and level sets of Gorban's H for equilibrium witn number NNN
1. LevelsBHNNN is partial equilibria lines and level sets of Boltzmann's H for equilibrium witn number NNN
1. TrajectN_M is partial equilibria lines and trajectory for equilibrium N and set of reaction rate constants M for system with complex balance
1. TrajectDetN_M is partial equilibria lines and trajectory for equilibrium N and set of reaction rate constants M for system with detailed balance
1. TrajectHN_M presents graphs of partial and total Gorban's H for full reaction time equilibrium N and set of reaction rate constants M for system with complex balance
1. TrajectDetHN_M presents graphs of partial and total Gorban's H for full reaction time equilibrium N and set of reaction rate constants M for system with detailed balance
1. TrajectHTN_M presents graphs of partial and total Gorban's H for initial time interval equilibrium N and set of reaction rate constants M for system with complex balance
1. TrajectDetHTN_M presents graphs of partial and total Gorban's H for initial time interval  equilibrium N and set of reaction rate constants M for system with detailed balance
1. TrajectHTBN_M presents graphs of Boltzmann's H and Gorban's H for initial time interval equilibrium N and set of reaction rate constants M for system with complex balance
1. TrajectDetHTBN_M presents graphs of Boltzmann's H and Gorban's H for initial time interval  equilibrium N and set of reaction rate constants M for system with detailed balance

## Nonlinear isomerisation model. Folder name is "Nonlinear Isomerisation"
This model does not have submodels. There are figures for three equilibria and four reaction reate constants sets.

Main file to start is ModelFigures.m
EqSel.mat contains three selected equilibria. Script "ModelFigures.m" allows you to generate another equilibria if necessary.

There are following types of figures:
1. directNNN is partial equilibria lines and projections from point c for equilibrium witn number NNN
1. levelsNNN is partial equilibria lines and level sets of Gorban's H for equilibrium witn number NNN
1. levelsBHNNN is partial equilibria lines and level sets of Boltzmann's H for equilibrium witn number NNN
1. trajectN_M is partial equilibria lines and trajectory for equilibrium N and set of reaction rate constants M for system with complex balance
1. trajectHN_M presents graphs of partial and total Gorban's H for full reaction time equilibrium N and set of reaction rate constants M for system with complex balance
1. trajectHTN_M presents graphs of partial and total Gorban's H for initial time interval equilibrium N and set of reaction rate constants M for system with complex balance
1. TrajectHTBN_M presents graphs of Boltzmann's H and Gorban's H for initial time interval equilibrium N and set of reaction rate constants M for system with complex balance

## Hydrogen Chloride production.  Folder name is "HCl"
There are two submodels with abstract parameters (prefix HClMod) and with real Hydrogen Chloride production parameters (prefix HCl).

Main file to start is HClFigures.m
mat files contain results of system integrations for different cases:
ModXTim is full time solution for artificial model X.
ModXTimStart is solution for artificial model X for short initial fragment.
ModXTimStartStart is solution for artificial model X for very short initial fragment.
timX are seven files with solution of real HCl model for full time (each file contains one partial time interval)
timStart is file with solution of real HCl model for initial time interval.
To recalculate one or more files it is necessary to remove them from work directory and restart HClPrepare for real HCl System or HClModPrepare for abstract model. Abstract model must be specified in "modified" variable in HClModPrepare for calculations and in HClFigure to form figures. Value modified=0 in HClFigure corresponds to real system.

There are following types of figures:
1. BLevelXXXPosN is level set with Boltzmann's H function value XXX (for example -0.5) with  position of camera N (0 or 1).
1. GLevelXXXPosN is level set with Gorban's H function value XXX (for example -0.5) with  position of camera N (0 or 1).
1. Conc present graphs of concentrations for full reaction time
1. ConcStart present graphs of concentrations for initial fragment of reaction time
1. H presents graphs of partial and total Gorban's H for full reaction time
1. HStart presents graphs of partial and total Gorban's H for initial time interval
1. HGHB presents graphs of Boltzmann's H and Gorban's H for full reaction time
1. HGHBStart presents graphs of Boltzmann's H and Gorban's H for initial time interval
1. LevelSections presentes horizontal sections of level set for Gorban's H=-0.9
1. Partial_N presents partial equilibria surface for reaction N
1. Polyhedron presents the reaction polyhedron.
1. Trajectory presents trajectory in the reaction polyhedron.

## Water Gas Shiift reaction.  Folder name is "WGS"
There are two submodels: abstract (prefix WGSMod) and real (prefix WGS) WGS reactions. To select abstract reaction specify modified = 1 in WGSFigures and modified = 0 for real WGS reaction.

Main file to start is WGSFigures.m
File WGSParamFound contains results of Nelder mead search of parameters for real WGS system. To restart search of parameters it is necessary to remove this file and run WGSFigures.m.
All other mat files contain results of system integrations for different cases:
WGSRes (WGSModRes) contains solution for WGS (abstract WGS) model for 0.59sec (3.2sec).
WGSResTrajectory (WGSModResTrajectory) contains solution for WGS (abstract WGS) model for 1.77sec (9.6sec).
WGSResShortTime (WGSModResShortTime) contains solution for WGS (abstract WGS) model for 5.9microseconds (1.6sec).

There are following types of figures:
1. BLevels and GLevels present level sets for Boltzmann's and Gorban's H correspondingly.
1. Concentrations present graphs of all consentrations as function of time.
1. ConcentrationsLOg present graphs of all consentrations as function of time in the logarithmic time scale.
1. H presents graphs of partial and total Gorban's H for full reaction time
1. HStart presents graphs of partial and total Gorban's H for initial time interval
1. HGHB presents graphs of Boltzmann's H and Gorban's H for full reaction time
1. HGHBStart presents graphs of Boltzmann's H and Gorban's H for initial time interval
1. Partial presents two partial equilibria lines and trajectory of solution in the reaction polygon.
1. Partial presents trajectory of solution in the reaction polygon.


