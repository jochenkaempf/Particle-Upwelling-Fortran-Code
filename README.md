# Particle-Upwelling-Fortran-Code
Modified Versions of the COHERENS model (Luyten et al., 1999) uses to simulate particle upwelling on continental slopes. The original model was written in Fortran 77 and has been converted to Fortran 90 by this author.

In order to run this model, one must create a folder "dat" for the data output. When using the G95 Fortran compiler, the model code can be compiled and run with the instructions:

g95 -c param.f90 Functions.f90 cohini.f90 cohrun.f90

g95 -o runB.exe bossNEW.f90 param.o Functions.o cohini.o cohrun.o

runB.exe 

The uploaded PDF entitled "OriginalSubmissionSmall" includes descriptions of all model settings implemented in the model code. Note that this code refers to the experiment "B-narrow" that uses the bathymetry file "topo1.dat" as input file. The FORTRAN code "canyon1.f95" created this bathymetry. It can also be used to produce the other bathymetries describes in the scientific text.   

The reference for the CORERENS model is: Luyten, P. J., Jones, J. E., Proctor, R., Tabor, A., Tett, P., & Wild-Allen, K. (1999), COHERENS - A Coupled Hydrodynamical-Ecological Model for Regional and Shelf Seas: User Documentation; MUMM Report; Management Unit of the North Sea: Brussels, Belgium, 914p. A PDF of the user documentation is available at: https://www.icbm.de/fileadmin/user_upload/icbm/ag/physoz/download/from_emil/COHERENS/print/userguide.pdf
