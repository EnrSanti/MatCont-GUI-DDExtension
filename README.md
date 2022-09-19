# MatCont-GUI-DDExtension
The extension of MatCont user interface allows users to insert and study systems containing delay differential equations.
Project done as Curricular Internship. 

CURRENT VERSION:

The current version of the project should allow the user to insert systems of DDEs and analyse their different aspects (i.e. do all the operations previously doable on ODE systems) through the GUI.

Currently mutliple DELAY PARAMETERs should be accepted, constants are also allowed as delays (e.g. x[t-1]).
Please note that when inserting the DDE use square brackets [] to specify the delay (e.g. y'=parameter\*y\*y[t-TAU]).
Renewal eqautions are also supported (note that delay equations must be written before renewal equations).

TODO:

-~~extend the system to accept multiple delay parameters~~

-~~extend the system to accept non constant delay parameters~~

-~~check wether in systems file use min/max~~

-~~continue code documentation~~

- eventual bug fixes & optimization

- code refactoring (i.e. the code may look a bit messy)
