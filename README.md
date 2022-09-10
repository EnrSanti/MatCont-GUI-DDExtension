# MatCont-GUI-DDExtension
The extension of MatCont user interface to systems of delay differential equations. Project done as Curricular Internship. 

CURRENT VERSION:

The current version of the project SHOULD allow the user to insert systems of DDEs and analyse their different aspects (i.e. do all the operations previously doable on ODE systems).

Currently mutliple DELAY PARAMETERs should be ACCEPTED, constants are also allowed as delays (e.g. x[t-1]).
Please note that when inserting the DDE use square brackets [] to specify the delay (e.g. y'=parameter\*y\*y[t-TAU]).
Delay equations must be written before renewal equations.

TODO:

-~~extend the system to accept multiple delay parameters~~

-~~extend the system to accept non constant delay parameters~~

-~~check wether in systems file use min/max~~

-~~continue code documentation~~

- eventual bug fixes & optimization

- code refactoring (i.e. the code may look a bit messy)
