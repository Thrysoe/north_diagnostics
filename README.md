# north_diagnostics
NORTH diagnostics reader and processing script

This program was designed to run with python 3 in a Spyder 6 environnement, with usual packages plus nptdms, cv2 and local libraries 
(probe.py, diagnostics.py, file_reader.py and dau.py) found in the folder diagnostics and utils of the repository. To access them,
it is necessary to include their location in the PYTHONPATH with the sys library and sys.path.append('The:\\location\\of\\mylibraries')
or add it permanently in the PYTHONPATH (see the settings of the interpretor).

To install modules on Spyder 6, the command pip install my_library can encounter difficulties. A way to avoid that issue is to 
install the libraries with the pip command on the python application and change the python interpretor of Spyder following:
Tools > Preferences > Python interpretor and then give the location of the python.exe file.

Each executable program are at the root of the repository and where tested in a Spyder 6 environment. A short description of each program
is given at its beginning, along with the input parameters. Only the imput parameters of each program should be modified.
