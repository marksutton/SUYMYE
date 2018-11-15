
This package distributed under GPL v3, Copyright (C) 2017 Mark Sutton 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
                         
You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.

1. Installation
MBL2016 requires the Qt library (www.qt.io). It has been tested with Qt versions 5.3 and 5.5 and 5.8, but should be compatible with all Qt version greater than 5.3. It has been tested only Windows and MacOS, but will in theory work on any Qt-capable system. To compile and run,  (a) install Qt, (b) open the .pro file MBL2016.pro in Qt Creator, (c) specify build kit, (d) ensure a randomnumbers.dat file exists in the build directory. This should contain a sequence of random bytes, and should be at least 65536 bytes in length. Suitable data can be downloaded from https://qrng.anu.edu.au/


MBL2016 uses QCustomPlot (http://www.qcustomplot.com/) v 1.3.1 for visualisation - source code included in this package, see comments for license


2. Usage

To initiate an MBL run, change parameters in GUI then use Start from the file menu. A simulation can be stopped prematurely using Stop. Use 'Chart to PDF' from the file menu to generate a PDF version of the onscreen chart.


Main parameters (always visible in main window)

Extinction probability: the probability of a lineage going extinct in a time increment (0=never, 1=always)

Speciation probability: the probability of a lineage speciating in a time increment (0=never, 1=always). NOTE - the program does not check that these sum to 1 or less.

FDT cut-off: How many increments old the most recent common ancestor of two species has to be for them to be considered to be different genera under the FDT algorithm

Time increments per tree - how many increments will each tree be run for.

Maximum leaf-count - when a tree reaches this number of leaves the tree-simulation will stop for that tree. If the program is in 'throw away trees hitting leaf limit' mode (Mode menu - recommended) these trees will be discarded. Otherwise they will be included in taxonomy with the 'present' considered to be the timestep in which they halted. To avoid any such behavior, set this to a very large value (100000 or more). Note that trees hitting the limit are reported in the log window even if the program is not in 'output results per tree' mode

Trees to run - how many separate MBL trees to simulate

Filename stub - if tree files are output (see below), the program prepends any text in this box to their name


Tree Export Menu

Tick any of the three export items to have MBL2016 generate a tree file with the selected taxonomy for each tree it successfully generates. The tree format is selected in this menu as well. More than one tree type can be exported at once, but only a single format can be selected. Tree export will only take place if an output folder has been chosen with the final item on this menu

Mode menu

Contains checkable menu-items to control the behaviour of the program.

Include RTDT/FDT results table - if checked these items generate tables of genus-count frequencies suitable for export (as csv files) via copy/paste to other software

Log scale - if checked the x (genus size) scale of the graph is logarithmic

Throw away trees hitting leaf limit - see maximum leaf count above

Continue until correct number of trees. In this mode, the program will continue generating trees until it has successfully generated the required number, rather than simply trying the specified number of generation runs. Trees that have 0 extant leaves at the end of the number of increments specified OR that hit the leaf limit will hence be re-tried. Be aware that with poorly chosen parameters, a simulation set up in this way may never finish, though it can always be halted with 'stop' from the file menu.

Output results per tree - if this is checked the program operates in verbose mode, logging the size of each tree found (and more besides).


Please contact m.sutton@ic.ac.uk for further details.
