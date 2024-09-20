| Revision # | 	1141 |
| ==========	 | ====|
|09/20/2024|	new keyword timeout. This is used for jobs taking more than a certain time fail, and known to fail.|
|06/24/2024|	Fixed a bug in IRC.sh and added keyword MaxBO to set the maximum bond order for breaking a bond |
|06/05/2024|      mopacamk changed to make it compatible with ASE's newest version. FutureWarnings from locate_barrless to barrless.err|
05/23/2024      Bugfixes in amk_parallel and sel_mol. locate_barrless modified to allow multiple searchs
03/06/2024	Removed dependence with zenity and gnuplot
01/15/2024	Fixed some bugs with llcalcs.sh and improved ChemKnow algorithm (new keyword CK_minima)
11/07/2023	Fixed some issues with the HL calcs of fragments
05/12/2022	Fixed some issues with qcore calcs
04/26/2022	Bugfix in utils and improved performace hl calcs
04/13/2022	Removed STOP in diag.f.
03/17/2022	Bugfix in final.sh.
03/11/2022	Interface with Gaussian 16.
02/03/2022	Min and Max temperatures for the kinetics set to 100 K and 9999 K, respectively.
12/10/2021	Bugfixes in Python scripts that read inputfile (charge and long mopac inputs were not read correctly).
12/01/2021	Bugfix in LocateTS.py.
11/18/2021	Bugfix in ChemKnow. Improved torsional search
11/16/2021	ChemKnow and barrierless search improved. Tutorial updated
11/10/2021	Tutorial updated and amk-tools described
10/17/2021	update tables with the geometry orientation for normal mode calcls in FINALDIR (low-level one). Tutorial updated
10/16/2021	Tutorial updated
10/02/2021	Tutorial updated
09/25/2021	References updated in tutorial and wiki
09/17/2021	tutorial updated and select.sh now exectuted from WRKDIR
09/06/2021	Max no. of species in population.pdf set to 20
07/26/2021	Bugfixes in FINAL and max no. of TSs in Energy Profile set to 100
07/24/2021	Bugfixes in screening
07/23/2021	Barrierless process are now included in the search (not in the kinetics)
06/25/2021	Screening and rxnetwork build optimized. Adaptive selection of temp in MD optimized.
06/18/2021	Bugfix in select_assoc script
06/17/2021	Systems with charge can now be modeled with ChemKnow and BXDE
04/06/2021	Bugfixes in irc_analysis.sh (mopac2021) and Heuristics (crossed bonds)
05/28/2021	New keyword Use_LET for mopac TS optimizations
05/27/2021	New keywords (recalc and Hookean) and improved efficiency
05/14/2021	ChemKnow improved and calcfc in g09 ts opt calcs
05/10/2021	keyword tight_ts and bug in bxde
05/02/2021	kinetics module improved and bugs corrected
04/07/2021	in tors.sh do not consider rotation about bonds that belong to a ring
03/16/2021	mopac calculator updated to ase3.21.1
03/14/2021	Bugfixes
02/15/2021	2021 release. Revision 1007--> qcore high-level calcs available
10/20/2020	linked_paths in python and 1D rotors in rrkm
10/18/2020	bug fix in linked_paths.sh
10/14/2020	tutorial updated and bug fixes in DVV.py
10/05/2020	Density matrix read in freq calc (MOPAC) and new xtb parameters for better convergence
07/15/2020	Implemented interface with Entos Qcore
04/21/2020	Implemented ExtForce and fixed some bugs
01/28/2020	Simplified adjacency matrix from XYZ
01/21/2020	2020 version
12/02/2019	Bug in FINAL.sh
11/22/2019	Maximum number of paths set to 50 (in bbfs.f90).
11/16/2019	Check of input structure in amk.sh. Fragmented molecules are no longer valid.
09/16/2019	pdfs are now also generated in FINAL_HL
07/09/2019	if name of working dir is too long, name--->wrkdir
06/30/2019	amk acronym replaces old tsscds acronym
04/18/2019	MIT license
04/17/2019	The label of the starting min in the kmc simulations is in tsdirll/KMC/starting_minimum
04/15/2019	A bug in get_energy_g09_MP2.sh was corrected
04/01/2019	threads=1 has been added to the input files in the examples folder. The use of this keyword is highly recommented to avoid multhreading in MOPAC calculations for much better performance
03/01/2019	A bug in the kmc.f90 source file was corrected
		
