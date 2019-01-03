#### SpliceHetero: An information theoretic model for measuring spliceomic intratumor heterogeneity from bulk-tumor RNA-seq #####
** Version 1.0. [2019-01-03]


1. Requirements
- Python 2.7
- numpy
- argparse
- scipy


2. Usage: python SpliceHetero.py [arguments]

Essemtial Arguments:
'-csp', '--Case_Sample_path' [PATH]             List of case-sample input (see web-page for more input format)
'-rsp', '--Reference_Sample_path' [PATH]                List of reference-sample input (see web-page for more input format)
'-odp', '--Out_Dir_path' [PATH]         Directory path for output'

Optional Arguments:
'-prn', '--Process_number' [INT]                Number of processes to use (default: 1)
'-slb', '--Stranded_Library_bool' [BOOL]                If it is stranded library (options: True/False, default: False)'


3. Test-data: please type following command for testing
sh run.test.sh

