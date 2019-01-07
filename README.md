#### SpliceHetero: An information theoretic model for measuring spliceomic intratumor heterogeneity from bulk-tumor RNA-seq #####
** Version 1.0. [2019-01-03]


1. Requirements
- Python 2.7
- numpy
- argparse
- scipy


2. Usage: python SpliceHetero.py [arguments]

Essemtial Arguments:

'-csp', '--Case_Sample_path' [PATH]             List of case-sample input (see CANCER_SAMPLES.txt for input format)

'-rsp', '--Reference_Sample_path' [PATH]                List of reference-sample input (see NORMAL_SAMPLES.txt for input format)

'-odp', '--Out_Dir_path' [PATH]         Directory path for output'

Optional Arguments:

'-prn', '--Process_number' [INT]                Number of processes to use (default: 1)

'-slb', '--Stranded_Library_bool' [BOOL]                If it is stranded library (options: True/False, default: False)'


#### Please see run.test.sh for more detailed usage

3. Contact: mdy89@snu.ac.kr


4. Web-page: http://biohealth.snu.ac.kr/software/SpliceHetero



