## A command for calculating sITH of samples in CANCER_SAMPLES.txt against the samples in NORMAL_SAMPLES.txt
## "-prn 3" indicates the program will use three processes in parallel
## "-slb False" indicates the RNA-seq library of input is not strand-specific (if you don't know which, it is safe to choose False)
## The input format is BED. See CANCER_SAMPLES.txt and NORMAL_SAMPLES.txt for input format

python SpliceHetero.py -csp CANCER_SAMPLES.txt -rsp NORMAL_SAMPLES.txt -odp TEST_OUT.Dir -prn 3 -slb False
