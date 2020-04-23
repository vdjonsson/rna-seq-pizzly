import pandas as pd

### read sample sheet ###
samples = pd.read_csv(config["samples"], sep="\t", dtype=str)
# set sample as index of dataframe
samples = samples.set_index("sample", drop=False)
samples.index.names = ["sample_id"]

workdir: config['workdir']
report: "report/workflow.rst"

### Set constants for references ###
REFERENCE = config['reference']['name']
REFPATH = config['refdir'] + REFERENCE + '/'
GTF = REFPATH + config['reference']['gtf']
FASTA = REFPATH + config['reference']['fasta']
INDEX = REFPATH + "{0}.idx".format(config['reference']['prefix'])
T2G = REFPATH + "t2g.pkl"
PROJECT =  config['project']

### helper functions ###
def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    u = samples.loc[wildcards.sample, :].dropna()
    return [ f"{u.fqdir}"+f"{u.fq1}", f"{u.fqdir}"+f"{u.fq2}" ]
