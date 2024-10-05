# Set the trace
set -uex

# SRR number
SRR=SRR519926

# Number of reads to sample
N=10000

# The output read names
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq

# Trimmed read names
T1=reads/${SRR}_1.trimmed.fastq
T2=reads/${SRR}_2.trimmed.fastq

# The adapter sequence
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

# The reads directory
RDIR=reads

# The reports directory
PDIR=reports

# Micromamba environment name
ENV_NAME=menv

# ----- actions below ----

# Initialize micromamba shell for bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate ${ENV_NAME}

# Ensure fastp is installed
if ! command -v fastp &> /dev/null
then
    echo "fastp not found, installing..."
    micromamba install -y -n ${ENV_NAME} fastp
fi

# Make the necessary directories
mkdir -p ${RDIR} ${PDIR}

# Download the FASTQ file
fastq-dump -X ${N} --split-files -O ${RDIR} ${SRR}

# Run fastqc
fastqc -q -o ${PDIR} ${R1} ${R2}

# Run fastp and trim for quality
fastp --adapter_sequence=${ADAPTER} --cut_tail \
      -i ${R1} -I ${R2} -o ${T1} -O ${T2}

# Run fastqc on the trimmed files
fastqc -q -o ${PDIR} ${T1} ${T2}

# Run multiqc on the reports directory
micromamba run -n ${ENV_NAME} multiqc -o ${PDIR} ${PDIR}
