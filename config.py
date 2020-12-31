import glob
INDIR="/mnt/chr4/mikrobiomy-2/Wyniki_sekwencjonowania/demultiplexed"
#OUTDIR="/mnt/chr7/data/julia"

OUTDIR="/home/julia/studia"
MINCONTIGLENGTH=200

THREADS=8

#SAMPLES=glob.glob(INDIR+"/*_depl_1.fq.gz")+glob.glob(INDIR+"/*_depl_2.fq.gz")
SAMPLES = ["6683_16-06-2015", "6685_04-06-2015", "6685_16-06-2015", "6690_04-06-2015", "6690_16-06-2015", "6695_04-06-2015", "6695_16-06-2015", "6704_04-06-2015", "6704_16-06-2015"]
SUFFIX1 = "_depl_1.fq.gz"
SUFFIX2 = "_depl_2.fq.gz"

KALLISTO="/home/julia/kallisto_kod/src/kallisto"
BOWTIEDIR="/home/julia/lib/bowtie2-2.2.9"
VELVETDIR="/home/julia/velvet_1.2.10"
OASESDIR="/home/julia/oases"
METAVELVETDIR="/home/julia/MetaVelvet-1.2.02"
MEGAHITDIR="/home/julia/megahit"
IDBADIR="/home/julia/lib/idba-1.1.3/bin"
SGADIR="/home/julia/sga"
