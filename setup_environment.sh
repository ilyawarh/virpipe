# Step 1: Create a new environment
mamba create -n virpipe_short2 -c bioconda -c conda-forge -c default -c ursky -y \
    python=3.9 \
    snakemake \
    fastqc \
    trim-galore \
    kaiju \
    spades \
    centrifuge \
    samtools \
    git

# Activate the environment
mamba activate virpipe_short2

#Make needed directories
mkdir results
mkdir tools

# Step 2: Install KrakenTools (requires git and manual installation)
cd tools
git clone https://github.com/jiarong/VirSorter2.git
cd VirSorter2
pip install -e ./


# Step 5: Install additional Python dependencies if needed
pip install pandas numpy
