# Step 1: Create a new environment
mamba create -c bioconda -c conda-forge -c default -c ursky -n virpipe_short1 -y \
    python=3.9 \
    snakemake \
    fastqc \
    trim-galore \
    kraken2 \
    spades \
    centrifuge \
    samtools \
    git

# Activate the environment
mamba activate virpipe_short1

# Step 2: Install KrakenTools (requires git and manual installation)
mkdir tools
cd tools
git clone https://github.com/jenniferlu717/KrakenTools.git
cd KrakenTools
chmod +x *.py
cd ../

# Step 3: Install VirSorter2 via conda
git clone https://github.com/jiarong/VirSorter2.git
cd VirSorter2
pip install -e ./
cd ../

# Step 4: Install DeepVirFinder
git clone https://github.com/jessieren/DeepVirFinder.git
cd ../


# Step 5: Install additional Python dependencies if needed
pip install pandas numpy
