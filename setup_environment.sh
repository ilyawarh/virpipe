# Step 1: Create a new environment
mamba create -n viral_pipeline -y \
    python=3.9 \
    snakemake \
    fastqc \
    trimgalore \
    kraken2 \
    spades \
    centrifuge \
    samtools \
    git

# Activate the environment
mamba activate viral_pipeline

# Step 2: Install KrakenTools (requires git and manual installation)
git clone https://github.com/jenniferlu717/KrakenTools.git
cd KrakenTools
chmod +x *.py
# Add KrakenTools directory to PATH for easier use
export PATH=$(pwd):$PATH
cd ..

# Step 3: Install VirSorter2 via conda
mamba install -c bioconda -c conda-forge virsorter2 -y

# Step 4: Install DeepVirFinder (requires pip)
pip install deepvirfinder

# Step 5: Install additional Python dependencies if needed
pip install pandas numpy

# Step 6: Verify installations
echo "Verifying tool installations..."
echo "FastQC: $(fastqc --version)"
echo "Trim Galore: $(trim_galore --version)"
echo "Kraken2: $(kraken2 --version)"
echo "SPAdes: $(spades.py --version)"
echo "Centrifuge: $(centrifuge --version)"
echo "VirSorter2: $(virsorter run --help | head -n 1)"
echo "DeepVirFinder: $(run_deepvirfinder.py --version 2>&1 | grep 'DeepVirFinder')"

echo "Environment setup is complete. Activate it with 'mamba activate viral_pipeline'."

