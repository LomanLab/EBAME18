wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh -b
echo 'PATH="$HOME/miniconda3/bin:$HOME/bin:$HOME/.local/bin:$PATH"' >> .profile
source .profile
conda config --add channels defaults
conda config --add channels bioconda
conda install -y kraken2 minimap2 miniasm bracken
sudo mkdir /mnt/ebame18
chown ubuntu /mnt/ebame18
cd /mnt/ebame18
wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz
tar xvfz minikraken_20171019_8GB.tgz
echo 'export KRAKEN_DEFAULT_DB=/mnt/ebame18/minikraken_20171019_8GB' >> .profile
mkdir kraken2-microbial-fatfree/
cd kraken2-microbial-fatfree/
wget https://refdb.s3.climb.ac.uk/kraken2-microbial/hash.k2d
wget https://refdb.s3.climb.ac.uk/kraken2-microbial/opts.k2d
wget https://refdb.s3.climb.ac.uk/kraken2-microbial/taxo.k2d
wget https://refdb.s3.climb.ac.uk/kraken2-microbial/database.kraken
wget https://refdb.s3.climb.ac.uk/kraken2-microbial/database2500mers.kraken
wget https://refdb.s3.climb.ac.uk/kraken2-microbial/database2500mers.kmer_distrib
echo 'echo 'export KRAKEN2_DEFAULT_DB=/mnt/ebame18/kraken2-microbial-fatfree' >> .profile
cd ..
wget http://nanopore.s3.climb.ac.uk/Kefir_RBK.fastq
conda install -y anvio racon nanoplot
