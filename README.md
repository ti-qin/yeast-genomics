# Usage
## Assembly
### 00 installation
```shell
mamba create -n Genome_as
mamba activate Genome_as
mamba install seqkit -c bioconda

## flye install 
mamba install flye -c bioconda
## necat install 
git clone https://github.com/xiaochuanle/NECAT.git
cd NECAT/src/
make
cd ../Linux-amd64/bin
export PATH=$PATH:$(pwd)
## nextdenovo
pip install paralleltask
wget https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz
tar -vxzf NextDenovo.tgz && cd NextDenovo
export PATH=$PATH:$(pwd)

```
### 01 flye
```shell
./01_genome_asm_fyle.sh \
-r ./raw_data/ont \
-o ./result/assembly/flye \
-t 32 \
-j 4 \
--overwrite
```
### 02 Necat
```shell
./02_genome_asm_necat.sh \
-r ./raw_data/ont \
-o ./result/assembly/necat \
-t 32 \
-j 4 \
--overwrite
```
### 03 Nextdenovo
```shell
./03_genome_asm_nextdenovo.sh \
-r ./raw_data/ont \
-w ./result/assembly/denovo \
-t 16 \
-j 2 \
-g 12000000
```     