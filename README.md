K-mer counting and [KAT](https://github.com/TGAC/KAT)-inspired k-mer specturm plot for metagenome assembly evaluation. 

## Getting Started

```sh
# Download and compile
git clone https://github.com/xfengnefx/yam
cd yam && make

# build k-mer hash table for assembly and reads
yam ht -t 48 -o sample asm.fa reads.fa
# build k-mer hash table for only assembly
yam ht -t 48 -2 -o sample asm.fa -
# build k-mer hash table for only reads
yam ht -t 48 -1 -o sample - reads.fa 

# collect k-mer spectrum
yam compare -o kspec sample.asm.hashtables sample.reads.hashtables

# plot, will generate a pdf file (dependency: numpy)
python misc/plotkmerspectrum.py kspec.mat kspec.pdf
```
