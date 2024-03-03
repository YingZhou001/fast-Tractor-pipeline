# Install htslib

Download lastest htslib and install to local directory

```bash
wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2
tar xvf htslib-1.19.1.tar.bz2
mkdir htslib-1.19.1-install
# show full path of installation location
ls htslib-1.19.1-install
cd htslib-1.19.1
./configure prefix==fullpath/to/htslib-1.19.1-install
make
make install
```
