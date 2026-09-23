#!/bin/bash
echo "user=$(whoami) distro=$(lsb_release -ds 2>/dev/null)"
echo "--- net test ---"
curl -sI -m 15 https://ftp.ensembl.org/ | head -1 || echo "NET FAIL"
echo "--- tools ---"
for t in apt samtools bedtools subread-align featureCounts rna-star STAR hisat2 bowtie2 python3 pip3 curl wget gzip make gcc; do
  command -v "$t" >/dev/null 2>&1 && echo "OK $t" || echo "MISS $t"
done
echo "--- sudo ---"
sudo -n true 2>/dev/null && echo "SUDO passwordless OK" || echo "SUDO needs password or unavailable"
