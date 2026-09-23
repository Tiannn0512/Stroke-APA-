#!/bin/bash
awk -F'\t' '$3=="UTR" || $3=="three_prime_UTR" {print $3}' /home/taylor/reference_v2/gencode.vM25.annotation.gtf | sort | uniq -c
