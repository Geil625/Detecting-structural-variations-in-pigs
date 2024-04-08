#!/bin/bash
pig=... # the pig individual under consideration
callset=... # the callset in evaluation
complex_region=... # the specified complex region

### overall calculation
truvari bench -p 0 -r 1000 \
-b ${pig}/benchmark_set.vcf.gz \
-c ${callset} \
-o unreference/${callset}/

### specific calculation
truvari bench -p 0 -r 1000 \
-b ${pig}/benchmark_set.vcf.gz \
-c ${callset} \
--includebed ${complex_region} \
-o unreference/${callset}/
