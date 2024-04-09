#!/bin/bash

pig=... # the pig individual under consideration
depths=(1 5 10 20 30)
randomseed=22345

for depth in ${depths[*]};
do
	rasusa -i .../${pig}/HiFi/integrated.fastq.gz -c $depth -g 2501912388 -s $randomseed -o .../${pig}/HiFi/${depth}x/subset.fastq.gz # HiFi data downsampling
	rasusa -i .../${pig}/ONT/combined.fastq.gz -c $depth -g 2501912388 -s $randomseed -o .../${pig}/ONT/${depth}x/subset.fastq.gz # ONT data downsampling
done
