#!/bin/bash

mkdir -p allfasta
awk '{filename = "allfasta/file" int((NR-1)/2); print > filename".txt"}' 4RNP_formatted.fa

