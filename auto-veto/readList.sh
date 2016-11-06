#!/bin/bash

# I'm just tired of typing this long command over and over.

cat runs/p3jdy-complete.txt | while read -r line; do echo $line; qsub submit.sh $line; done

