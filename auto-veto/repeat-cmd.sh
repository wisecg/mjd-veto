#!/bin/bash

make --quiet
for run in {1..10}
do
  ./skim-coins -l 3 runs/test.txt skimout
done

