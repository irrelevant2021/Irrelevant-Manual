#!/bin/bash

for file in transformations/*.json; do
  relpath=${file:16}  # strip off "transformations/"
#  echo $relpath
  dirpath=${relpath%.*}  # strip off final ".json"
#  echo $dirpath
  openfe quickrun $file -o results/$relpath -d results/$dirpath
done
