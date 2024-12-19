#!/bin/bash
sed -i  's/sampling_per_window = 4/sampling_per_window = 0.4/'  ./ties/ties*/*/TIES.cfg
sed -i  's/equili_per_window = 2/equili_per_window = 0.1/'  ./ties/ties*/*/TIES.cfg
sed -i  's/total_reps = 5/total_reps = 12/'  ./ties/ties*/*/TIES.cfg
sed -i  's/split_run = 1/split_run = 0/'  ./ties/ties*/*/TIES.cfg
