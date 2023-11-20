#!/bin/bash

source_dir="/home/users/jcfan/MyPrograms/Hill_energy_per_cluster/src/"
output_dir="/home/users/jcfan/MyPrograms/Hill_energy_per_cluster/build"
build_name="hill_sherlock"

gcc "${source_dir}"/*.c -o "${output_dir}/${build_name}" -lm
