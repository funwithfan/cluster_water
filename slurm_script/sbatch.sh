#!/bin/bash

P=540
numFrame=2501
save_every=100
resultDir="/home/users/jcfan/sc_water/Hill_result"
simulationDir="/scratch/users/jcfan/sc_water/MD"

# 使用 seq 命令生成 T1 和 T2 数组的元素
T1=($(seq 300 25 1200))
T2=($(seq 600 10 800))

# 合并数组并排序去重
T_list=($(echo "${T1[@]}" "${T2[@]}" | tr ' ' '\n' | sort -nu))

num_jobs=${#T_list[@]}

for ((i=0; i<$num_jobs; i++)); do
    T=${T_list[i]}
    job_name="${P}_${T}"  
    output_file="${P}_${T}_%j.out" 
    error_file="${P}_${T}_%j.err" 
    time_limit="5:00:00"
    partition="normal"
    cpus="1"
    memory="2GB"

    # 编写 SLURM 脚本文件
    echo "#!/bin/bash" > slurm_script_${i}.sh
    echo "#SBATCH --job-name=${job_name}" >> slurm_script_${i}.sh
    echo "#SBATCH --output=${output_file}" >> slurm_script_${i}.sh
    echo "#SBATCH --error=${error_file}" >> slurm_script_${i}.sh
    echo "#SBATCH --time=${time_limit}" >> slurm_script_${i}.sh
    echo "#SBATCH -p ${partition}" >> slurm_script_${i}.sh
    echo "#SBATCH -c ${cpus}" >> slurm_script_${i}.sh
    echo "#SBATCH --mem=${memory}" >> slurm_script_${i}.sh

    echo "/home/users/jcfan/MyPrograms/Hill_energy_per_cluster/build/hill_sherlock ${P} ${T} ${numFrame} ${save_every} ${resultDir} ${simulationDir}" >> slurm_script_${i}.sh
    
    # 提交作业
    sbatch slurm_script_${i}.sh

    # 删除生成的 SLURM 脚本文件
    rm slurm_script_${i}.sh
done
