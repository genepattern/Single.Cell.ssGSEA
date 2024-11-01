#!bin/bash

#echo "/home/ubuntu/run_sc_ssgsea.py --input_file $1 --gene_set_database_file $2 --output_file_name $3 --n_threads $4 --cluster_data_label $5 --chip_file $6"

#source /home/ubuntu/sc_ssgsea_venv/bin/activate && python3 /home/ubuntu/run_sc_ssgsea.py \
#	--input_file $1 \
#	--gene_set_database_file $2 \
#	--output_file_name $3 \
#	--n_threads $4 \
#	--cluster_data_label $5 \
#	--chip_file $6

cmd="python3 /home/ubuntu/run_sc_ssgsea.py --input_file $1 --gene_set_database_file $2 --output_file_name $3 --n_threads $4 --cluster_data_label $5"

if [[ $6 -ne 0 ]]; then
	cmd="${cmd} --chip_file $6"
fi

echo cmd

source /home/ubuntu/sc_ssgsea_venv/bin/activate && eval $cmd