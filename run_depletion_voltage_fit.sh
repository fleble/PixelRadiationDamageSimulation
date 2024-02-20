
simulation_code_path=src/simulation/

sector=BPix_BmI_SEC1_LYR1
#period_type=per_year  # per_year or per_phase
#period=2022
period_type=per_phase  # per_year or per_phase
period=phase1

output_simulation_directory=data/simulation
output_simulation_file_name=${output_simulation_directory}/simulation_${sector}_${period}.root

output_depletion_voltage_directory=data/depletion_voltage_fit
output_depletion_voltage_file_name=${sector}_${period}
depletion_voltage_file_name=data/depletion_voltage_measurements/${period}.py

python src/depletion_voltage_fit/fit_depletion_voltage_data.py -i ${output_simulation_file_name} -dv ${depletion_voltage_file_name} -od ${output_depletion_voltage_directory} -o ${output_depletion_voltage_file_name} -ymin 0 -ymax 800

