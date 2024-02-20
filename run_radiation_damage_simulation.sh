
simulation_code_path=src/simulation/

# Command to compile the code, need to do this after each code change
g++ ${simulation_code_path}/run_radiation_damage_simulation.cpp -o ${simulation_code_path}/run_radiation_damage_simulation -std=c++11 -I${LCG_BASE_PATH}/include -L${LCG_BASE_PATH}/lib -lboost_program_options -Wall `root-config --cflags --libs`

if [ "$?" != "0" ]; then
    echo -e "\n\033[1;31mCompilation FAILED!\033[0m"
    exit 1
else
    echo -e "\n\033[1;32mCompilation SUCCESSFUL!\033[0m"
fi

pixel_monitoring_path=../PixelMonitoring
sector=BPix_BmI_SEC1_LYR1
#period_type=per_year  # per_year or per_phase
#period=2022
period_type=per_phase  # per_year or per_phase
period=phase1

profile_input_path=${pixel_monitoring_path}/data/radiation_simulation/profiles/${period_type}/${sector}/profile_${sector}_${period}.txt
output_simulation_directory=data/simulation
output_simulation_file_name=${output_simulation_directory}/simulation_${sector}_${period}.root

if [ ! -d ${output_directory} ]; then
    mkdir -p ${output_directory}
fi

${simulation_code_path}/run_radiation_damage_simulation -i ${profile_input_path} -o ${output_simulation_file_name}

