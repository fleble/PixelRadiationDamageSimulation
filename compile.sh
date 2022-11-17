g++ run_radiation_damage_simulation.cpp -o run_radiation_damage_simulation -std=c++11 -I${LCG_BASE_PATH}/include -L${LCG_BASE_PATH}/lib -lboost_program_options -Wall `root-config --cflags --libs`
