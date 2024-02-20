python Vdepl_data.py -i sim_r1d1rog1BpO_vdepl_clean_jul.root -o Vdepl_sim_r1d1rog1BpO_clean_all.root -l d1 -d Vdepl_data/CMS_PIXELS_VDEPL.csv -b --ring 1 --disk 1 --all &
python Vdepl_data.py -i sim_r1d3rog1BpO_vdepl_clean_jul.root -o Vdepl_sim_r1d3rog1BpO_clean_all.root -l d2 -d Vdepl_data/CMS_PIXELS_VDEPL.csv -b --ring 1 --disk 3 --all &
python Vdepl_data.py -i sim_r1d1rog1BpO_vdepl_clean_jul.root -o Vdepl_sim_r1d1rog1BpO_clean_log.root -l d1 -d Vdepl_data/CMS_PIXELS_VDEPL.csv -b --ring 1 --disk 1 --log &
python Vdepl_data.py -i sim_r1d3rog1BpO_vdepl_clean_jul.root -o Vdepl_sim_r1d3rog1BpO_clean_log.root -l d2 -d Vdepl_data/CMS_PIXELS_VDEPL.csv -b --ring 1 --disk 3 --log &
python Vdepl_data.py -i sim_r1d1rog1BpO_vdepl_clean_jul.root -o Vdepl_sim_r1d1rog1BpO_clean.root -l d1 -d Vdepl_data/CMS_PIXELS_VDEPL.csv -b --ring 1 --disk 1 &
python Vdepl_data.py -i sim_r1d3rog1BpO_vdepl_clean_jul.root -o Vdepl_sim_r1d3rog1BpO_clean.root -l d2 -d Vdepl_data/CMS_PIXELS_VDEPL.csv -b --ring 1 --disk 3 &
