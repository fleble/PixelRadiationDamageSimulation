lxplus7_regex="lxplus7"
lxplus9_regex="lxplus9"

if [[ $HOSTNAME =~ ${lxplus9_regex} ]]; then
  export LCG_BASE_PATH=/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-clang16-dbg
elif  [[ $HOSTNAME =~ ${lxplus7_regex} ]]; then
  export LCG_BASE_PATH=/cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt
else
  echo "Unsupported hostname ${HOSTNAME}". Only supporting lxplus 7 and 9.
fi

source ${LCG_BASE_PATH}/setup.sh

export PIXEL_RADSIM_DIR=${PWD}

export PYTHONPATH=${PYTHONPATH}:${PWD}
export PYTHONPATH=${PYTHONPATH}:${PWD}/src

