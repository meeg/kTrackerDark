#Set locations
export KTRACKER_ROOT="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
export KTRACKER_LIB=$KTRACKER_ROOT/lib-opt  #todo test for opt/dbg

#make the lib directory
if [ ! -d "$KTRACKER_LIB" ]; then
  mkdir "$KTRACKER_LIB"
fi

#Set libs
#todo remove existing references to kTracker
export LD_LIBRARY_PATH=$KTRACKER_LIB:$LD_LIBRARY_PATH
