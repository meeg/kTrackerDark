#Set locations
export KTRACKER_ROOT="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
export KTRACKER_LIB=$KTRACKER_ROOT/lib-opt  #todo test for opt/dbg

#make the lib directory
if [ ! -d "$KTRACKER_LIB" ]; then
  mkdir "$KTRACKER_LIB"
fi

#Set libs
#remove existing references to kTracker
if [ -f $SEAQUEST_SETUP_ROOT/setup_bash_utils.sh ]; then
  source $SEAQUEST_SETUP_ROOT/setup_bash_utils.sh
  export LD_LIBRARY_PATH=`minidropit $LD_LIBRARY_PATH seaquest/ktracker`
fi
export LD_LIBRARY_PATH=$KTRACKER_LIB:$LD_LIBRARY_PATH
