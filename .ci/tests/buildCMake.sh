#!/bin/sh
help()
{
  echo "./buildCMake.sh as_host workingdir [options] [-- <hostenv.sh options>]"
  echo "  as_host                   First argument must be the host configuration to use for environment loading"
  echo "  workingdir                Second argument must be the working dir to immediate cd to"
  echo "  -c                        Configuration build type, piped directly into configure"
  echo "  -b                        Build command passed into compile"
  echo "  -r                        Clean command passed into cleanCmake"
  echo "  -e                        environment variables in comma-delimited list, e.g. var=1,foo,bar=0"
  echo "  -- <hostenv.sh options>   Directly pass options to hostenv.sh, equivalent to hostenv.sh <options>"
  echo "  -h                  Print this message"
  echo ""
  echo "If you wish to use an env var in your arg such as '-c \$PRESELECT -e PRESELECT=GNU',"
  echo "you will need to do '-c \\\$PRESELECT -e PRESELECT=32' to delay shell expansion"
}


echo "Input arguments:"
echo "$*"

AS_HOST=$1
shift
if [ $AS_HOST = "-h" ]; then
  help
  exit 0
fi

workingDirectory=$1
shift

cd $workingDirectory

# Get some helper functions
. .ci/env/helpers.sh

while getopts c:b:r:e:h opt; do
  case $opt in
    c)
      configuration="$OPTARG"
    ;;
    b)
      buildCommand="$OPTARG"
    ;;
    r)
      cleanCommand="$OPTARG"
    ;;
    e)
      envVars="$envVars,$OPTARG"
    ;;
    h)  help; exit 0 ;;
    *)  help; exit 1 ;;
    :)  help; exit 1 ;;
    \?) help; exit 1 ;;
  esac
done

shift "$((OPTIND - 1))"

# Everything else goes to our env setup
. .ci/env/hostenv.sh $*

# Now evaluate env vars in case it pulls from hostenv.sh
if [ ! -z "$envVars" ]; then
  setenvStr "$envVars"
fi

# Re-evaluate input values for delayed expansion
eval "cleanCommand=\"$cleanCommand\""
eval "configuration=\"$configuration\""
eval "buildCommand=\"$buildCommand\""

echo "./cleanCMake.sh -a $cleanCommand"
./cleanCMake.sh -a $cleanCommand
echo "Clean done"

echo "./configure_new $configuration"
./configure_new $configuration
result=$?
echo "Configure done"

if [ $result -ne 0 ]; then
  echo  "Failed to configure, command returned non-zero exit status"
  exit 1
fi

echo "./compile_new $buildCommand"
./compile_new $buildCommand
result=$?
echo "Compile done"



if [ $result -ne 0 ]; then
  echo "Failed to compile, command returned non-zero exit status"
  exit 1
fi

echo "TEST $(basename $0) PASS"
