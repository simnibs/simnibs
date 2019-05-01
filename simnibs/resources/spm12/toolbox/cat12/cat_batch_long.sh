#! /bin/sh

########################################################
# global parameters
########################################################
# $Id: cat_batch_long.sh 1270 2018-02-07 13:32:14Z gaser $

matlab=matlab     # you can use other matlab versions by changing the matlab parameter
defaults_file=""
LOGDIR=$PWD
nojvm=""
output_surface=0
fg=

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  check_matlab
  check_files
  run_batch

  exit 0
}


########################################################
# check arguments and files
########################################################

parse_args ()
{
  local optname optarg
  count=0

  if [ $# -lt 1 ]; then
    help
    exit 1
  fi

  while [ $# -gt 0 ]
  do
	optname="`echo $1 | sed 's,=.*,,'`"
	optarg="`echo $2 | sed 's,^[^=]*=,,'`"
	case "$1" in
        --matlab* | -m*)
			exit_if_empty "$optname" "$optarg"
			matlab=$optarg
			shift
			;;
        --defaults_file* | -d*)
            exit_if_empty "$optname" "$optarg"
            defaults_file=$optarg
            shift
            ;;
        --surface* | -s*)
            exit_if_empty "$optname" "$optarg"
            output_surface=1
            ;;
        --nojvm | -n*)
            exit_if_empty "$optname" "$optarg"
            nojvm=" -nojvm "
            ;;
        --fg* | -fg*)
            exit_if_empty "$optname" "$optarg"
            fg=1
            ;;
        --logdir* | -l*)
            exit_if_empty "$optname" "$optarg"
            LOGDIR=$optarg
            if [ ! -d $LOGDIR ] 
            then
              mkdir -p $LOGDIR
            fi
            shift
            ;;
		-h | --help | -v | --version | -V)
			help
			exit 1
			;;
		-*)
			echo "`basename $0`: ERROR: Unrecognized option \"$1\"" >&2
			;;
        *)
            ARRAY[$count]=$1
            ((count++))
            ;;
	esac
    shift
  done

}

########################################################
# check arguments
########################################################

exit_if_empty ()
{
  local desc val

  desc="$1"
  shift
  val="$*"

  if [ -z "$val" ]
  then
	echo 'ERROR: No argument given with \"$desc\" command line argument!' >&2
	exit 1
  fi
}

########################################################
# check files
########################################################

check_files ()
{
  
  SIZE_OF_ARRAY="${#ARRAY[@]}"
  if [ "$SIZE_OF_ARRAY" -eq 0 ]
  then
      echo 'ERROR: No files given!' >&2
      help
      exit 1
  fi

  i=0
  while [ "$i" -lt "$SIZE_OF_ARRAY" ]
  do
    if [ ! -f "${ARRAY[$i]}" ]; then
      if [ ! -L "${ARRAY[$i]}" ]; then
        echo ERROR: File ${ARRAY[$i]} not found
        help
        exit 1
      fi
    fi
    ((i++))
  done

}

########################################################
# run batch
########################################################

run_batch ()
{
    cwd=`dirname $0`
    pwd=$PWD
    
    # we have to go into toolbox folder to find matlab files
    cd $cwd
    
    spm12=`dirname $cwd`
    spm12=`dirname $spm12`

    if [ "${LOGDIR}" == "" ]; then
        LOGDIR=`dirname ${ARRAY[0]}`
    fi
    
    export MATLABPATH=$spm12

    SIZE_OF_ARRAY="${#ARRAY[@]}"
	
    # argument empty?
    if [ ! "${defaults_file}" == "" ]; then
        # check wether absolute or relative names were given
        if [ ! -f ${defaults_file} -a -f ${pwd}/${defaults_file} ]; then
            defaults_file=${pwd}/${defaults_file}
        fi
    
        # check whether defaults file exist
        if [ ! -f ${defaults_file} ];  then
            echo $defaults_file not found.
            exit
        fi
    fi

    TMP=/tmp/cat_$$
    i=0
    ARG_LIST=""
    while [ "$i" -lt "$SIZE_OF_ARRAY" ]
    do
        # check wether absolute or relative names were given
        if [ ! -f ${ARRAY[$i]} ];  then
            FILE=${pwd}/${ARRAY[$i]}
        else
            FILE=${ARRAY[$i]}
        fi
        if [ -z "${ARG_LIST}" ]; then
            ARG_LIST="$FILE"
        else
            ARG_LIST="${ARG_LIST} $FILE"
        fi
        ((i++))
    done
    
    echo ${ARG_LIST} >> ${TMP}

	time=`date "+%Y%b%d_%H%M"`
    vbmlog=${LOGDIR}/cat_${HOSTNAME}_${time}.log
	echo Check $vbmlog for logging information
	echo

    COMMAND="cat_batch_long('${TMP}','${output_surface}','${defaults_file}')"
	echo Running ${ARG_LIST}
	echo > $vbmlog
	echo ---------------------------------- >> $vbmlog
	date >> $vbmlog
	echo ---------------------------------- >> $vbmlog
	echo >> $vbmlog
	echo $0 $ARG_LIST >> $vbmlog
	echo >> $vbmlog

    if [ -z "$fg" ]; then
      nohup ${matlab} "$nojvm" -nodisplay -nosplash -r "$COMMAND" >> $vbmlog 2>&1 &
    else
      nohup ${matlab} "$nojvm" -nodisplay -nosplash -r "$COMMAND" >> $vbmlog 2>&1
    fi
    
	exit 0
}

########################################################
# check if matlab exist
########################################################

check_matlab ()
{
  found=`which ${matlab} 2>/dev/null`
  if [ ! -n "$found" ];then
    echo $matlab not found.
    exit 1
  fi
}

########################################################
# help
########################################################

help ()
{
cat <<__EOM__

USAGE:
   cat_batch_long.sh file1.nii file2.nii ... filex.nii [-d] [-m matlabcommand]
   
   -m       matlab command
   -d       optional default file
   -fg      do not run matlab process in background
   -surface enable surface and thickness estimation
   -nojvm   supress call of jvm using the -nojvm flag

   Only one batch filename is allowed. Optionally you can set the matlab command 
   with the "-m" option. As default no display is used (via the -nodisplay option 
   in matlab). However sometimes the batch file needs a graphical output and the 
   display should be enabled with the option "-d".

PURPOSE:
   Command line call of SPM12 batch files

EXAMPLE
   cat_batch_long.sh all_files*.nii -m /Volumes/UltraMax/MATLAB_R2010b.app/bin/matlab
   This command will process all given files in the longitudinal pipeline. As matlab command 
   /Volumes/UltraMax/MATLAB_R2010b.app/bin/matlab will be used.
   
INPUT:
   filenames

OUTPUT:
   ${LOGDIR}/spm_${HOSTNAME}_${time}.log for log information

USED FUNCTIONS:
   SPM12

SETTINGS
   matlab command: $matlab
   
This script was written by Christian Gaser (christian.gaser@uni-jena.de).

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}
