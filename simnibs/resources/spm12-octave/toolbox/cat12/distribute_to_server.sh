 #! /bin/sh

########################################################
# global parameters
########################################################
version='distribute_to_server.sh $Id: distribute_to_server.sh 764 2015-11-17 13:11:53Z gaser $'

COMMAND=""
SERVER=localhost
PATTERN=""
USER=`whoami`
DIR=""

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  
  check_files
  distribute

  exit 0
}

########################################################
# check arguments and files
########################################################

parse_args ()
{
  local optname optarg files_to_calculate all_files
  count=0
  files_to_calculate=0
  all_files=0
  while [ $# -gt 0 ]
  do
    optname="`echo $1 | sed 's,=.*,,'`"
    optarg="`echo $2 | sed 's,^[^=]*=,,'`"
    case "$1" in
        --command* | -c*)
            exit_if_empty "$optname" "$optarg"
            COMMAND=$optarg
            shift
            ;;
        --server* | -s*)
            exit_if_empty "$optname" "$optarg"
            SERVER=$optarg
            shift
            ;;
        --pattern* | -p*)
            exit_if_empty "$optname" "$optarg"
            PATTERN=$optarg
            shift
            ;;
        --dir* | -u*)
            exit_if_empty "$optname" "$optarg"
            USER=$optarg
            shift
            ;;
        --dir* | -d*)
            exit_if_empty "$optname" "$optarg"
            DIR=$optarg
            shift

            if [ -z "$PATTERN" ]; then
              echo Pattern have to be defined first to use that function.
              exit 0
            fi

            # exclude that patterns from search
            list=`find $DIR -name "*.[in][mi][gi]" \! -name "*wrp[0-3]*.nii"  \! -name "*wp[0-3]*.nii" \! -name "wm*.nii"   \! -name "wrm*.nii"  \! -name "bf*.nii"  \! -name "p[0-3]*.nii"  \! -name "iy_*.nii"  \! -name "y_*.nii"  \! -name "rp[0-3]*.nii"`

            for i in ${list} ; do
              ext="${i##*.}"
              name="${i%.*}"
              # remove leading "./"
              name=`echo $name|sed -e 's/\.\///g'`
              bname="${name##*/}"
              dname="${name%/*}"
              
              for j in ${dname}/${PATTERN}${bname}*.nii ; do
                if [ ! -f "$j" ]; then
                  ARRAY[$count]=$i
                  ((count++))
                fi
              done
            done
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

  if [ "$count" == "0" ] && [ -z "PATTERN" ] ; then
    echo All files are already processed.
    exit 0
  fi
  
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
    echo ERROR: "No argument given with \"$desc\" command line argument!" >&2
    exit 1
  fi
}

########################################################
# check files
########################################################

check_files ()
{
  if [ -z "$COMMAND" ];
  then
    echo "$FUNCNAME ERROR - no command defined."
      help
    exit 1
  fi
  
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
    if [ ! -f "${ARRAY[$i]}" -a ! -L "${ARRAY[$i]}" ] && [ ! -d "${ARRAY[$i]}" ]; then
      echo ERROR: File or directory ${ARRAY[$i]} not found
      help
      exit 1
    fi
    ((i++))
  done

}

########################################################
# run distribute
########################################################

distribute ()
{

  NUMBER_OF_SERVERS=0
  for k in ${SERVER}; do
    ((NUMBER_OF_SERVERS++))
  done

  SIZE_OF_ARRAY="${#ARRAY[@]}"
  BLOCK=$((10000* $SIZE_OF_ARRAY / $NUMBER_OF_SERVERS ))

  i=0
  while [ "$i" -lt "$SIZE_OF_ARRAY" ]
  do
    count=$((10000* $i / $BLOCK ))
    if [ -z "${ARG_LIST[$count]}" ]; then
      ARG_LIST[$count]="${ARRAY[$i]}"
    else
      ARG_LIST[$count]="${ARG_LIST[$count]} ${ARRAY[$i]}"
    fi
    ((i++))
  done
    
  i=0
  for x in ${SERVER};
  do
    if [ ! "${ARG_LIST[$i]}" == "" ]; then
      j=$(($i+1))
      echo job ${j}/"$NUMBER_OF_SERVERS":
      echo $COMMAND ${ARG_LIST[$i]}
      if [ "$x" == "localhost" ]; then
        $COMMAND ${ARG_LIST[$i]}
      else
        bash -c "ssh ${USER}@${x} $COMMAND ${ARG_LIST[$i]}"
      fi
    fi
    ((i++))
  done

}

########################################################
# help
########################################################

help ()
{
cat <<__EOM__

USAGE:
  distribute_to_server.sh [-s server] [-p pattern] -c command_to_distribute_to_server.sh filename|filepattern|-d directory
  
   -c   command that should be distributed
   -s   server list (if empty the command runs on the local machine)
   -p   pattern to for search of already processed files that is prepended. 
   -u   user

   Only one filename or pattern or disrectory using the -d flag is allowed. This can be either a single file or a pattern
   with wildcards to process multiple files or even a directory. For the latter case you also have to define a pattern that
   is used for the search for already processed files. 

PURPOSE:
   distribute_to_server.sh a job or command

OUTPUT:

EXAMPLE
   distribute_to_server.sh -c "niismooth -v -fwhm 8" sTRIO*.nii
   smoothing with fwhm of 8mm for all files sTRIO*.nii. Use verbose mode to see diagnostic output.
   
   distribute_to_server.sh -s "141.35.68.68 141.35.68.72 141.35.68.73 141.35.68.74 141.35.68.75" -c "/Volumes/UltraMax/spm12/toolbox/cat12/cat_batch_cat.sh -p 8 -d /Volumes/UltraMax/cat_defaults_p0123.m -m /Volumes/UltraMax/MATLAB_R2010b.app/bin/matlab" /Volumes/UltraMax/SVE.LPBA40.testdata/S*.nii
   CAT12 batch for all files in /Volumes/UltraMax/SVE.LPBA40.testdata/S*.nii with 8 parallel jobs and optional default file

   distribute_to_server.sh -s "141.35.68.68 141.35.68.73 141.35.68.74 141.35.68.75" -c "/Volumes/UltraMax/spm12/toolbox/cat12/cat_batch_cat.sh -p 8 -d /Volumes/UltraMax/cat_defaults_m0wrp12.m -m /Volumes/UltraMax/MATLAB_R2010b.app/bin/matlab" -p m0wrp1 -d /Volumes/UltraMax/SVE.LPBA40.testdata
   CAT12 batch with 8 parallel jobs and optional default file. Only those files in /Volumes/UltraMax/SVE.LPBA40.testdata/ are processed where no prepended m0wrp1 pattern can be found. All other files are skipped.

   Using of MATLAB, SPM and VBM commands like the NLM filter function "cat_vol_sanlm({'file1','file2'})". CFILES contain the files for each job.  
     distribute_to_server.sh -s "141.35.68.96" -c "/Volumes/vbmDB/MRData/batches/cat_batch_cat.sh -m /Volumes/Ultramax/MATLAB_R2010b.app/bin/matlab -v /Volumes/Ultramax/spm12/toolbox/cat12/ -c \"cat_vol_sanlm\(CFILES\)\"" -p 8 -u local /Volumes/vbmDB/MRData/release20140211/pre/vbm8/INDI/HC/NYa/sub44*/INDI_*.nii

   Quality assurance for files processed by VBM8 preprocessing. Requires external noise correction by sanlm-call with 'n' prefix
     sh distribute_to_server.sh -s "141.35.68.96 141.35.68.95 141.35.68.75 141.35.68.74" -c "/Volumes/vbmDB/MRData/batches/spm12/toolbox/cat12/cat_batch_cat.sh -m /Volumes/Ultramax/MATLAB_R2010b.app/bin/matlab -v /Volumes/vbmDB/MRData/batches/spm12/toolbox/cat12 -p 4 -c \"cat_vol_sanlm\(struct\('data',CFILES,'prefix',\'n'\)\)\"" -u local /Volumes/vbmDB/MRData/release20140211/pre/vbm8/IXI/HC/*/*/p0*.nii
     sh distribute_to_server.sh -s "141.35.68.96 141.35.68.95 141.35.68.75 141.35.68.74" -c "/Volumes/vbmDB/MRData/batches/spm12/toolbox/cat12/cat_batch_cat.sh -m /Volumes/Ultramax/MATLAB_R2010b.app/bin/matlab -v /Volumes/vbmDB/MRData/batches/spm12/toolbox/cat12 -l /Volumes/vbmDB/MRData/log -p 4 -c \"cat_tst_qa2\(\'p0\',CFILES,struct\(\'prefix\',\'vbm8_\',\'write_csv\',0\)\)\"" -u local /Volumes/vbmDB/MRData/release20140211/pre/vbm8/IXI/HC/*/*/p0*.nii
    
This script was written by Christian Gaser (christian.gaser@uni-jena.de).
This is ${version}.

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}

