#! /bin/bash -e

PATH_TO_SELF="$(pwd $( dirname "${BASH_SOURCE[0]}" ) )" 
cd "$( dirname "${BASH_SOURCE[0]}" )" 

if [ -z "$1" ]; then
	INSTALL_DIR=$PATH_TO_SELF/miniconda2
	echo ""
	echo "Select the path where to install the miniconda instance:"
	echo ""
	read -r -p "[$INSTALL_DIR] >>> " response
	if [[ ! -z "$response" ]]; then
		INSTALL_DIR=$response
	fi
else
	INSTALL_DIR=$1
fi




echo ''
echo 'Installing Miniconda. Please review the EULA at https://docs.continuum.io/anaconda/eula'
read -r -p "Do you agree? [y/N] " response
case $response in
    [yY][eE][sS]|[yY] ) 
    echo ''
    ;;
    * )
    echo 'Cannot Proced'
    exit 1
    ;;
esac



unset PYTHONPATH
unset PYTHONHOME
case `uname` in
    Linux )
        if [ `uname -m` != "x86_64" ]; then
            echo 'Installation only automated on 64 bit systems'
            exit 1
        fi

        if hash wget 2>/dev/null; then 
            wget "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
        elif hash curl 2>/dev/null; then
            curl "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh" \
                -o "Miniconda3-latest-Linux-x86_64.sh"
        else 
            echo "ERROR: SimNIBS installer needs either the curl of wget packages"
            exit 1; 
        fi

        bash Miniconda3-latest-Linux-x86_64.sh -b -p $INSTALL_DIR

        ;;
    Darwin )
        curl "https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh" \
            -o "Miniconda3-latest-MacOSX-x86_64.sh"
        bash Miniconda3-latest-MacOSX-x86_64.sh -b -f -p $INSTALL_DIR
        ;;
    * )
        echo "SimNIBS is only avaliable for Linux and MacOSX Systems"
        exit 1
        ;;
esac


rm Miniconda*
