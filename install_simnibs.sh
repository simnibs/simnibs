#!/bin/bash
set -e
VERSION=3.0.0
DEV=false
if [ "$1" = "develop" ]; then
	DEV=true
fi
cd "$( dirname "${BASH_SOURCE[0]}" )" 
if [ "$DEV" = true ]; then
	echo ''
	echo "Installing in development mode"

	TARGET_DIR="$(pwd $( dirname "${BASH_SOURCE[0]}" ) )" 
	MINICONDADIR=$TARGET_DIR/miniconda3
	echo ""
	echo "Select the path where to install the miniconda instance:"
	echo ""
	read -r -p "[$MINICONDADIR] >>> " response
	if [[ ! -z "$response" ]]; then
		MINICONDADIR=$response
	fi
	MINICONDADIR="$(cd "$(dirname "$MINICONDADIR")"; pwd)/$(basename "$MINICONDADIR")"
else
	echo ''
	echo "The default installation directory is:"
	echo "$HOME/simnibs_$VERSION"
	echo ""
	echo "  -  Press ENTER to confirm the location"
	echo "  -  Press CTR-C to abort the installation"
	echo "  -  Or specify another location below"
	echo ""

	read -r -p "[$HOME/simnibs_$VERSION/] >>> " response

	if [ -z "$response" ]; then
		TARGET_DIR="$HOME/simnibs_$VERSION/"
	    echo ''
	else
		TARGET_DIR=$response
	fi
	TARGET_DIR="$(cd "$(dirname "$TARGET_DIR")"; pwd)/$(basename "$TARGET_DIR")"
	if [ ! -d "$TARGET_DIR" ]; then
		mkdir -p $TARGET_DIR
	fi
	cp -r ./ $TARGET_DIR
	MINICONDADIR=$TARGET_DIR/miniconda3
fi

/bin/bash install_miniconda_instance.sh $MINICONDADIR
if [ "$DEV" = true ]; then
	$MINICONDADIR/bin/conda env create -f dev_environment.yml
else
	$MINICONDADIR/bin/conda env create -f environment.yml
fi
cd $TARGET_DIR/Python_modules/src
$MINICONDADIR/envs/simnibs_env/bin/python3 setup.py develop
cd ../pygpc
$MINICONDADIR/envs/simnibs_env/bin/python3 setup.py develop
cd ../../

/bin/bash create_wrappers.sh $MINICONDADIR/envs/simnibs_env/bin/python3 
ln -s -f $MINICONDADIR/envs/simnibs_env/bin/python3 bin/simnibs_python
if [ "$DEV" = true ]; then
	ln -s -f $MINICONDADIR/envs/simnibs_env/bin/py.test bin/simnibs_test
fi

echo 'Changing Gmsh default configuration to SimNIBS standard'
cp gmsh-options_simnibsdefault ~/.gmsh-options

case `uname` in
	Linux )
		cd bin
		ln -s -f `which expr` expr
		ln -s -f `which getopt` getopt
		ln -s -f linux/gmsh gmsh
		ln -s -f linux/meshfix meshfix
		cd ..
	;;
	Darwin )
		cd bin
		ln -s -f osx/expr expr
		ln -s -f osx/getopt getopt
		ln -s -f osx/gmsh gmsh
		ln -s -f osx/meshfix meshfix
		cd ..
	;;
esac

echo ""
case `uname` in
    Linux )
		echo 'Modifing ~/.bashrc'
		echo 'Backing up ~/.bashrc to ~/.bashrc.simnibs.bk'
		if [ -f ~/.bashrc ]; then
			cp ~/.bashrc ~/.bashrc.simnibs.bk
		fi
		if grep -Fq "SIMNIBSDIR" ~/.bashrc ; then 
			echo "Another SIMNIBS installation detected, overwriting \$SIMNIBSDIR in the .bashrc file"
			sed -i.simnibs.bak -e "s|export SIMNIBSDIR=.*|export SIMNIBSDIR=$TARGET_DIR|" ~/.bashrc
			sed -i.simnibs.bak -e "s|source \$SIMNIBSDIR/.*|source \$SIMNIBSDIR/simnibs_conf.sh|" ~/.bashrc
		else
		    echo "">> ~/.bashrc
		    echo "export SIMNIBSDIR=$TARGET_DIR">> ~/.bashrc
		    echo "source \$SIMNIBSDIR/simnibs_conf.sh">> ~/.bashrc
		fi
    ;;
    Darwin )
		echo 'Modifing ~/.bash_profile'
		echo 'Backing up ~/.bash_profile to ~/.bash_profile.simnibs.bk'
		if [ -f ~/.bash_profile ]; then
			cp ~/.bash_profile ~/.bash_profile.simnibs.bk
		fi
		if grep -Fq "SIMNIBSDIR" ~/.bash_profile ; then 
			echo "Another SIMNIBS installation detected, overwriting \$SIMNIBSDIR in the .bash_profile file"
			sed -i.simnibs.bak -e "s|export SIMNIBSDIR=.*|export SIMNIBSDIR=$TARGET_DIR|" ~/.bash_profile
			sed -i.simnibs.bak -e "s|source \$SIMNIBSDIR/.*|source \$SIMNIBSDIR/simnibs_conf.sh|" ~/.bash_profile
		else
			echo "">> ~/.bash_profile
			echo "export SIMNIBSDIR=$TARGET_DIR">> ~/.bash_profile
			echo "source \$SIMNIBSDIR/simnibs_conf.sh">> ~/.bash_profile
		fi
    ;;
esac

if [ -z "$FSL_DIR" ]; then
    echo ""
    echo "Could not find FSL installation"
    echo "mri2mesh will not work"

fi


if [ -z "$FREESURFER_HOME" ]; then
    echo ""
    echo "Could not find FreeSurfer installation"
    echo "mri2mesh will not work"

fi

echo ""
echo "SimNIBS $VERSION successfully installed"
echo "You can start the GUI by opening a new terminal window and typing simnibs_gui"

exit 0
