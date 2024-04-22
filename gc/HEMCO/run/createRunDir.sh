#!/bin/bash

# createRunDir.sh: Create HEMCO standalone run directory
#
# Optional argument: run directory name
#
# If optional run directory name argument is not passed then the user
# will be prompted to enter a name interactively, or choose to use the
# default name gchp_{simulation}/
#
# Usage: ./createRunDir.sh [rundirname]
#
# Initial version: E. Lundgren,10/5/2018

srcrundir=$(pwd)
cd ..
hemcodir=$(pwd)
cd ${srcrundir}

# Define separator lines
thickline="\n===========================================================\n"
thinline="\n-----------------------------------------------------------\n"

printf "${thickline}HEMCO STANDALONE RUN DIRECTORY CREATION${thickline}"

#-----------------------------------------------------------------
# Export data root path in ~/.geoschem/config if file exists
#-----------------------------------------------------------------
if [[ -f ${HOME}/.geoschem/config ]]; then
    source ${HOME}/.geoschem/config
    if [[ ! -d ${GC_DATA_ROOT} ]]; then
	printf "\nWarning: Default root data directory does not exist!"
        printf "\nSet new path below or manually edit ${HOME}/.geoschem/config.\n"
    fi
else
    printf "\nDefine path to ExtData."
    printf "\nThis will be stored in ${HOME}/.geoschem/config for future automatic use.\n"
    mkdir -p ${HOME}/.geoschem
fi

#-----------------------------------------------------------------
# One-time configuration of data root path in ~/.geoschem/config
#-----------------------------------------------------------------
if [[ -z "${GC_DATA_ROOT}" ]]; then
    printf "${thinline}Enter path for ExtData:${thinline}"
    valid_path=0
    while [ "$valid_path" -eq 0 ]; do
	read -e extdata
	if [[ ${extdata} = "q" ]]; then
	    printf "\nExiting.\n"
	    exit 1
	elif [[ ! -d ${extdata} ]]; then
            printf "\nError: ${extdata} does not exist. Enter a new path or hit q to quit.\n"
	else
	    valid_path=1
	    echo "export GC_DATA_ROOT=${extdata}" >> ${HOME}/.geoschem/config
            source ${HOME}/.geoschem/config
	fi
    done
fi

#-----------------------------------------------------------------
# Ask user to select meteorology source
#-----------------------------------------------------------------
printf "${thinline}Choose meteorology source:${thinline}"
printf "  1. MERRA-2 (Recommended)\n"
printf "  2. GEOS-FP\n"
printf "  3. GISS ModelE2.1 (GCAP 2.0)\n"
valid_met=0
while [ "${valid_met}" -eq 0 ]; do
    read met_num
    valid_met=1
    if [[ ${met_num} = "1" ]]; then
	met_name='MERRA2'
	met_name_lc="merra2"
	met_dir='MERRA2'
	met_native='0.5x0.625'
	met_latres='05'
	met_lonres='0625'
	met_extension='nc4'
	met_cn_year='2015'
	pressure_unit='Pa '
	pressure_scale='0.01'
    elif [[ ${met_num} = "2" ]]; then
	met_name='GEOSFP'
	met_name_lc="geosfp"
	met_dir='GEOS_FP'
	met_native='0.25x0.3125'
	met_latres='025'
	met_lonres='03125'
	met_extension='nc'
	met_cn_year='2011'
	pressure_unit='hPa'
	pressure_scale='1.0 '
    elif [[ ${met_num} = "3" ]]; then
	met_name='ModelE2.1'
	met_name_lc='modele2.1'
	met_dir='E21'
	met_resolution='2x2.5'
	met_native='2x2.5'
	met_latres='20'
	met_lonres='25'
	met_extension='nc4'
	met_cn_year='1950'
	pressure_unit='Pa '
	pressure_scale='0.01'
    else
	printf "Invalid meteorology option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to select horizontal resolution
#-----------------------------------------------------------------
printf "${thinline}Choose horizontal resolution:${thinline}"
if [[ ${met_name} = "ModelE2.1" ]] || [[ ${met_name} = "ModelE2.2" ]]; then
    printf "  1. 4.0  x 5.0 *\n"
    printf "  2. 2.0  x 2.5\n"
    printf "  3. 0.5  x 0.625 *\n"
    printf "  4. 0.25 x 0.3125 *${thinline}"
    printf "             * Will be interpolated online from native 2.0 x 2.5 resolution\n"
else
    printf "  1. 4.0 x 5.0\n"
    printf "  2. 2.0 x 2.5\n"
    printf "  3. 0.5 x 0.625\n"
    printf "  4. 0.25 x 0.3125\n"
    printf "  5. Custom\n"
fi

valid_res=0
while [ "${valid_res}" -eq 0 ]; do
    read res_num
    valid_res=1
    if [[ ${res_num} = "1" ]]; then
	grid_res='4x5'
	grid_res_long='4.0x5.0'
	grid_dir=$grid_res
	grid_file='HEMCO_sa_Grid.4x5.rc'
    elif [[ ${res_num} = "2" ]]; then
	grid_res='2x25'
	grid_res_long='2.0x2.5'
	grid_dir='2x2.5'
	grid_file='HEMCO_sa_Grid.2x25.rc'
    elif [[ ${res_num} = "3" ]]; then
	grid_res='05x0625'
	grid_res_long='0.5x0.625'
	grid_dir=$grid_res_long
	grid_file='HEMCO_sa_Grid.05x0625.rc'
    elif [[ ${res_num} = "4" ]]; then
	grid_res='025x03125'
	grid_res_long='0.25x0.3125'
	grid_dir=$grid_res_long
	grid_file='HEMCO_sa_Grid.025x03125.rc'
    elif [[ ${res_num} = "5" ]]; then
	printf "You will need to provide your own HEMCO_sa_Grid.rc file.\n"
	printf "See the HEMCO standalone guide for more information:\n"
	printf "http://wiki.seas.harvard.edu/geos-chem/index.php/HEMCO_standalone\n"
	valid_res=1
    else
	printf "Invalid resolution option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to provide path to HEMCO_Config.template file
#-----------------------------------------------------------------
printf "${thinline}Enter the file path to a HEMCO_Config.rc with your \n"
printf "emissions settings.\n\n"
printf " - This should be a HEMCO_Config.rc file from a \n"
printf "   pre-generated GEOS-Chem run directory and not a\n"
printf "   template config file from the GEOS-Chem repository.\n"
printf "\n"
printf " - If you do not have a pre-generated HEMCO_Config.rc file,\n"
printf "   type ./HEMCO_Config.rc.sample at the prompt below.\n"
printf "   This will copy a sample configuration file into your\n"
printf "   run directory.  You can then edit this configuration\n"
printf "   file with your preferred emission settings.${thinline}"

valid_path=0
while [ "$valid_path" -eq 0 ]; do
    read -e hco_config_file

    # Test for quitting
    if [[ "x${hco_config_file}" == "xq" ]]; then
	printf "\nExiting.\n"
	exit 1
    fi

    # Replace ~ with the user's home directory
    # NOTE: This is a safe algorithm.
    if [[ "${hco_config_file}" =~ '~' ]]; then
       hco_config_file="${hco_config_file/#\~/$HOME}"
       echo "Expanding to: ${hco_config_file}"
    fi

    if [[ ! -f ${hco_config_file} ]]; then
        printf "\nError: ${hco_config_file} does not exist. Enter a new file path or hit q to quit.\n"
    else
	valid_path=1
    fi

    if [[ "$hco_config_file" == *".rc"* ]]; then
	hco_config_dir=$(dirname $hco_config_file)
    fi
    
done

#-----------------------------------------------------------------
# Ask user to define path where directoy will be created
#-----------------------------------------------------------------
printf "${thinline}Enter path where the run directory will be created:${thinline}"

valid_path=0
while [ "$valid_path" -eq 0 ]; do
    read -e rundir_path

    # Test for quitting
    if [[ "x${rundir_path}" == "xq" ]]; then
	printf "\nExiting.\n"
	exit 1
    fi

    # Replace ~ with the user's home directory
    # NOTE: This is a safe algorithm.
    if [[ "${rundir_path}" =~ '~' ]]; then
       rundir_path="${rundir_path/#\~/$HOME}"
       echo "Expanding to: ${rundir_path}"
    fi

    # If this is just a new directory within an existing one,
    # give the user the option to proceed
    if [[ ! -d ${rundir_path} ]]; then
        if [[ -d $(dirname ${rundir_path} ) ]]; then
            printf "\nWarning: ${rundir_path} does not exist,\nbut the parent directory does.\nWould you like to make this directory? (y/n/q)\n"
            read mk_rundir
            if [[ "x${mk_rundir}" == "xy" ]]; then
                mkdir $rundir_path
	    elif [[ "x${mk_rundir}" == "xq" ]]; then
		printf "\nExiting.\n"
		exit 1
            fi
        fi
    fi

    # Ask user to supply a new path again
    if [[ ! -d ${rundir_path} ]]; then
        printf "\nERROR: ${rundir_path} does not exist. Enter a new path or hit q to quit.\n"
    else
	valid_path=1
    fi
done

#-----------------------------------------------------------------
# Ask user to define run directoy name if not passed as argument
#-----------------------------------------------------------------
if [ -z "$1" ]; then
    printf "${thinline}Enter run directory name, or press return to use default:\n\n"
    printf "NOTE: This will be a subfolder of the path you entered above.${thinline}"
    
    read -e rundir_name
    if [[ -z "${rundir_name}" ]]; then
	rundir_name=hemco_${grid_res}_"${met_name,,}"
	printf "  -- Using default directory name ${rundir_name}\n"
    fi
else
    rundir_name=$1
fi

#-----------------------------------------------------------------
# Ask user for a new run directory name if specified one exists
#-----------------------------------------------------------------
rundir=${rundir_path}/${rundir_name}
valid_rundir=0
while [ "${valid_rundir}" -eq 0 ]; do
    if [[ -d ${rundir} ]]; then
	printf "\nWarning: ${rundir} already exists.\n"
        printf "Enter a different run directory name, or q to quit:\n"
	read -e new_rundir
	if [[ ${new_rundir} = "q" ]]; then
	    printf "Exiting.\n"
	    exit 1
	else
	    rundir=${rundir_path}/${new_rundir}
	fi
    else
        valid_rundir=1
    fi
done

#-----------------------------------------------------------------
# Create run directory
#-----------------------------------------------------------------
mkdir -p ${rundir}

# Copy run directory files and subdirectories
cp -r ./OutputDir ${rundir}
cp ./HEMCO_sa_Config.template       ${rundir}/HEMCO_sa_Config.rc
cp ./HEMCO_sa_Time.rc               ${rundir}
cp ./HEMCO_sa_Spec.rc               ${rundir}
cp ./${grid_file}                   ${rundir}
cp ./runHEMCO.sh                    ${rundir}
cp ./README                         ${rundir}
cp ${hco_config_dir}/HEMCO_Config.* ${rundir}
if  [[ -f ${hco_config_dir}/HEMCO_Diagn.rc ]]; then
    cp ${hco_config_dir}/HEMCO_Diagn.rc ${rundir}
else
    printf "\nCould not find a HEMCO_Diagn.rc file corresponding to HEMCO_Config.rc!\n"
    printf "A sample HEMCO_Diagn.rc will be copied to the run directory.\n"
    cp ./HEMCO_Diagn.rc.sample ${rundir}/HEMCO_Diagn.rc
fi

# Create symbolic link to code directory
ln -s ${hemcodir} ${rundir}/CodeDir

# Create build directory
mkdir ${rundir}/build
printf "To build HEMCO type:\n   cmake ../CodeDir\n   make -j\n   make install\n" >> ${rundir}/build/README

#=============================================================================
#### Replacement for `sed -i -e` that works on both MacOS and Linux
#=============================================================================
function sed_ie() {
    REGEX=${1}
    FILE=${2}
    if [[ "x$(uname -s)" == "xDarwin" ]]; then
	sed -i '' -e "${REGEX}" "${FILE}"          # MacOS/Darwin
    else
	sed -i -e "${REGEX}" "${FILE}"             # GNU/Linux
    fi
}

#=============================================================================
#### Define function to replace values in config files
#=============================================================================
function replace_colon_sep_val() {
    KEY=${1}
    VALUE=${2}
    FILE=${3}

    # Debug print (leave commented out)
    # printf '%-30s : %-20s %-20s\n' "${KEY//\\}" "${VALUE}" "${FILE}"

    # Replace value in line starting with 'whitespace + key + whitespace + : +
    # whitespace + value' where whitespace is variable length including none
    #
    # MacOS sed does not allow you to use \t for tab.  The quick fix is
    # to use printf to save a tab character to a variable, and to use
    # that in the regular expression everywhere you would have used \t.
    # See the Github issue geoschem/geos-chem #617. (bmy, 2/23/21)
    TAB=$(printf "\t")
    REGEX="s|^\([${TAB} ]*${KEY}[${TAB} ]*:[${TAB} ]*\).*|\1${VALUE}|"
    sed_ie "${REGEX}" "${FILE}"
}

#--------------------------------------------------------------------
# Navigate to run directory and set up input files
#--------------------------------------------------------------------
cd ${rundir}

# Replace token strings in certain files
sed_ie "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   HEMCO_sa_Config.rc
sed_ie "s|{GRID_FILE}|${grid_file}|"      HEMCO_sa_Config.rc
sed_ie "s|{MET_NAME}|${met_name}|"        HEMCO_sa_Config.rc
sed_ie "s|{GRID_RES}|${grid_res}|"        HEMCO_sa_Config.rc
sed_ie "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   HEMCO_Config.rc
sed_ie "s|{GRID_DIR}|${grid_dir}|"        HEMCO_Config.rc
sed_ie "s|{MET_DIR}|${met_dir}|"          HEMCO_Config.rc

#--------------------------------------------------------------------
# Navigate back to source code directory
#--------------------------------------------------------------------
cd ${srcrundir}

#----------------------------------------------------------------------
# Archive repository version in run directory file rundir.version
#----------------------------------------------------------------------
version_log=${rundir}/rundir.version
echo "This run directory was created with ${srcrundir}/createRunDir.sh." > ${version_log}
echo " " >> ${version_log}
echo "HEMCO repository version information:" >> ${version_log}
cd ${hemcodir}
remote_url=$(git config --get remote.origin.url)
code_branch=$(git rev-parse --abbrev-ref HEAD)
last_commit=$(git log -n 1 --pretty=format:"%s")
commit_date=$(git log -n 1 --pretty=format:"%cd")
commit_user=$(git log -n 1 --pretty=format:"%cn")
commit_hash=$(git log -n 1 --pretty=format:"%h")
cd ${srcrundir}
printf "\n  Remote URL: ${remote_url}" >> ${version_log}
printf "\n  Branch: ${code_branch}"    >> ${version_log}
printf "\n  Commit: ${last_commit}"    >> ${version_log}
printf "\n  Date: ${commit_date}"      >> ${version_log}
printf "\n  User: ${commit_user}"      >> ${version_log}
printf "\n  Hash: ${commit_hash}"      >> ${version_log}

#-----------------------------------------------------------------
# Ask user whether to track run directory changes with git
#-----------------------------------------------------------------
printf "${thinline}Do you want to track run directory changes with git? (y/n)${thinline}"
valid_response=0
while [ "$valid_response" -eq 0 ]; do
    read enable_git
    if [[ ${enable_git} = "y" ]]; then
	cd ${rundir}
	printf "\n\nChanges to the following run directory files are tracked by git:\n\n" >> ${version_log}
	git init
	git add *.rc *.sh
	printf " " >> ${version_log}
	git commit -m "Initial run directory" >> ${version_log}
	cd ${curdir}
	valid_response=1
    elif [[ ${enable_git} = "n" ]]; then
	valid_response=1
    else
	printf "Input not recognized. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Done!
#-----------------------------------------------------------------
printf "\nCreated ${rundir}\n"

exit 0
