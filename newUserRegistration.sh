#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#                    Mofidied for WRF-GC by Ao Ding                           ! 
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: newUserRegistration.sh
#
# !DESCRIPTION: Defines utility functions for first-time GEOS-Chem
#  user registration.  This code has been been abstracted out of
#  run/GCClassic/createRunDir.sh and run/GCHP/createRundir.sh.
#\\
#\\
# !REVISION HISTORY:
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

####=========================================================================
#### Common variables
####==========================================================================

# User prompt for reading in data
USER_PROMPT=">>> "

####=========================================================================
#### Send registration details to atmoschem wordpress website via curl
####==========================================================================
function postRegistration() {
    #### Sends registration details to API ####
    curl --location --request POST "https://www.download.atmoschem.org.cn/wp-json/contact-form-7/v1/contact-forms/20/feedback" \
        --form 'your-name='"${2}" \
        --form 'your-email='"${1}" \
        --form 'institution='"${3}" \
        --form '_wpcf7_unit_tag="wpcf7-f20-p47-o1"' \
        --form 'name_of_pi='"${4}" \
        --form 'site='"${5}" \
        --form 'git_username='"${6}" \
        --form 'research_interest='"${7}" \
        --form 'env_type='"${8}" 
}

####=========================================================================
#### Read an entry from the command line
####==========================================================================
function userInput() {

    # Keep asking for the input until user enters a non-blank input
    while [[ -z "${val}" ]]; do
        IFS='\n' read -r -p "${USER_PROMPT}" val
    done

    # Return the answer
    printf "${val}"
}

####=========================================================================
#### Query user to provide registration information
####
#### NOTE: See https://www.baeldung.com/linux/ifs-shell-variable
#### for a description of what the IFS variable below does
####==========================================================================
function registerNewUserwrapper() {
    # Check if user is registered
    if [[ -f ${HOME}/.wrfgc/config ]]; then 
        source ${HOME}/.wrfgc/config
    else 
        printf " @echo \nThis will be stored in ${HOME}/.wrfgc/config for future automatic use.\n" ;
        mkdir -p ${HOME}/.wrfgc; 
    fi
    # Check if user is registered

    if [[ -z "${GC_USER_REGISTERED}" ]];then 
        registerNewUser
        fi
}
function registerNewUser() {
    

    # Ask user several questions
    printf "\nInitiating User Registration:\n"
    printf "You will only need to fill this information out once.\n"
    printf "Please respond to all questions.\n"

    printf "${thinline}What is your name?${thinline}"
    name=$(userInput)

    # Ask for email
    printf "${thinline}What is your email address?${thinline}"
    email=$(userInput)

    printf "${thinline}What is the name of your research institution?${thinline}"
    institution=$(userInput)

    printf "${thinline}What is the name of your principal invesigator?\n"
    printf "(Enter 'self' if you are the principal investigator.)${thinline}"
    name_of_pi=$(userInput)

    printf "${thinline}Please provide the web site for your institution\n"
    printf "(group website, company website, etc.)?${thinline}"
    site=$(userInput)

    printf "${thinline}Please provide your github username (if any) so that we\n"
    printf "can recognize you in submitted issues and pull requests.${thinline}"
    git_username=$(userInput)

    printf "${thinline}Where do you plan to run WRF-GC?\n"
    printf "(e.g. local compute cluster, AWS, other supercomputer)?${thinline}"
    env_type=$(userInput)

    printf "${thinline}Please briefly describe how you plan on using WRF-GC\n"
    printf "${thinline}"
    research_interest=$(userInput)

    # Send information to database on AWS
    postRegistration "${email}"             "${name}"       "${institution}"  \
                     "${name_of_pi}"        "${site}"       "${git_username}" \
                     "${research_interest}"  "${env_type}"

    # Update the .wrfgc/config file and apply settings
    echo "export GC_USER_REGISTERED=true" >> "${HOME}/.wrfgc/config"
    . ${HOME}/.wrfgc/config
}


