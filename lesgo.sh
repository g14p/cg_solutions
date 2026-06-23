#!/bin/bash

# Author: Georg Pernice
# Purpose:      A set of bash functions to make easier the comparison of reference and results.
#               Assumes, you are using a terminal workflow with following programs:
#                       - vim
#                       - ycm plugin for vim
#                       - feh to view images
# Date : Looking forward to such things

#cd ~/git/cg_exercises/cg_exercise_01/01_colors
#cd ~/git/cg_exercises/cg_exercise_02/02_whitted
#cd ~/git/cg_exercises/cg_exercise_03/03_path_tracing
cd ~/git/cg_exercises/cg_exercise_04/04_textures

ls assignment_references > /dev/null 
if [[ $? -eq 0 ]]; 
    then echo 
    echo "Your call :) "
    echo "'showimages' the images" 
    echo "or 'update' to render them again"
    echo "'comparescene' to see the difference in red .. or directly"
    echo "run 'init-ycm' if its the first time to initialize autocompletion for this exercise";
    else echo No images available.; 
fi
echo "vim into the source code with 'v'"
echo
echo Last thing i did: 
echo --------------------
echo $(git log -n 1 --decorate=no | tail -n 1)
echo Some time may need to 
echo .  .  . fix those artefacts by perfecting kugelschnitt
echo .  .  . fix the pixelerrors in ex03 as their rand actuallz uses seed? but unsure if possible
echo .  .  . fix in ex03 at least the cornellbox -- the reference has more wild pixels -- so update 1h and see if gets better

function init-ycm {
    echo Initialize Ycm by generating compile_commands.json
    sed -i "s/project (cg)/project(cg)\nset( CMAKE_EXPORT_COMPILE_COMMANDS ON ) # ONLY FOR YouCompleteMe/" CMakeLists.txt 
    cmake .

    
}

function update {
    ./cg --create-images ; 
}
function showimages {
    if [[ -z $1 ]]; 
    then echo "To use this command run it like showimages after_<TEILAUFGABE>"
        echo  
        echo Teilaufgaben are a, b, c . ;
    else
        PROGRESS=$1
        # compare a certain scene at a certain progress -> pipe to diffimage -> view diff image
        feh assignment_references/*${PROGRESS}* & feh assignment_images/*${PROGRESS}* ;
    fi
}

function comparescene {
    if [[ -z $1 ]]; 
    then echo "To use this command run it like comparescene <SCENE> after_<TEILAUFGABE>"
        echo  
        echo Scenes are box, spheres. 
        echo Teilaufgaben are a, b, c . ;
    else
        SCENE=$1
        PROGRESS=$2
        # compare a certain scene at a certain progress -> pipe to diffimage -> view diff image
        mkdir -p assignment_differences
        compare assignment_references/${SCENE}_${PROGRESS}.png assignment_images/${SCENE}_${PROGRESS}.png assignment_differences/${SCENE}_${PROGRESS}.png;
        feh assignment_differences/${SCENE}_${PROGRESS}.png ;
     fi

}
function v {
    vim -p src/exercise_* ~/lesgo.sh .  ;
}


