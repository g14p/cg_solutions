#!/bin/bash

# Author: Georg Pernice
# Purpose: A set of bash functions to make easier the comparison of reference and results
# Date : Looking forward to such things
cd ~/git/cg_exercises/cg_exercise_02/02_whitted
    #cd ~/git/cg_exercises/cg_exercise_03/03_path_tracing

echo "Your call :) "
echo "'showimages' the images" 
echo "or 'update' to render them again"
echo "or vim into the source code with 'v'"
echo 'comparescene' to see the difference in red
echo
echo Last thing i did: 
echo --------------------
echo $(git log -n 1 --decorate=no | tail -n 1)
echo Next time need to fix those artefacts appearing since b. I think they come from implementation of b.

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
    vim -p src/exercise_02.cpp ~/lesgo.sh .  ;
}
