#!/bin/bash

# Author: Georg Pernice
# Purpose: A set of bash functions to make easier the comparison of reference and results
# Date : Looking forward to such things

echo "Your call :) "
echo "'showimages' the images" 
echo "or 'update' to render them again"
echo "or vim into the source code with 'v'"
echo compare_spheres to see the difference of spheres

cd ~/git/cg_exercises/cg_exercise_02/02_whitted
PROGRESS="after_b"
SCENE=spheres
SCENE=box

function update {
    ./cg --create-images ; 
}
function showimages {
    feh assignment_references/*${PROGRESS}* & feh assignment_images/*${PROGRESS}* ;
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
