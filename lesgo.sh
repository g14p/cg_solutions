#!/bin/bash

# Author: Georg Pernice
# Purpose:      A set of bash functions to make easier the comparison of reference and results.
#               Assumes, you are using a terminal workflow with following programs:
#                       - vim
#                       - ycm plugin for vim
#                       - feh to view images
# Date : Looking forward to such things

BASE_PATH=~/git/cg_solutions
cd ${BASE_PATH}/cg_exercise_01/01_colors
SCENES="<no scenes no reference imgs>"
#cd ${BASE_PATH}/cg_exercise_02/02_whitted
#SCENES="box, spheres"
#cd ${BASE_PATH}/cg_exercise_03/03_path_tracing
#SCENES="box, spheres"
#cd ${BASE_PATH}/cg_exercise_04/04_textures
#SCENES="box, spheres"

# Terminal Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
ORANGE='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # Reset (No Color)

echo "If its the first time since cloning u need the rest of the files"
echo "download the .zip from ilias and unzip it here.."
echo -e "${BLUE}"
echo -e " Q U I C K   O P T I O N S"
echo -e " _________________________"
echo
echo "   --- 'v' into the source code to open with nvim"
echo -e "${NC}"
ls assignment_references 2>/dev/null
if [[ $? -eq 0 ]]; then
  echo
  echo "Your call :) "
  echo -e "${BLUE}"
  echo "--- 'showimages' the images"
  echo "--- 'update' to render them again"
  echo "--- 'comparescene' to see the difference in red .. or directly"
else
  echo -e "   ${ORANGE}No images available.${NC}"
fi
echo
echo " O V E R V I E W   O F   P R O G R E S S "
echo
echo " +---------------------------------------------------------"
echo " | "
echo -e " | ${NC}Last Commit:${NC}"
echo -e " | ${GREEN}$(git log -n 1 --decorate=no | tail -n 1)${NC}"
echo " | "
echo -e " | ${NC}Next TODOs${NC}"
echo -e " | ${ORANGE}. . . continue ex02 ${NC}"
echo -e " | ${RED}. . . fix generate grid in ex01${NC}"
echo " | "
echo " +---------------------------------------------------------"

function update {
  ./cg --create-images
}
function showimages {
  if [[ -z $1 ]]; then
    echo "To use this command run it like showimages after_<TEILAUFGABE>"
    echo
    echo Teilaufgaben are a, b, c .
  else
    PROGRESS=$1
    # compare a certain scene at a certain progress -> pipe to diffimage -> view diff image
    feh assignment_references/*${PROGRESS}* &
    feh assignment_images/*${PROGRESS}*
  fi
}

function comparescene {
  if [[ -z $1 ]]; then
    echo "To use this command run it like comparescene <SCENE> <TEILAUFGABE>"
    echo
    echo Scenes are ${SCENES}
    echo Teilaufgaben are a, b, c .
    echo "For your current uebung the available assignment_references are:"
    ls assignment_references
  else
    echo comparing..
    SCENE=$1
    PROGRESS="after_${2}"
    # compare a certain scene at a certain progress -> pipe to diffimage -> view diff image
    mkdir -p assignment_differences
    compare assignment_references/${SCENE}_${PROGRESS}.png assignment_images/${SCENE}_${PROGRESS}.png assignment_differences/${SCENE}_${PROGRESS}.png
    feh assignment_differences/${SCENE}_${PROGRESS}.png
  fi

}
function v {
  nvim -p src/exercise_* ../../lesgo.sh .
}
