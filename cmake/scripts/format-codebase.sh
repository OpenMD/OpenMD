#!/bin/bash
#
# This script was originally taken from the format-codebase.sh file in the
# Drychem libraries. The code has been modified to match the OpenMD coding
# style.
#
# We have retained the DryChem copyright and MIT license on this file:
#
# Copyright (c) 2020-2021 Cody R. Drisko. All rights reserved.
# Licensed under the MIT License. See the LICENSE file in the project root for more information.
#
# Name: format-codebase.sh - Version 1.2.0
# Author: crdrisko
# Date: 10/07/2020-07:20:58
# Description: Use clang-format on all files in the repository with the option to ignore specified files


#@ DESCRIPTION: Print the format-code program's help message
#@ USAGE: printHelpMessage
printHelpMessage() {
  printf "\nUSAGE: format-codebase [-h] [-f path] [-i fileName]\n\n"
  printf "  -h  Prints help information about the format-codebase program.\n\n"
  printf "  -f  REQUIRED: Absolute path to clang-format.\n"
  printf "  -i  OPTIONAL: Filename and path (relative to project root) to ignore\n"
  printf "        when formatting.\n\n"
  printf "EXAMPLE: format-codebase -f /opt/local/libexec/llvm-10/bin/clang-format -i \"src/antlr/*\"\n\n"
}

#@ DESCRIPTION: Use clang-format to format each file in the repository
#@ USAGE: formatFiles LIST
formatFiles() {
  for elem in "$@"; do
    if [[ -f "$elem" && ("${elem##*.}" == cpp || "${elem##*.}" == hpp) ]]; then
      for file in "${ignoreFiles[@]}"; do
        if [[ "$PWD/$elem" == "$file" ]]; then
          continue 2;
        fi
      done

      printf "Formatting: %s\n" "$elem"
      "${formatterPath:?Path to clang-format is required.}" -i -style=file "$elem"
    elif [[ -d $elem ]]; then
      cd "$elem" || ( printf "Could not change into required directory.\n" && exit 1 )

      formatFiles -- *

      cd ../
    fi
  done
}

#@ DESCRIPTION: Execute the main portion of the script's code
#@ USAGE: main LIST
main() {
  declare -a ignoreFiles
  declare formatterPath

  local OPTIND opt

  while getopts f:i:h opt; do
    case $opt in
      f) formatterPath="${OPTARG}" ;;
      i) if [ "${OPTARG:$(( ${#OPTARG} - 1 )):1}" == '*' ]; then
           for file in $PWD/$OPTARG; do
             ignoreFiles+=( "$file" )
           done
         else
           ignoreFiles+=( "$PWD/$OPTARG" )
         fi ;;
      h) printHelpMessage && exit 0 ;;
      *) printf "Invalid option flag passed to program.\n" && exit 1 ;;
    esac
  done

  shift $((OPTIND-1))

  formatFiles -- *
}

main "$@"
