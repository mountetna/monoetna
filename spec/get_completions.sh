#!/usr/bin/env bash

source /app/*.completion

COMP_LINE=$*
COMP_POINT=${#COMP_LINE}

eval set -- "$@"

COMP_WORDS=("$@")

# add '' to COMP_WORDS if the last character of the command line is a space
[[ ${COMP_LINE[@]: -1} = ' ' ]] && COMP_WORDS+=('')

COMP_CWORD=$(( ${#COMP_WORDS[@]} - 1 ))

completion=$(complete -p "$1" 2>/dev/null | awk '{print $(NF-1)}')

[[ -n $completion ]] || return 1

"$completion"

if [ -z ${COMPREPLY+x} ]; then
  echo "DEFAULT"
else
  printf '%s\n' "${COMPREPLY[@]}" | LC_ALL=C sort
fi
