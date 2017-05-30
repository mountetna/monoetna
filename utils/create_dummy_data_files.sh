#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
  cp "../source/small_file_x" "/data1/Ipi/$line"
done < "$1"
