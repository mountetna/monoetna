#!/usr/bin/env bash

set -ex

docker run --rm -v $(realpath $1):/work/input.yml:ro bgruber/yaml2json input.yml | docker run --rm -i imega/jq -c --stream  | python -c "
import json, sys

for line in sys.stdin.readlines():
  obj = json.loads(line)
  if len(obj) > 1:
    path, value = obj
    print('eval \"\${APP_NAME}__' + '__'.join([p[1:].upper() for p in path[1:]]) + '=' + str(value) + '\"')
"
