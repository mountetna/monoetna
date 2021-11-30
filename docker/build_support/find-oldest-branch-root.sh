#!/bin/sh

set -e

if [[ -z "$IS_CI" ]]; then
  echo "found branch root at HEAD~1" >&2
  echo HEAD~1
  exit 0
fi

if ! [[ -z "$TARGET_ROOT_BRANCH "]]; then
  echo "found branch root at $TARGET_ROOT_BRANCH" >&2
  echo $TARGET_ROOT_BRANCH
  exit 0
fi

result=$(git log --pretty='format:%h %an' | grep GithubAction | head -n 2 | cut -d' ' -f1)
echo "found branch root at $result" >&2
echo "$result"


