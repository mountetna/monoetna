#!/bin/sh

set -e

if ! test -z "$TARGET_ROOT_BRANCH"; then
  echo "found branch root at $TARGET_ROOT_BRANCH" >&2
  echo $TARGET_ROOT_BRANCH
  exit 0
fi

echo "found branch root at HEAD~1" >&2
echo HEAD~1
