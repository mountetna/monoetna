#!/bin/sh

set -e

git log --pretty='format:%h %an' | grep GithubAction | head -n 1 | cut -d' ' -f1

