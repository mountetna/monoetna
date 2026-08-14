#!/bin/bash

CONFIG_YAML=$1
# get groups from yaml
readarray -t groups < <(yq ".\":$METIS_ENV\".\":storage\" | keys[] | sub(\":\",\"\")" $CONFIG_YAML)

# for each group in yaml, add the user
for group in "${groups[@]}"
do
	GID=$(yq ".\":$METIS_ENV\".\":storage\".\":$group\".\":gid\"" $CONFIG_YAML)
	addgroup $group --gid $GID
	usermod -a -G $group daemon
done
