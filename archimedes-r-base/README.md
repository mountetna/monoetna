# archimedes-r-base
Build base docker image for archimedes-r, and manage pushes to dockerhub

# Adding new packages
1) Add target packages to the `Imports:` section of the DESCRIPTION file, using `Remotes:` as well for packages not in CRAN.
2) Commit / PR to master. Then the triggered build attempt, the jetpack::install() call will pull any new packages and update the renv.lock file.
3) Run `bash extract_package_files.sh` from this repo's base folder to pull the new renv.lock file.  Then commit this in too.
