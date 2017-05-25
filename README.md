# The Metis File Repository Software.

## Why?

Here at the Computational Biology Core there is a need to store large files from labratory experiements. These large files could be kept on Amazon S3 but that is undesirable for two reasons.

1. It cost a lot of money to get our files back if we use a service like Amazon S3.
2. The files need to be near our HPC Cluster for processing and analysis.

## What?

This file server software allows us to upload and download files using HTTP, much in the same way Box does. The server also allows us to do these operations via AJAX calls. To ensure security the server generates HMAC tokens that are used to verify the user and data being transfered. The server also has a small front end to view the files that were uploaded.

## How?

We deploy this server on an Linux system. The application itself is a Ruby/Rack application. The web server is Apache using Phusion Passenger. The file system for saving the actual files is a large disk that we mount on the linux file system. We keep metadata about the files on disk in a Postgres DB.

## Other info?

This system has a layer of user authentication. This file server is part of a larger system called Mt. Etna and this file server relys upon our other project, called Janus, for user authentication. The UI is built with Bable JS and Webpack.

## Depenencies?

  System Packages:
```
  'git', 'tmux', 'openssl', 'curl', 'python-pip', 'httpd', 'httpd-devel', 'mod_ssl', 'openssl-devel', 'readline-devel', 'zlib-devel', 'postgresql-devel', 'nodejs', 'postgresql-server', 'postgresql-contrib', 'pygpgme'
```

  Ruby Gems: 
```
  'rack', 'pg', 'sequel', 'bundler'
```

  Node NPM:
```
  'babel-cli', 'webpack'
```

## Local Development?

I run VMs that look like the STAGE and PROD machines. I then mount my local project folder onto that machine. Here is the mount command...

  `$ sudo mount -t vboxsf -o rw,uid=1000,gid=1000 metis /var/www/metis`

When you first start up a VM for local development we need to enable symlinks.

  `$ VBoxManage setextradata metis-dev VBoxInternal2/ SharedFoldersEnableSymlinksCreate/metis 1`
