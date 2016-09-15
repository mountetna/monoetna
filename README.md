# The Metis Repository Software.

## Starting this server

  $ sudo thin -p 80

## Generic Setup for Metis

Create a machine with the user "developer" and password "developer"

### If you are only in development mode it may be handy to allow passwordless sudo
### Add this line near the end of the "/etc/sudoers" file

  `developer ALL = (ALL) NOPASSWD: ALL`

#### update the OS
  
  ```
  $ sudo apt-get update -y
  $ sudo apt-get upgrade -y
  ```

### Install Openssh (if not already installed)
### Install the "developer" ssh pub key (Should be inlcuded in the folder)
  
  `$ sudo apt-get install openssh-server`

### Install other dependences needed for Ruby and Friends

  ```
  $ sudo apt-get install \
    git \
    autoconf \
    bison \
    build-essential \
    libssl-dev \
    libyaml-dev \
    libreadline6-dev \
    zlib1g-dev \
    libncurses5-dev \
    libffi-dev \
    libgdbm3 \
    libgdbm-dev \
    libpq-dev \
    libpcap-dev;
  ```

### Install secondary packages needed for general development

  ```
  $ sudo apt-get install -y \
    dkms \
    linux-headers-generic \
    linux-headers-$(uname -r);
  ```

### Install some system monitoring tools if you need profiling

  ```
  $ sudo apt-get install \
    iftop \
    iotop ;

  $ cd /opt
  $ sudo git clone https://github.com/raboof/nethogs.git
  $ cd ./nethogs
  $ sudo make
  $ cd /usr/local/bin
  $ sudo ln -s /opt/nethogs/src/nethogs ./nethogs
  ```

### If you plan on running rails or sinatra with phusion passenger you need these deps too.
  
  ```
  $ sudo apt-get install \
    apache2 \
    libcurl4-openssl-dev \
    apache2-threaded-dev \
    libapr1-dev \
    libaprutil1-dev;
  ```

### restart the machine to capture any kernal updates

  ```
  $ sudo shutdown -r now
  ```

### Install "rbenv" (Ruby Environment Manager)
  
  https://github.com/rbenv/rbenv
  
  ```
  $ cd ~
  $ git clone https://github.com/rbenv/rbenv.git ~/.rbenv
  ```

  Follow the instructions on the "rbenv" page to finish it's install
  https://github.com/rbenv/rbenv

### Install the "install" untility for "rbenv"
  https://github.com/rbenv/ruby-build###readme

  ```
  $ cd ~
  $ git clone https://github.com/rbenv/ruby-build.git ~/.rbenv/plugins/ruby-build 
  ```

### Install the Ruby version you need. (Rails requires min 2.2.2)

  `$ rbenv install -v 2.2.2`

### Set ruby version 2.2.0 as the system default and that is set correctly

  ```
  $ rbenv global 2.2.2
  $ ruby -v
  ```

### Create the appropriate project folders in "/var"

  ```
  $ cd /var
  $ sudo chown -R developer:developer ./www
  $ cd www
  ```

### Mount the project directory from the host machine if necessary

  `$ sudo mount -t vboxsf -o rw,uid=1000,gid=1000 metis-thin /var/www/metis`
