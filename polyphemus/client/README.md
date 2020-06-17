## The Front End for the Polyphemus Network Manager

## Dependencies, unfortunately we have dependencies.

### These items should be installed globablly as command line utilties.
  * nodejs
  * npm
    * babel - this will convert our ES6 code to ES5, but in the future we can omit this. This will also convert our jsx code into js. So we do not need a jsx util.
    * webpack - this will minify and condense our JS code
    * jest - javascript unit testing

  ```
  $ sudo apt-get install nodejs;
  $ sudo apt-get install npm;
  $ sudo npm install -g babel-cli;
  $ sudo npm install -g webpack;
  $ sudo npm install -g jest;
  ```

### Symlinking Node
  
  Sometimes the package manager for the system will have 'nodejs' rather than 'node' installed. In this case you want to symlink 'nodejs' to 'node' so babel and webpack can use it properly.

  `$ sudo ln -s `which nodejs` /usr/bin/node`
 
### These items should be installed within the project folder.

  ```
  $ npm install --save react
  $ npm install --save react-dom
  $ npm install --save react-redux
  $ npm install --save redux
  $ npm install --save babel-preset-es2015
  $ npm install --save babel-preset-react
  ```

  OR you can just run...

  `$ npm install`

  npm will read the 'package.json' file and install the appropriate dependencies.

### A note on the Bootstrap Lib.
  
  We are only using the glyphicons from bootstrap. So don't bother using any other Bootstrap graphical elements.

## Active Development.

  There are two utilities that you will need to run while you are developing.

  * babel
  * webpack

  I like running two terminals (using tmux) side by side and running each util in a separate pane so I can see the build process as it happens

  First start up babel to watch the 'jsx' files and build normal 'js' files. In my case even regular Javascript code is labeled with the .jsx extension so I can tell a processed from unprocessed file.

  `$ babel ./client/jsx --watch --out-dir ./client/js`

  If the `webpack.config.js` file is present in the parent folder you can just run:

  `$ webpack --watch`

## NPM scripted tasks. 

  In the 'package.json' file there are 'tasks' that can be run, which are just short cuts for commands.
  
  ```
  $ npm run clean; // delete the contents of the './js' folder
  $ npm run build; // run babel and webpack one time to 'build' the UI
  $ npm run babel; // run babel on the './jsx' folder and output to the './js' folder
  $ npm run webpack; // run webpack on the entry script and produce a packaged JS app 
  ```

  Again you can look inside the 'package.json' file to see the details of these operations.

## Static Files

  Since we don't keep our images and fonts in the repo we just create symlinks to the static resources.

  ```
  $ sudo -i -u [USER] ln -s /var/www/polyphemus-static/img /var/www/polyphemus/client/img
  $ sudo -i -u [USER] ln -s /var/www/polyphemus-static/fonts /var/www/polyphemus/client/fonts
  ```
