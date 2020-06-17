## The Front End for the Janus Auth server.

  This is a really simple front end. There should be no interaction (all services are via AJAX) so all of this code is mostly placeholding.

## Static Files

  We do have a few static resrouces that need to be linked.

  ```
  $ sudo -i -u [USER] ln -s /var/www/janus-static/img /var/www/janus/client/img
  $ sudo -i -u [USER] ln -s /var/www/janus-static/fonts /var/www/janus/client/fonts
  ```
