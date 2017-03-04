/*
 * This file is used by webpack to package different components into separate
 * end points needed by the application. In this way we can keep all of our JS
 * modular and just drop application components where we need them and add a
 * line here for packaging.
 *
 * This file can also be viewed as describing the different JS sub applications
 * that are used in totality.
 */

module.exports = {

  'entry': {

    'basic-view': './client/js/basic-view-controller.js',
    'network-utils': './client/js/network-utils-controller.js',
    'user-admin': './client/js/user-admin-controller.js'
  },
  'output': {

    'path': './client/js',
    'filename': '[name].bundle.js'
  }
}