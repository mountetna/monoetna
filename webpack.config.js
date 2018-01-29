/*
 * This file is used by webpack to package different components into separate
 * end points needed by the application. In this way we can keep all of our JS
 * modular and just drop application components where we need them and add a
 * line here for packaging.
 *
 * This file can also be viewed as describing the different JS sub applications
 * that are used in totality.
 */

var path = require('path');

module.exports = {
  context: path.resolve(__dirname),
  resolve: {
    extensions: [ '.js', '.jsx' ]
  },
  entry: {
    'metis-main': './lib/client/jsx/metis-uploader-controller.jsx',
    'utils': './lib/client/jsx/libraries/utils.jsx'
  },
  output: {
    path: __dirname,
    filename: 'public/js/[name].bundle.js'
  },
  module: {
    rules: [
      {
        loader: 'babel-loader',
        include: [ path.resolve(__dirname, 'lib/client/jsx'), ],
        test: /\.jsx?$/,
        query: {
          presets: ['es2015', 'stage-0', 'react'],
        }
      },

      {
        loader: 'file-loader',
        test: /\.(jpe?g|png|svg)$/i,
        include: [
          path.resolve(__dirname, 'lib/client/images'),
        ],
        options: {
          name: '[name].[ext]',
          outputPath: 'public/images/',
          publicPath: function(url) { return url.replace(/public/,'') }
        }
      }
    ]
  }
}
