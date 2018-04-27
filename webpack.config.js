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
var ExtractTextPlugin = require('extract-text-webpack-plugin');

module.exports = {
  context: path.resolve(__dirname),
  resolve: {
    extensions: [ '.js', '.jsx' ],
    alias: {
      'font-awesome': path.join(__dirname, 'node_modules/font-awesome')
    }
  },
  entry: {
    'metis-main': './lib/client/jsx/metis.jsx',
    'utils': './lib/client/jsx/libraries/utils.jsx',
    'metis-stylesheet': './lib/client/css/metis.scss'
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
      },

      {
        test: /\.(png|jpg|jpeg|gif|svg|woff|woff2|ttf|eot)(\?.*$|$)/,

        include: [
          path.resolve(__dirname, 'node_modules/font-awesome'),
        ],

        loader: 'file-loader',

        options: {
          name: '[name].[ext]',
          outputPath: 'public/fonts/',
          publicPath: function(url) { return url.replace(/public/,'') }
        }
      },

      {
        // sass / scss loader for webpack
        test: /\.(sass|scss)$/,

        loader: ExtractTextPlugin.extract(['css-loader', 'sass-loader'])
      }
    ]
  },
  plugins: [
    new ExtractTextPlugin({ // define where to save the file
      filename: 'public/css/metis.bundle.css',
      allChunks: true,
    }),
  ]
}
