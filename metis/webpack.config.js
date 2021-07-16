var path = require('path');
var ExtractTextPlugin = require('extract-text-webpack-plugin');
const EventHooksPlugin = require('event-hooks-webpack-plugin');
const fs = require('fs-extra');

module.exports = {
  context: path.resolve(__dirname),
  resolve: {
    extensions: ['.js', '.jsx', '.ts', '.tsx'],
    alias: {
      'font-awesome': path.join(
        __dirname,
        'node_modules/@fortawesome/fontawesome-free'
      ),
      react: path.join(__dirname, 'node_modules/react'),
      'react-dom': path.join(__dirname, 'node_modules/react-dom'),
      'react-redux': path.join(__dirname, 'node_modules/react-redux')
    },
    symlinks: false
  },
  entry: {
    'metis-main': './lib/client/jsx/metis.jsx',
    'metis-stylesheet': './lib/client/css/metis.scss'
  },
  output: {
    path: __dirname,
    filename: 'public/js/[name].bundle.js',
    publicPath: '/js/'
  },
  module: {
    rules: [
      {
        loader: 'babel-loader',
        include: [
          path.resolve(__dirname, 'node_modules/etna-js/'),
          path.resolve(__dirname, 'node_modules/downzip/'),
          path.resolve(__dirname, 'lib/client/jsx')
        ],
        test: /\.(js|ts)x?$/
      },

      {
        loader: ['style-loader', 'css-loader'],
        include: [
          path.resolve(__dirname, 'node_modules/etna-js/'),
          path.resolve(__dirname, 'node_modules/animate.css/'),
          path.resolve(__dirname, 'node_modules/react-notifications-component'),
          '/etna/packages/etna-js'
        ],
        test: /\.css$/
      },

      {
        loader: 'file-loader',
        include: [
          path.resolve(__dirname, 'node_modules/etna-js/'),
          '/etna/packages/etna-js'
        ],
        test: /\.(jpe?g|png|svg)$/i,

        options: {
          name: '[name].[ext]',
          outputPath: 'public/images/',
          publicPath: '/images'
        }
      },

      {
        loader: 'file-loader',
        test: /\.(jpe?g|png|svg)$/i,
        include: [path.resolve(__dirname, 'lib/client/images')],
        options: {
          name: '[name].[ext]',
          outputPath: 'public/images/',
          publicPath: function (url) {
            return url.match(/^public/)
              ? url.replace(/public/, '')
              : `/images/${url}`;
          }
        }
      },

      {
        test: /\.(png|jpg|jpeg|gif|svg|woff|woff2|ttf|eot)(\?.*$|$)/,

        include: [
          path.resolve(__dirname, 'node_modules/@fortawesome/fontawesome-free')
        ],

        loader: 'file-loader',

        options: {
          name: '[name].[ext]',
          outputPath: 'public/fonts/',
          publicPath: function (url) {
            return url.match(/^public/)
              ? url.replace(/public/, '')
              : `/fonts/${url}`;
          }
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
    new ExtractTextPlugin({
      // define where to save the file
      filename: 'public/css/metis.bundle.css',
      allChunks: true
    }),
    new EventHooksPlugin({
      'after-emit': (compilation, done) => {
        console.log('\n\nCopying source files to compiled\n\n');
        fs.copy('downzip-sw.js', 'public/js/downzip-sw.js', done);
      }
    })
  ]
};
