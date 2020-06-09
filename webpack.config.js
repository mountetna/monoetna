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
const MiniCssExtractPlugin = require('mini-css-extract-plugin');

module.exports = {
  context: path.resolve(__dirname),
  resolve: {
    extensions: [ '.js', '.jsx' ],
    alias: {
      'font-awesome': path.join(__dirname, 'node_modules/@fortawesome/fontawesome-free')
    }
  },
  entry: {
    'main': './src/index.js',
    'stylesheet': './src/css/etna.scss'
  },
  output: {
    path: __dirname,
    filename: './lib/dist/[name].js'
  },
  module: {
    rules: [
      {
        loader: 'babel-loader',
        include: [ path.resolve(__dirname, 'src/jsx'), ],
        test: /\.jsx?$/,
        query: {
          presets: ['env', 'stage-0', 'react'],
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
          publicPath: function(url) { return url.match(/^public/) ? url.replace(/public/,'') : `/images/${url}` }
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
          publicPath: function(url) { return url.match(/^public/) ? url.replace(/public/,'') : `/fonts/${url}` }
        }
      },
      {
        test: /\.s[ac]ss$/i,
        include: [ path.resolve(__dirname, 'src/css'), ],
        use: [
          'style-loader',
          'css-loader',
          {
            loader: 'sass-loader',
            options: {
              implementation: require('node-sass'),
            },
          },
        ],
      },
      {
        test: /\.css$/,
        include: [ path.resolve(__dirname, 'src/css'), ],
        use: [MiniCssExtractPlugin.loader, 'css-loader'],
      }
    ]
  },
  plugins: [
    new MiniCssExtractPlugin({ // define where to save the file
      filename: './lib/[name].css'
    }),
  ],
  optimization: {
    splitChunks: {
      cacheGroups: {
        styles: {
          name: 'styles',
          test: /\.css$/,
          chunks: 'all',
          enforce: true,
        },
      },
    },
  },
}
