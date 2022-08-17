var path = require('path');
var MiniCssExtractPlugin = require('mini-css-extract-plugin');
const EventHooksPlugin = require('event-hooks-webpack-plugin');
const fs = require('fs-extra');
var webpack = require('webpack');

module.exports = (env) => ({
  mode: process.env.NODE_ENV || 'development',
  context: '/app',
  devtool: 'source-map',
  resolve: {
    modules: ['/etna/node_modules'],
    extensions: [
      '.js',
      '.jsx',
      '.ts',
      '.tsx',
      '.scss',
      '.png',
      '.jpg',
      '.jpeg',
      '.svg',
      '.css'
    ],
    alias: {
      'code-mirror': path.join(__dirname, 'node_modules/codemirror/lib'),
      'font-awesome': path.join(
        __dirname,
        'node_modules/@fortawesome/fontawesome-free'
      ),
      react: path.join(__dirname, 'node_modules/react'),
      'react-dom': path.join(__dirname, 'node_modules/react-dom'),
      'react-redux': path.join(__dirname, 'node_modules/react-redux'),
      stream: path.join(__dirname, 'node_modules/stream-browserify')
    },
    symlinks: false
  },
  entry: [
    `./lib/client/jsx/${env.APP_NAME}.jsx`,
    './lib/client/scss/application.scss'
  ],
  output: {
    filename: `public/js/${env.APP_NAME}.bundle.js`,
    path: '/app',
    publicPath: '/js/'
  },
  module: {
    rules: [
      {
        test: /\.js$/,
        loader: require.resolve('@open-wc/webpack-import-meta-loader')
      },
      {
        loader: 'babel-loader',

        // Skip any files outside of your project's `src` directory
        include: [
          path.resolve('/app', 'lib/client/jsx'),
          path.resolve('/etna', 'node_modules/etna-js/'),
          path.resolve('/etna', 'node_modules/downzip/'),
          path.resolve('/app', 'stories/'),
          '/etna/packages/etna-js'
        ],

        // Only run `.js`, `.jsx`, `.ts`, and `.tsx` files through Babel
        test: /\.(js|ts)x?$/
      },

      {
        loader: ['style-loader', 'css-loader'],
        include: [
          path.resolve('/etna', 'node_modules/etna-js/'),
          path.resolve('/etna', 'node_modules/animate.css/'),
          path.resolve('/etna', 'node_modules/react-notifications-component'),
          '/etna/packages/etna-js'
        ],
        test: /\.css$/
      },

      {
        loader: 'file-loader',
        include: [
          path.resolve('/etna', 'node_modules/etna-js/'),
          '/etna/packages/etna-js',
          path.resolve('/app', 'lib/client/img/')
        ],
        test: /\.(jpe?g|png|svg)$/i,

        options: {
          name: '[name].[ext]',
          outputPath: 'public/images/',
          publicPath: '/images'
        }
      },

      {
        // sass / scss loader for webpack
        test: /\.(sass|scss)$/,
        include: [
          path.resolve(
            '/app',
            'node_modules/react-loader-spinner/dist/loader/css'
          ),
          path.resolve(__dirname, 'node_modules/etna-js/'),
          '/etna/packages/etna-js',
          path.resolve(__dirname, 'node_modules/react-notifications-component'),
          path.resolve(__dirname, 'node_modules/codemirror/'),
          path.resolve('/app', 'lib/client/scss')
        ],
        use: [MiniCssExtractPlugin.loader, 'css-loader', 'sass-loader']
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
      }
    ]
  },
  plugins: [
    new MiniCssExtractPlugin({
      filename: `public/css/${env.APP_NAME}.bundle.css`
    }),
    new webpack.DefinePlugin({
      'process.env': {
        NODE_ENV: JSON.stringify(process.env.NODE_ENV || 'development')
      }
    }),
    new EventHooksPlugin({
      'after-emit': (compilation, done) => {
        if ('metis' === env.APP_NAME) {
          console.log('\n\nCopying source files to compiled\n\n');
          fs.copy('downzip-sw.js', 'public/js/downzip-sw.js', done);
        }
      }
    })
  ]
});
