var path = require('path');
var MiniCssExtractPlugin = require('mini-css-extract-plugin');
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
      react: path.join(__dirname, 'node_modules/react'),
      'react-dom': path.join(__dirname, 'node_modules/react-dom'),
      'react-redux': path.join(__dirname, 'node_modules/react-redux')
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
        loader: 'babel-loader',

        // Skip any files outside of your project's `src` directory
        include: [
          path.resolve('/app', 'lib/client/jsx'),
          path.resolve('/etna', 'node_modules/etna-js/'),
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
          path.resolve('/app', 'lib/client/scss')
        ],
        use: [MiniCssExtractPlugin.loader, 'css-loader', 'sass-loader']
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
    })
  ]
});
