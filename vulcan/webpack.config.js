var path = require('path');
var MiniCssExtractPlugin = require('mini-css-extract-plugin');
var webpack = require('webpack');

module.exports = (env) => ({
  mode: env.NODE_ENV || 'development',
  context: path.resolve(__dirname),
  resolve: {
    extensions: ['.js', '.jsx', '.ts', '.tsx', '.scss', '.png', '.jpg', '.jpeg', '.svg', '.css'],
    alias: {
      'code-mirror': path.join(__dirname, 'node_modules/codemirror/lib'),
      react: path.join(__dirname, 'node_modules/react'),
      'react-dom': path.join(__dirname, 'node_modules/react-dom'),
      'react-redux': path.join(__dirname, 'node_modules/react-redux')
    },
    symlinks: false
  },
  entry: ['./lib/client/jsx/vulcan.jsx', './lib/client/scss/application.scss'],
  output: {
    filename: 'public/js/vulcan.bundle.js',
    path: __dirname,
    publicPath: '/js/'
  },
  module: {
    rules: [
      {
        loader: 'babel-loader',

        // Skip any files outside of your project's `src` directory
        include: [
          path.resolve(__dirname, 'lib/client/jsx'),
          path.resolve(__dirname, 'node_modules/etna-js/'),
          path.resolve(__dirname, 'stories/'),
          '/etna/packages/etna-js'
        ],

        // Only run `.js`, `.jsx`, `.ts`, and `.tsx` files through Babel
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
        test: /\.css$/,
      },

      {
        loader: 'file-loader',
        include: [
          path.resolve(__dirname, 'node_modules/etna-js/'),
          '/etna/packages/etna-js',
          path.resolve(__dirname, 'lib/client/img/')
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
            __dirname,
            'node_modules/react-loader-spinner/dist/loader/css'
          ),
          path.resolve(__dirname, 'node_modules/etna-js/'),
          '/etna/packages/etna-js',
          path.resolve(__dirname, 'node_modules/react-notifications-component'),
          path.resolve(__dirname, 'lib/client/scss')
        ],

        // loader: ExtractTextPlugin.extract(['css-loader', 'sass-loader']),
        // loader: env.NODE_ENV === 'storybook' ? ['css-loader', 'sass-loader'] : undefined,
        use: [MiniCssExtractPlugin.loader, 'css-loader', 'sass-loader'],
      }
    ]
  },
  plugins: [
    // new ExtractTextPlugin({
    //   // define where to save the file
    //   filename: 'public/css/vulcan.bundle.css',
    //   allChunks: true
    // }),
    new MiniCssExtractPlugin({
      filename: 'public/css/vulcan.bundle.css',
    }),
    new webpack.DefinePlugin({
      'process.env': {
        NODE_ENV: JSON.stringify(env ? env.NODE_ENV : 'development')
      }
    })
  ]
});