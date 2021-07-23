var webpack = require('webpack');
var path = require('path');
const baseWebpackConfig = require('../webpack.config')({...process.env, NODE_ENV: 'storybook'});

module.exports = {
  stories: [
    "../stories/**/*.stories.mdx",
    "../stories/**/*.stories.@(js|jsx|ts|tsx)"
  ],
  addons: [
    "@storybook/addon-links",
    "@storybook/addon-essentials"
  ],
  webpackFinal: (config) => {
    console.log(config.module.rules);
    return {
      ...config,
      module: {
        ...config.module,
        rules: [
          {
            test: /\.(stories|story)\.[tj]sx?$/,
            loader: '/app/node_modules/@storybook/source-loader/dist/cjs/index.js',
            options: { injectStoryParameters: true, inspectLocalDependencies: true },
            enforce: 'pre'
          },

          {
            test: /\.(svg|ico|jpg|jpeg|png|apng|gif|eot|otf|webp|ttf|woff|woff2|cur|ani|pdf)(\?.*)?$/,
            loader: '/app/node_modules/@storybook/builder-webpack4/node_modules/file-loader/dist/cjs.js',
            options: { esModule: false, name: 'static/media/[path][name].[ext]' }
          },

          ...baseWebpackConfig.module.rules,
          // ...config.module.rules
        ],
      },
      plugins: [
        ...config.plugins,
        ...baseWebpackConfig.plugins,
      ],
      resolve: {
        ...config.resolve,
        symlinks: baseWebpackConfig.resolve.symlinks,
        extensions: [...config.resolve.extensions, ...baseWebpackConfig.resolve.extensions],
        alias: {...config.resolve.alias, ...baseWebpackConfig.resolve.alias},
      }
    };
  },
}