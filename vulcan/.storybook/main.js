const baseWebpackConfig = require('../webpack.config')(process.env);

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
    return {
      ...config,
      module: {
        ...config.module,
        rules: [...config.module.rules, ...baseWebpackConfig.module.rules],
      },
      plugins: [...config.plugins, ...baseWebpackConfig.plugins],
      resolve: {
        ...config.resolve,
        symlinks: baseWebpackConfig.resolve.symlinks,
        extensions: [...config.resolve.extensions, ...baseWebpackConfig.resolve.extensions],
        alias: {...config.resolve.alias, ...baseWebpackConfig.resolve.alias},
      }
    };
  },
}