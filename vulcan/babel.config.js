module.exports = {
  presets: ['@babel/env', '@babel/react', '@babel/typescript'],
  plugins: ['@babel/plugin-proposal-class-properties', 'emotion', ['@babel/plugin-proposal-private-methods', { 'loose': false }], ["@babel/plugin-proposal-private-property-in-object", { "loose": false }]]
};
