module.exports = {
  presets: [
    [
      '@babel/preset-env',
      {
        useBuiltIns: 'entry',
        corejs: 3
      }
    ],
    '@babel/react',
    '@babel/typescript'
  ],
  plugins: [
    '@babel/plugin-proposal-class-properties',
    '@emotion',
    '@babel/transform-runtime',
    ['@babel/plugin-proposal-private-methods', {loose: false}],
    ['@babel/plugin-proposal-private-property-in-object', {loose: false}]
  ]
};
