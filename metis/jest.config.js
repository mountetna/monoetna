module.exports = {
  globals: {
    fetch: require('node-fetch'),
    CONFIG: {
      project_name: 'labors',
      polyphemus_host: 'https://polyphemus.test'
    }
  },
  transformIgnorePatterns: [
    // "node_modules/(?!(etna-js)/)"
  ],
  moduleNameMapper: {
    '^service-worker-loader!': '<rootDir>/__mocks__/service-worker-loader.js',
    '^@material-ui/core/styles':
      '<rootDir>/node_modules/@material-ui/core/styles',
    '^color$': '<rootDir>/node_modules/color',
    '^color-string$': '<rootDir>/node_modules/color-string'
  },
  testMatch: ['**/__tests__/**/?(*.)(spec|test).(j|t)s?(x)'],
  collectCoverageFrom: ['**/*.js?(x)'],
  setupFilesAfterEnv: ['./spec/setup.js'],
  setupFiles: ['raf/polyfill']
};
