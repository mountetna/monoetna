module.exports = {
  globals: {
    fetch: require('node-fetch'),
    CONFIG: {
      project_name: 'labors'
    }
  },
  transformIgnorePatterns: [
    // "node_modules/(?!(etna-js)/)"
  ],
  moduleNameMapper: {
    "^service-worker-loader\!": "<rootDir>/__mocks__/service-worker-loader.js"
  },
  testMatch: [
    "**/__tests__/**/?(*.)(spec|test).js?(x)"
  ],
  collectCoverageFrom: [
    "**/*.js?(x)"
  ],
  setupFilesAfterEnv: [
    "./spec/setup.js"
  ],
  setupFiles: [
    "raf/polyfill"
  ]
};