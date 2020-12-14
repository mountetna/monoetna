module.exports = {
  globals: {
    __TEST__: true,
    CONFIG: {
      project_name: 'labors'
    }
  },
  moduleNameMapper: {
    "\\.(css|scss)$": "identity-obj-proxy"
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