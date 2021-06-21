module.exports = {
  globals: {
    __TEST__: true,
    CONFIG: {
      project_name: 'labors',
      magma_host: "https://magma.test"
    }
  },
  moduleNameMapper: {
    "\\.(css|scss)$": "identity-obj-proxy"
  },
  testMatch: [
    '**/__tests__/**/?(*.)(spec|test).(j|t)s?(x)'
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