module.exports = {
  globals: {
    CONFIG: {
      project_name: 'labors',
      baseURL: 'http://localhost'
    },
  },
  transformIgnorePatterns: [
    // "node_modules/(?!(etna-js)/)"
  ],
  moduleNameMapper: {
    "^service-worker-loader\!": "<rootDir>/__mocks__/service-worker-loader.js",
    '^react$': '<rootDir>/node_modules/react',
    '^react-redux$': '<rootDir>/node_modules/react-redux',
    '^react-dom$': '<rootDir>/node_modules/react-dom',
    '^react-modal$': '<rootDir>/node_modules/react-modal',
    '^enzyme$': '<rootDir>/node_modules/enzyme',
    '^enzyme-adapter-react-16$':
      '<rootDir>/node_modules/enzyme-adapter-react-16'
  },
  testURL: 'http://localhost',
  testMatch: [ "**/__tests__/**/?(*.)(spec|test).(j|t)s?(x)" ],
  collectCoverageFrom: [ "**/*.js?(x)" ],
  setupFilesAfterEnv: [ "./spec/setup.js" ],
  setupFiles: [ "raf/polyfill" ]
};
