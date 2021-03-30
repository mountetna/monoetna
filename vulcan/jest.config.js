module.exports = {
  globals: {
    fetch: require('node-fetch'),
    CONFIG: {
      project_name: 'labors',
      magma_host: 'https://magma.test',
      vulcan_host: 'https://vulcan.test'
    }
  },
  testURL: 'http://localhost',
  transformIgnorePatterns: ['node_modules/(?!(etna-js)/)'],
  moduleNameMapper: {
    '^service-worker-loader!': '<rootDir>/__mocks__/service-worker-loader.js',
    '^.*[.](css|CSS)$': 'identity-obj-proxy',
    '^react$': '<rootDir>/node_modules/react',
    '^react-redux$': '<rootDir>/node_modules/react-redux',
    '^react-dom$': '<rootDir>/node_modules/react-dom',
    '^enzyme$': '<rootDir>/node_modules/enzyme',
    '^enzyme-adapter-react-16$':
      '<rootDir>/node_modules/enzyme-adapter-react-16'
  },
  testMatch: ['**/__tests__/**/?(*.)(spec|test).(j|t)s?(x)'],
  collectCoverageFrom: ['**/*.js?(x)'],
  setupFilesAfterEnv: ['./lib/client/jsx/spec/setup.js'],
  setupFiles: ['raf/polyfill']
};
