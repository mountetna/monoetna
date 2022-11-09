module.exports = {
  globals: {
    fetch: require('isomorphic-fetch'),
    CONFIG: {
      project_name: 'labors',
      magma_host: 'https://magma.test',
      polyphemus_host: 'https://polyphemus.test',
      vulcan_host: 'https://vulcan.test',
      baseURL: 'http://localhost'
    }
  },
  testURL: 'http://localhost',
  transformIgnorePatterns: ['node_modules/(?!(etna-js|@codemirror|downzip)/)'],
  moduleNameMapper: {
    '^service-worker-loader!': '/app/__mocks__/service-worker-loader.js',
    '^.*[.](svg|css|CSS)$': 'identity-obj-proxy',
    '^react$': '/etna/node_modules/react',
    '^react-redux$': '/etna/node_modules/react-redux',
    '^react-codemirror2$': '/etna/node_modules/react-codemirror2',
    '^react-dom$': '/etna/node_modules/react-dom',
    '^react-modal$': '/etna/node_modules/react-modal',
    '^enzyme$': '/etna/node_modules/enzyme',
    '^enzyme-adapter-react-16$': '/etna/node_modules/enzyme-adapter-react-16',
    '^@material-ui/core/styles': '/etna/node_modules/@material-ui/core/styles',
    '^color$': '/etna/node_modules/color',
    '^color-string$': '/etna/node_modules/color-string',
  },
  testMatch: [
    '**/test/**/?(*.)(spec|test).(j|t)s?(x)',
    '**/__tests__/**/?(*.)(spec|test).(j|t)s?(x)'
  ],
  collectCoverageFrom: ['**/*.js?(x)'],
  setupFilesAfterEnv: ['./setup-jest.js'],
  setupFiles: ['raf/polyfill']
};
