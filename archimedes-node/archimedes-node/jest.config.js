export default {
  testMatch: ['**/__tests__/**/?(*.)(spec|test).(j|t)s?(x)'],
  collectCoverageFrom: ['**/*.js?(x)'],
  setupFiles: ['raf/polyfill']
};
