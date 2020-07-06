// suppresses a specific React warning (comment raf out to bring it back)
import 'regenerator-runtime/runtime';
import nock from 'nock';
const raf = (global.requestAnimationFrame = (cb) => setTimeout(cb, 0));

const Enzyme = require('enzyme');
const EnzymeAdapter = require('enzyme-adapter-react-16');
// Setup enzyme's react adapter
Enzyme.configure({ adapter: new EnzymeAdapter() });

global.fetch = require('node-fetch');
nock.disableNetConnect();
nock.emitter.on('no match', req => {
  console.error('got request with no match', { headers: req.options.headers, href: req.options.href, method: req.options.method });
})
