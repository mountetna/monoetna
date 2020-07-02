import 'regenerator-runtime/runtime';
// suppresses a specific React warning (comment raf out to bring it back)
import nock from 'nock';
import { stringifyRequest } from 'nock/lib/common';
const raf = (global.requestAnimationFrame = (cb) => setTimeout(cb, 0));

const Enzyme = require('enzyme');
const EnzymeAdapter = require('enzyme-adapter-react-16');
// Setup enzyme's react adapter
Enzyme.configure({ adapter: new EnzymeAdapter() });

global.fetch = require('node-fetch');

nock.disableNetConnect();
nock.emitter.on('no match', (...args) => {
  // nock has inconsistent parameters for this event based on how the match failed, which results in us having complex
  // unpacking logic.
  let [ _, options, bodyString ] = args;
  if (options == null) {
    ({ options } = args[0]);
    if (options == null) {
      options = args[0];
    }
  }

  const reqString = stringifyRequest(options, bodyString);
  console.error('got request with no match', reqString);
})

// nodejs equivalent
global.FormData = URLSearchParams;