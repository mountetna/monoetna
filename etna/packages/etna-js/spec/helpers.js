import nock from 'nock';
import thunk from 'redux-thunk';
import configureMockStore from 'redux-mock-store';
import {appSubscription} from "../utils/subscription";

export const mockStore = (state, middleware = []) => configureMockStore([thunk, ...middleware])(state);

function nockDebug(...args) {
  if (process.env['NOCK_DEBUG'] === '1') {
    console.log(...args);
  }
}

export const stubUrl = ({
  verb = 'get',
  path,
  response,
  request,
  status = 200,
  headers = {},
  host = 'http://localhost',
  url,
}) => {
  let nocked;

  if (url) {
    if (!(url instanceof URL)) {
      url = new URL(url);
    }

    host = url.protocol + "//" + url.hostname;
    path = url.pathname;
    if (url.search) {
      path += url.search;
    }
  }

  nockDebug('Stubbing for', {request, path, verb, response, host}, '\nactive')

  return new Promise((resolve, reject) => {
    const base = nock(host)[verb](path, request);
    nocked = response instanceof Function
      ? base.reply(function (uri, body, cb) {
        response(uri, body, function () {
          nockDebug('responding to', { uri, body }, '\nusing custom response cb')
          nock.removeInterceptor(base);
          resolve();
          return cb.apply(this, arguments);
        })
      })
      : base.reply(function () {
        nockDebug('responding to', { path, request }, '\nreplying with', { status, response })
        nock.removeInterceptor(base);
        resolve();
        return [status, response, {
          'Access-Control-Allow-Origin': '*',
          'Content-type': 'application/json',
          ...headers
        }]
      })
  });
};

export async function joinedDeferredPromises(...promiseChains) {
  const result = [];
  for (let promiseChain of promiseChains) {
    const deferredChain = await Promise.all(promiseChain);
    result.push(deferredChain.reduce((p, n) => p.then(n), Promise.resolve()));
  }
  return result;
}

// This remains here for backwards compatibility, but moving forward it is invoked automatically
// by the bellow for each blocks as well, so doing this in individual tests should not be required.
export const cleanStubs = () => nock.cleanAll();

beforeEach(() => {
  appSubscription.addCleanup(cleanStubs);
})

afterEach(() => {
  appSubscription.end();
})

export function delay(ms) {
  return new Promise((resolve) => {
    setTimeout(resolve, ms);
  });
}

global.expect.extend({
  toResolve: async (promise) => {
    const resolved = await promise.then(() => true, () => false);
    if (resolved) {
      return {pass: true, message: 'Expected promise to resolve, and it did'};
    }

    return {pass: false, message: 'Expected promise to resolve, but it rejected'};
  }
})

export const mockDate = () => {
  const currentDate = new Date();
  global.Date = jest.fn(() => currentDate);
}

export const mockFetch = () => global.fetch = fetch;
