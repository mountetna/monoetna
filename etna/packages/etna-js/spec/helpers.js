import nock from 'nock';
import thunk from 'redux-thunk';
import configureMockStore from 'redux-mock-store';

export const mockStore = (state, middleware = []) => configureMockStore([thunk, ...middleware])(state);

export const stubUrl = ({
  verb = 'get',
  path,
  response,
  request,
  status = 200,
  headers = {},
  host = 'http://localhost'
}) => {
  let nocked;

  console.log('Stubbing for', {request, path, verb, response}, '\nactive')

  return new Promise((resolve, reject) => {
    const base = nock(host)[verb](path, request);
    nocked = response instanceof Function
      ? base.reply(function (uri, body, cb) {
        response(uri, body, function () {
          console.log('responding to', { uri, body }, '\nusing custom response cb')
          nock.removeInterceptor(base);
          resolve();
          return cb.apply(this, arguments);
        })
      })
      : base.reply(function () {
        console.log('responding to', { path, request }, '\nreplying with', { status, response })
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

export const cleanStubs = () => nock.cleanAll();

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
