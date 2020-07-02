import nock from 'nock';
import thunk from 'redux-thunk';
import configureMockStore from 'redux-mock-store';

export const mockStore = configureMockStore([thunk]);

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

  return new Promise((resolve, reject) => {
    nocked = nock(host)
      [verb](path, request)
      .reply(status, () => {
        console.log('stubUrl: Received request matching', { verb, path, request });
        resolve(nocked);
        return response;
      }, {
        'Access-Control-Allow-Origin': '*',
        'Content-type': 'application/json',
        ...headers
      });
  });
};

export const cleanStubs = () => nock.cleanAll();

global.expect.extend({
  toResolve: async (promise) => {
    const resolved = await promise.then(() => true, () => false);
    if (resolved) {
      return { pass: true, message: 'Expected promise to resolve, and it did' };
    }

    return { pass: false, message: 'Expected promise to resolve, but it rejected' };
  }
})
