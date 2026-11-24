import {REQUEST_FAILED, requestError} from '../../../lib/client/jsx/components/cross_project_search/request_error';

describe('requestError', () => {
  it('reads the error a failed etna-js fetch rejects with', async () => {
    expect(await requestError(Promise.resolve({error: 'Label is required'}))).toEqual('Label is required');
    expect(await requestError(Promise.resolve({errors: ['a', {message: 'b'}]}))).toEqual('a, b');
    expect(await requestError(Promise.resolve('Bad gateway'))).toEqual('Bad gateway');
    expect(await requestError(Promise.resolve({unexpected: true}))).toEqual(REQUEST_FAILED);
  });

  it('handles a body that cannot be read and errors that are not responses', async () => {
    expect(await requestError(Promise.reject(new Error('unreadable')))).toEqual(REQUEST_FAILED);
    expect(await requestError(new TypeError('Failed to fetch'))).toEqual('Failed to fetch');
    expect(await requestError(undefined)).toEqual(REQUEST_FAILED);
  });
});
