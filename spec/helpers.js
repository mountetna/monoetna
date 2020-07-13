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
  host = 'http://www.fake.com'
}) => {
  nock(host)
    [verb](path, request)
    .reply(status, response, {
      'Access-Control-Allow-Origin': '*',
      'Content-type': 'application/json',
      ...headers
    });
};

export const cleanStubs = () => nock.cleanAll();
