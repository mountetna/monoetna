import React from 'react';
import nock from 'nock';
import thunk from 'redux-thunk';
import configureMockStore from 'redux-mock-store';
import fetch from 'isomorphic-fetch';

export const stubUrl = ({
  verb = 'get',
  path,
  response,
  request,
  status = 200,
  headers = {},
  host = 'http://localhost',
  times = 1
}: {
  verb: string;
  path: string;
  request: any;
  response: any;
  status: number;
  headers?: any;
  host: string;
  times?: number;
}) => {
  nock(host)
    [verb](path, request)
    .times(times)
    .reply(status, response, {
      'Access-Control-Allow-Origin': '*',
      'Content-type': 'application/json',
      ...headers
    });
};

export const cleanStubs = () => nock.cleanAll();

export const mockStore = configureMockStore([thunk]);

export const mockDate = () => {
  const currentDate = new Date();

  global.Date = jest.fn(() => currentDate) as any;
};

export const mockFetch = () => (global.fetch = fetch);
