import React from 'react';
import {Provider} from 'react-redux';
import { act, create } from 'react-test-renderer';
import {mockStore} from '../test-helpers';
import {stubUrl} from 'etna-js/spec/helpers';
import JanusSettings from '../janus-settings';

describe('JanusSettings', () => {
  let store;

  beforeEach(() => {
    store = mockStore({
      user: {
        email: "janus@two-faces.org",
        name: "Janus Bifrons",
        permissions: { }
      }
    });
  });

  const createNodeMock = element => {
    if (element.type === "textarea") {
      return document.createElement("textarea");
    } else {
      return null;
    }
  };

  it('renders', async () => {
    const initialStubs = [
      stubUrl({
        verb: 'get',
	      path: '/api/user/info',
        status: 200,
        response: {
          user: {
            email: "janus@two-faces.org",                                  
            name: "Janus Bifrons",
            public_key: '74:68:69:73:69:73:6e:6f:74:72:61:6e:64:6f:6d:21'
          }
        }
      })
    ];

    let component = create(
      <Provider store={store}>
        <JanusSettings/>
      </Provider>,
      { createNodeMock }
    );

    await act( async () => {
        await new Promise((resolve) => setTimeout(resolve, 15));
    });

    expect( component.toJSON() ).toMatchSnapshot()
  });
});
