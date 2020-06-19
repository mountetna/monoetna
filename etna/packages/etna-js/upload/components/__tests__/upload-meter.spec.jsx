import React from 'react';
import { mockStore } from '../../../spec/helpers';
import renderer from 'react-test-renderer';
import UploadMeter from '../upload-meter';

describe('UploadMeter', () => {
  let store;

  beforeEach(() => {
    store = mockStore({
      location: {
        path: '/labors/browse/monster/Nemean Lion'
      }
    });
  });

  it('matches snapshot when upload starting', () => {
    const upload = {
      status: 'paused',
      upload_speeds: [],
      file_size: 100,
      current_byte_position: 0
    };
    const tree = renderer
      .create(<UploadMeter upload={upload} store={store} />)
      .toJSON();

    expect(tree).toMatchSnapshot();
  });

  it('matches snapshot when upload is in progress', () => {
    const upload = {
      status: 'active',
      upload_speeds: [100, 500, 2000, 10],
      file_size: 100,
      current_byte_position: 50
    };
    const tree = renderer
      .create(<UploadMeter upload={upload} store={store} />)
      .toJSON();

    expect(tree).toMatchSnapshot();
  });

  it('matches snapshot when upload is completed', () => {
    const upload = {
      status: 'complete',
      upload_speeds: [100, 500, 2000, 10, 8, 25, 827, 42],
      file_size: 100,
      current_byte_position: 100
    };
    const tree = renderer
      .create(<UploadMeter upload={upload} store={store} />)
      .toJSON();

    expect(tree).toMatchSnapshot();
  });
});
