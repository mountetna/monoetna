import React from 'react';
import { mount, shallow } from 'enzyme';
import { mockStore } from '../../../spec/helpers';
import renderer from 'react-test-renderer';
import UploadControl from '../upload-control';
import {
  CANCEL_UPLOAD,
  CONTINUE_UPLOAD,
  PAUSE_UPLOAD
} from '../../actions/upload_actions';

describe('UploadControl', () => {
  let store;

  beforeEach(() => {
    store = mockStore({
      location: {
        path: '/labors/browse/monster/Nemean Lion'
      }
    });
  });

  it('renders the play and cancel buttons when status is paused', () => {
    const component = mount(
      <UploadControl upload={{ status: 'paused' }} store={store} />
    );

    const buttonGroup = component.find('.upload-control-group');
    expect(buttonGroup.exists()).toBeTruthy();

    const playButton = component.find('.fa-play');
    expect(playButton.exists()).toBeTruthy();

    const cancelButton = component.find('.fa-times');
    expect(cancelButton.exists()).toBeTruthy();

    const tree = renderer
      .create(<UploadControl upload={{ status: 'paused' }} store={store} />)
      .toJSON();

    expect(tree).toMatchSnapshot();
  });

  it('renders the pause and cancel buttons when status is active', () => {
    const component = mount(
      <UploadControl upload={{ status: 'active' }} store={store} />
    );

    const buttonGroup = component.find('.upload-control-group');
    expect(buttonGroup.exists()).toBeTruthy();

    const pauseButton = component.find('.fa-pause');
    expect(pauseButton.exists()).toBeTruthy();

    const cancelButton = component.find('.fa-times');
    expect(cancelButton.exists()).toBeTruthy();

    const tree = renderer
      .create(<UploadControl upload={{ status: 'active' }} store={store} />)
      .toJSON();

    expect(tree).toMatchSnapshot();
  });

  it('renders the retry and cancel buttons when status is failed', () => {
    const component = mount(
      <UploadControl upload={{ status: 'failed' }} store={store} />
    );

    const buttonGroup = component.find('.upload-control-group');
    expect(buttonGroup.exists()).toBeTruthy();

    const retryButton = component.find('.fa-retweet');
    expect(retryButton.exists()).toBeTruthy();

    const cancelButton = component.find('.fa-times');
    expect(cancelButton.exists()).toBeTruthy();

    const tree = renderer
      .create(<UploadControl upload={{ status: 'failed' }} store={store} />)
      .toJSON();

    expect(tree).toMatchSnapshot();
  });

  it('dispatches CANCEL_UPLOAD action when cancel button clicked', () => {
    expect(store.getActions()).toEqual([]);
    const component = mount(
      <UploadControl upload={{ status: 'active' }} store={store} />
    );

    const cancelButton = component.find('.fa-times').first();
    cancelButton.simulate('click');

    expect(store.getActions()).toEqual([
      {
        type: CANCEL_UPLOAD,
        upload: {
          status: 'active'
        }
      }
    ]);
  });

  it('dispatches PAUSE_UPLOAD action when pause button clicked', () => {
    expect(store.getActions()).toEqual([]);
    const component = mount(
      <UploadControl upload={{ status: 'active' }} store={store} />
    );

    const pauseButton = component.find('.fa-pause').first();
    pauseButton.simulate('click');

    expect(store.getActions()).toEqual([
      {
        type: PAUSE_UPLOAD,
        upload: {
          status: 'active'
        }
      }
    ]);
  });

  it('dispatches CONTINUE_UPLOAD action when play button clicked', () => {
    expect(store.getActions()).toEqual([]);
    const component = mount(
      <UploadControl upload={{ status: 'paused' }} store={store} />
    );

    const playButton = component.find('.fa-play').first();
    playButton.simulate('click');

    expect(store.getActions()).toEqual([
      {
        type: CONTINUE_UPLOAD,
        upload: {
          status: 'paused'
        }
      }
    ]);
  });

  // This functionality doesn't seem to currently exist?
  xit('dispatches SELECT_UPLOAD action when retry button clicked', () => {
    expect(store.getActions()).toEqual([]);
    const component = mount(
      <UploadControl upload={{ status: 'failed' }} store={store} />
    );

    const retweetButton = component.find('.fa-retweet').first();
    retweetButton.simulate('click');

    expect(store.getActions()).toEqual([
      {
        type: SELECT_UPLOAD,
        upload: {
          status: 'failed'
        }
      }
    ]);
  });
});
