import React from 'react';
import {render, screen, waitFor, fireEvent} from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';
import {ListEntryUpdatedColumn} from '../ListEntry';

describe('ListEntryUpdatedColumn', () => {
  it('renders correctly', async () => {
    const testObject = {
      updated_at: '2022-09-19T13:43:09+00:00',
      author: 'zeus@olympus.org|Zeus Almighty'
    };
    const testWidths = {};
    const {asFragment} = render(
      <ListEntryUpdatedColumn obj={testObject} widths={testWidths} />
    );

    await waitFor(() => screen.getByText(/by/));

    expect(screen.queryByText(/Zeus Almighty/)).toBeTruthy();
    expect(screen.queryByText(/zeus@olympus.org/)).toBeFalsy();
    expect(asFragment()).toMatchSnapshot();
  });
});
