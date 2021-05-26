import React from 'react';
import {rest} from 'msw';
import {setupServer} from 'msw/node';
import {
  render,
  fireEvent,
  waitFor,
  screen,
  within
} from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';
import UserTable from '../flags-user-table';
import {flagsSpecWrapper} from '../../spec/flags-helpers';

const mockUsers = [
  {
    email: 'janus@two-faces.org',
    name: 'Janus Bifrons',
    projects: ['doors', 'portals'],
    flags: ['inside', 'outside']
  },
  {
    email: 'portunus@two-faces.org',
    name: 'Portunus',
    projects: ['doors'],
    flags: ['underground']
  }
];

const mockProjects = [
  {project_name: 'doors', project_name_full: 'Oak doors'},
  {project_name: 'portals', project_name_full: 'Magic portals'}
];

const mockState = {
  projects: mockProjects,
  users: mockUsers
};

const handlers = [
  rest.get('/users', (req, res, ctx) => {
    return res(ctx.json({users: mockUsers}));
  }),
  rest.get('/projects', (req, res, ctx) => {
    return res(ctx.json({projects: mockProjects}));
  })
];

const server = setupServer(...handlers);

beforeAll(() => server.listen());
afterEach(() => server.resetHandlers());
afterAll(() => server.close());

describe('UserTable', () => {
  it('renders correctly', async () => {
    let {asFragment} = render(<UserTable />, {
      wrapper: flagsSpecWrapper(mockState)
    });

    await waitFor(() => screen.getByText(/Janus/));

    expect(screen.queryByText('doors,portals')).toBeTruthy();
    expect(screen.queryByText('inside')).toBeTruthy();
    expect(screen.queryByText('outside')).toBeTruthy();
    expect(screen.getByText(/Janus/)).toHaveTextContent('Janus Bifrons');
    expect(screen.getByText(/Portunus/)).toHaveTextContent('Portunus');

    expect(asFragment()).toMatchSnapshot();
  });

  it('shows Add/Remove card when users are selected and Manage button clicked', async () => {
    let {asFragment} = render(<UserTable />, {
      wrapper: flagsSpecWrapper(mockState)
    });

    await waitFor(() => screen.getByText(/Janus/));

    fireEvent.click(screen.getByText(/Janus/));
    fireEvent.click(screen.getByText(/manage flags/i));

    await waitFor(() => {
      expect(screen.queryByText('Add Flag')).toBeTruthy();
      expect(screen.queryByText('Remove Flag')).toBeTruthy();

      expect(asFragment()).toMatchSnapshot();
    });
  });

  it('filters users based on text search', async () => {
    let {asFragment} = render(<UserTable />, {
      wrapper: flagsSpecWrapper(mockState)
    });

    await waitFor(() => screen.getByText(/Janus/));

    fireEvent.change(screen.getByRole('textbox', {name: ''}), {
      target: {value: 'unus'}
    });

    await waitFor(() => {
      expect(screen.queryByText(/Janus/)).toBeFalsy();
      expect(screen.queryByText(/Portunus/)).toBeTruthy();

      expect(asFragment()).toMatchSnapshot();
    });
  });

  it('filters users based on project search', async () => {
    let {asFragment} = render(<UserTable />, {
      wrapper: flagsSpecWrapper(mockState)
    });

    await waitFor(() => screen.getByText(/Janus/));

    fireEvent.mouseDown(screen.getByRole('textbox', {name: 'Projects'}));
    const listbox = within(screen.getByRole('listbox'));
    fireEvent.click(listbox.getByText(/portals/));

    await waitFor(() => {
      expect(screen.queryByText(/Janus/)).toBeTruthy();
      expect(screen.queryByText(/Portunus/)).toBeFalsy();

      expect(asFragment()).toMatchSnapshot();
    });
  });

  it('filters users based on flag search', async () => {
    let {asFragment} = render(<UserTable />, {
      wrapper: flagsSpecWrapper(mockState)
    });

    await waitFor(() => screen.getByText(/Janus/));

    fireEvent.mouseDown(screen.getByRole('textbox', {name: 'Flags'}));
    const listbox = within(screen.getByRole('listbox'));
    fireEvent.click(listbox.getByText(/outside/));

    await waitFor(() => {
      expect(screen.queryByText(/Janus/)).toBeTruthy();
      expect(screen.queryByText(/Portunus/)).toBeFalsy();

      expect(asFragment()).toMatchSnapshot();
    });
  });

  it('filters users based on all 3 search params', async () => {
    let {asFragment} = render(<UserTable />, {
      wrapper: flagsSpecWrapper(mockState)
    });

    await waitFor(() => screen.getByText(/Janus/));

    fireEvent.change(screen.getByRole('textbox', {name: ''}), {
      target: {value: 'nus'}
    });

    fireEvent.mouseDown(screen.getByRole('textbox', {name: 'Projects'}));
    let listbox = within(screen.getByRole('listbox'));
    fireEvent.click(listbox.getByText(/portals/));

    fireEvent.mouseDown(screen.getByRole('textbox', {name: 'Flags'}));
    listbox = within(screen.getByRole('listbox'));
    fireEvent.click(listbox.getByText(/underground/));

    await waitFor(() => {
      expect(screen.queryByText(/Janus/)).toBeFalsy();
      expect(screen.queryByText(/Portunus/)).toBeFalsy();

      expect(asFragment()).toMatchSnapshot();
    });
  });
});
