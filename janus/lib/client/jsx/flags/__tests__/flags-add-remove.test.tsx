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
import AddRemove from '../flags-add-remove';
import {FlagsProvider} from '../flags-context';

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
  rest.post('/flag_user', (req, res, ctx) => {
    return res(ctx.json({}));
  }),
  rest.get('/users', (req, res, ctx) => {
    return res(ctx.json({users: mockUsers}));
  })
];

const server = setupServer(...handlers);

beforeAll(() => server.listen());
afterEach(() => server.resetHandlers());
afterAll(() => server.close());

describe('AddRemove', () => {
  it('renders correctly', async () => {
    const wrapper = ({children}: any) => (
      <FlagsProvider state={mockState}>{children}</FlagsProvider>
    );

    let {asFragment} = render(
      <AddRemove selectedUsers={mockUsers.slice(0, 1)} />,
      {wrapper}
    );

    await waitFor(() => screen.getByText(/1 user/));

    expect(screen.queryByText('Add Flag')).toBeTruthy();
    expect(screen.queryByText('Remove Flag')).toBeTruthy();

    expect(asFragment()).toMatchSnapshot();
  });

  it('appends flags to existing flags', async () => {
    const wrapper = ({children}: any) => (
      <FlagsProvider state={mockState}>{children}</FlagsProvider>
    );

    let {asFragment} = render(
      <AddRemove selectedUsers={mockUsers.slice(0, 1)} />,
      {wrapper}
    );

    await waitFor(() => screen.getByText(/1 user/));

    fireEvent.change(screen.getByRole('textbox', {name: ''}), {
      target: {value: 'newflag'}
    });
    fireEvent.click(screen.getByText('Add Flag'));

    await waitFor(() => {
      expect(screen.queryByText(/Saved/)).toBeTruthy();

      expect(asFragment()).toMatchSnapshot();
    });
  });

  it('removes flags from users', async () => {
    const wrapper = ({children}: any) => (
      <FlagsProvider state={mockState}>{children}</FlagsProvider>
    );

    let {asFragment} = render(
      <AddRemove selectedUsers={mockUsers.slice(0, 1)} />,
      {wrapper}
    );

    await waitFor(() => screen.getByText(/1 user/));

    fireEvent.change(screen.getByRole('textbox', {name: ''}), {
      target: {value: 'newflag'}
    });
    fireEvent.click(screen.getByText('Remove Flag'));

    await waitFor(() => {
      expect(screen.queryByText(/Saved/)).toBeTruthy();

      expect(asFragment()).toMatchSnapshot();
    });
  });

  it('reports errors to the user', async () => {
    server.use(
      rest.post('/flag_user', (req, res, ctx) => {
        // Respond with "422 Error" status for this test.
        return res(
          ctx.status(422),
          ctx.json({error: 'Must be an array of strings'})
        );
      })
    );

    const wrapper = ({children}: any) => (
      <FlagsProvider state={mockState}>{children}</FlagsProvider>
    );

    let {asFragment} = render(
      <AddRemove selectedUsers={mockUsers.slice(0, 1)} />,
      {wrapper}
    );

    await waitFor(() => screen.getByText(/1 user/));

    fireEvent.change(screen.getByRole('textbox', {name: ''}), {
      target: {value: 'bad-flag'}
    });
    fireEvent.click(screen.getByText('Add Flag'));

    await waitFor(() => {
      expect(screen.queryByText(/Must be an array of strings/)).toBeTruthy();

      expect(asFragment()).toMatchSnapshot();
    });
  });
});
