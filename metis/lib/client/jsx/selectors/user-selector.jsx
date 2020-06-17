const ROLES = { a: 'administrator', e: 'editor', v: 'viewer' };

export const selectUserPermissions = ({user: {permissions}}) => permissions || {};

export const selectProjectPermissions = (state) => selectUserPermissions(state)[CONFIG.project_name] || {};

export const selectUserRole = (state) => selectProjectPermissions(state).role;

export const selectUserName = ({user: {first,last}}) => ({first,last});
