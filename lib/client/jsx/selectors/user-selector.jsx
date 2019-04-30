const ROLES = { a: 'administrator', e: 'editor', v: 'viewer' };

const selectPermissions = ({user: {perm}}) => {
  // permissions are encoded as 'a:project1,project2;v:project3'
  let token_roles = perm.split(/;/).map(role_projects => {
    let [ token_role, projects ] = role_projects.split(/:/);
    let project_names = projects.split(/,/);
    return { token_role, project_names };
  })

  let { token_role: project_role } = token_roles.find(({project_names}) => project_names.find(project_name => project_name == CONFIG.project_name));

  return project_role ? {
    role: ROLES[project_role.toLowerCase()],
    restricted: (project_role < 'a')
  } : {};
}

export const selectUserRole = (state) => selectPermissions(state).role;

export const selectUserName = ({user: {first,last}}) => ({first,last});
