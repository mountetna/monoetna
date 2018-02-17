const roles = { a: 'administrator', e: 'editor', v: 'viewer' };

const parsePermissions = (perms) => {
  // permissions are encoded as 'a:project1,project2;v:project3'
  let permissions = perms.split(/;/).map(perm => {
    let [ role, projects ] = perm.split(/:/);
    role = roles[role.toLowerCase()];
    return projects.split(/,/).map(
      project_name => ({ role, project_name })
    )
  }).reduce((perms,perm) => perms.concat(perm), []);

  return permissions;
}

const parseToken = (token) => {
  let [ header, params, signature ] = token.split(/\./);

  let { email, first, last, perm } = JSON.parse(atob(params));

  return {
    userEmail: email,
    firstName: first,
    lastName: last,
    permissions: parsePermissions(perm)
  };
}

const user = (state, action) => {
  if (!state) state = {
    userEmail: '',
    firstName: '',
    lastName: '',
    permissions: [],

    masterPerms: false,

    loginStatus: false,
    loginError: false,
    loginErrorMsg: 'Invalid sign in.'
  };

  switch(action.type){
    case 'ADD_USER':
      return {
        ...state,
        ...parseToken(action.token)
      }

    default:
      return state;
  }
};

export default user;
