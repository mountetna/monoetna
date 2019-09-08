const ROLES = {a: 'administrator', e: 'editor', v: 'viewer'};

const parsePermissions = (perms) => {
  // Permissions are encoded as 'a:project1,project2;v:project3'
  let permissions = perms.split(/;/).map(perm => {
    let [ encoded_role, projects ] = perm.split(/:/);
    let role = ROLES[encoded_role.toLowerCase()];
    let privileged = encoded_role.toUpperCase() == encoded_role;
    return projects.split(/,/).map(
      project_name=>({role, privileged, project_name})
    )
  }).reduce((perms,projects) => {
    projects.forEach(perm=>perms[perm.project_name] = perm);
    return perms
  }, {});

  return permissions;
}

const parseToken = (token)=>{
  let [header, params, signature] = token.split(/\./);
  let {email, first, last, perm} = JSON.parse(atob(params));

  return {
    email,
    first,
    last,
    permissions: parsePermissions(perm)
  };
}

const userReducer = function(user, action) {
  if (!user) user = { }
  switch(action.type) {
    case 'ADD_USER':
      return parseToken(action.token);
    default:
      return user;
  }
}

export default userReducer;
