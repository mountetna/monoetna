const ROLES = {a: 'administrator', e: 'editor', v: 'viewer'};

const parsePermissions = (perms) => {
  // Permissions are encoded as 'a:project1,project2;v:project3'
  return perms.split(/;/).map(perm => {
    let [encoded_role, projects] = perm.split(/:/);
    let role = ROLES[encoded_role.toLowerCase()];
    let privileged = encoded_role.toUpperCase() == encoded_role;
    return projects.split(/,/).map(
      project_name => ({ role, privileged, project_name })
    )
  }).reduce((perms, projects) => {
    projects.forEach(perm => perms[perm.project_name] = perm);
    return perms
  }, {});
}

export function parseToken(token) {
  let [header, params, signature] = token.split(/\./);
  let {email, first, last, perm} = JSON.parse(atob(params));

  return {
    email,
    first,
    last,
    permissions: parsePermissions(perm)
  };
}

export const projectNameFull = (projects, project_name) => {
  if (!projects) return project_name;

  let matching_project = projects.find((project) => {
    return project.project_name === project_name;
  });

  return matching_project ? matching_project.project_name_full : project_name;
};
