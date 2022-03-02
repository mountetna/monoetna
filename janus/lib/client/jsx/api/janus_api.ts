import {json_get, json_post} from 'etna-js/utils/fetch';

export const fetchUsers = (): Promise<any> => json_get('/api/users');

export const updateUserFlags = ({
  email,
  flags
}: {
  email: string;
  flags: string[];
}): Promise<any> => json_post('/api/admin/flag_user', {email, flags});

export const fetchProjects = (): Promise<any> => json_get('/api/user/projects');

export const updateProject = (
  project_name: string,
  revisions: {[key: string]: any}
) => json_post(`/api/admin/${project_name}/update`, revisions);
