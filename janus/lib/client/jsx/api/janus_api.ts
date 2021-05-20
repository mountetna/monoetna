import {json_get, json_post} from 'etna-js/utils/fetch';

export const fetchUsers = (): Promise<any> => json_get('/users');

export const updateUserFlags = ({
  email,
  flags
}: {
  email: string;
  flags: string[];
}): Promise<any> => json_post('/flag_user', {email, flags});

export const fetchProjects = (): Promise<any> => json_get('/projects');
