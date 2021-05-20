import {json_get, json_post} from 'etna-js/utils/fetch';

export const fetchUsers = (): Promise<any> => json_get('/users');

export const updateUserFlags = ({
  email,
  flags
}: {
  email: string;
  flags: string[];
}) => json_post('/flag_user', {email, flags});
