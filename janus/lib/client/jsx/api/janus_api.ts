import {json_get} from 'etna-js/utils/fetch';

export const fetchUsers = (): Promise<any> => json_get('/users');
