import React, {useState, useEffect, useCallback} from 'react';

import {json_get} from 'etna-js/utils/fetch';

import {UserInterface} from './models/user_models';

const UsersView = () => {
  let [users, setUsers] = useState([] as UserInterface[]);

  useEffect(() => {
    json_get('/users').then(({users}) => {
      setUsers(users);
    });
  }, []);
  return (
    <div id='janus-users-main'>
      <Identity user={user} />
      <UserProjects user={user} projects={projects} />
    </div>
  );
};

export default UsersView;
