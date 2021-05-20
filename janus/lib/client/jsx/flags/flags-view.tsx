import React, {useState, useEffect, useCallback} from 'react';

import {fetchUsers} from '../api/janus_api';

import {UserFlagsInterface} from '../models/user_models';

const User = (user: UserFlagsInterface) => {
  return <div>{user.email}</div>;
};

const FlagsView = () => {
  let [users, setUsers] = useState([] as UserFlagsInterface[]);

  useEffect(() => {
    fetchUsers().then(({users}) => {
      setUsers(users);
    });
  }, []);
  return <div id='janus-users-main'></div>;
};

export default FlagsView;
