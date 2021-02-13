import React, {useState, useEffect} from 'react';
import Identity from './identity';
import UserProjects from './user-projects';
import {selectUser} from 'etna-js/selectors/user-selector';
import {json_get} from 'etna-js/utils/fetch';
import {useReduxState} from 'etna-js/hooks/useReduxState';

const JanusMain = () => {
  let user = useReduxState( state => selectUser(state) );

  let [ projects, setProjects ]  = useState([]);

  useEffect(
    () => {
      json_get('/projects')
        .then(({projects}) => setProjects(projects))
    }, []
  )
  return <div id='janus-main'>
    <Identity user={user}/>
    <UserProjects user={user} projects={projects}/>
  </div>;
}

export default JanusMain;
