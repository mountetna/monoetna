import React, {useState, useEffect} from 'react';
import Identity from './identity';
import UserProjects from './user-projects';
import {selectUser} from 'etna-js/selectors/user-selector';
import {checkStatus} from 'etna-js/utils/fetch';
import {useReduxState} from 'etna-js/hooks/useReduxState';

const JanusRoot = () => {
  let user = useReduxState( state => selectUser(state) );

  let [ projects, setProjects ]  = useState([]);

  useEffect(
    () => fetch('/projects?full=1').then(checkStatus).then(
      ({projects}) => setProjects(projects)
    ), []
  )
  return <div>
    <Identity user={user}/>
    <UserProjects projects={projects}/>
  </div>;
}

export default JanusRoot;
