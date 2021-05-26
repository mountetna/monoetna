import React, {useState, useEffect} from 'react';
import Identity from './identity';
import UserProjects from './user-projects';
import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {fetchProjects} from './api/janus_api'

const JanusMain = () => {
  let user = useReduxState( state => selectUser(state) );

  let [ projects, setProjects ]  = useState([]);

  useEffect(
    () => {
      fetchProjects()
      .then(({projects}) => {
        setProjects(projects);
        console.log({projects});
      })
    }, []
  )
  return <div id='janus-main'>
    <Identity user={user}/>
    <UserProjects user={user} projects={projects}/>
  </div>;
}

export default JanusMain;
