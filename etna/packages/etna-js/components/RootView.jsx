require('./RootView.css');
import React, {useEffect} from 'react';
import {useReduxState} from '../hooks/useReduxState';

// Module imports.
import {selectUserPermissions} from '../selectors/user-selector';
import {selectProjects} from '../selectors/janus-selector';
import {projectNameFull} from '../utils/janus';
import {fetchProjectsAction} from '../actions/janus-actions';
import {useActionInvoker} from '../hooks/useActionInvoker';
import {useFeatureFlag} from '../hooks/useFeatureFlag';

const Project = ({project_name, project_name_full, role, privileged}) => (
  <div className='project'>
    <div className='project_name'>
      <a href={`/${project_name}`}>{project_name}</a>
    </div>
    <div className='role'>
      {role || 'viewer'}
      {privileged ? ', privileged access' : ''}
    </div>
    <div className='full_name'>{project_name_full}</div>
  </div>
);

const ProjectList = ({title, projects}) => {
  return (
    <>
      <div className='title'>{title}</div>
      <div className='projects'>
        <div className='project header'>
          <div className='project_name'> project_name </div>
          <div className='role'> role </div>
          <div className='full_name'> title </div>
        </div>
        <div>
          <div>
            {projects.length == 0 ? (
              <span className='empty-list'>{'No available projects.'}</span>
            ) : (
              projects.map((project, i) => <Project key={i} {...project} />)
            )}
          </div>
        </div>
      </div>
    </>
  );
};

const RootView = () => {
  const invoke = useActionInvoker();
  const canCommunity = !useFeatureFlag('external');

  const {my_projects, resourceProjects, communityProjects} = useReduxState(
    (state) => {
      let permissions = selectUserPermissions(state);
      let projects = selectProjects(state);

      let my_projects = Object.values(permissions)
        .map(({project_name, ...perms}) => ({
          project_name,
          ...perms,
          project_name_full: projectNameFull(projects, project_name)
        }))
        .sort(({project_name: n1}, {project_name: n2}) =>
          n1 < n2 ? -1 : n1 > n2 ? 1 : 0
        );

      const resourceProjects = projects.filter(
        (proj) => proj.resource && !proj.requires_agreement
      );
      const communityProjects = projects.filter(
        (proj) => proj.resource && proj.requires_agreement
      );

      return {my_projects, resourceProjects, communityProjects};
    }
  );

  useEffect(() => {
    invoke(fetchProjectsAction());
  }, []);

  return (
    <div className='root-view'>
      <ProjectList title='Your Projects' projects={my_projects} />
      <ProjectList title='Resource Projects' projects={resourceProjects} />
      {canCommunity && <ProjectList title='Community Projects' projects={communityProjects} />}
    </div>
  );
};

export default RootView;
