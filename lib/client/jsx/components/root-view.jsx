import * as React from 'react';
import { connect } from 'react-redux';
import { selectUserPermissions } from '../selectors/user-selector';

export class RootView extends React.Component{
  render(){
    let { permissions } = this.props;

    if (Object.keys(permissions).length <= 0)
      return <div className='projects'>{'No available projects.'}</div>;

    return(
      <div className='root'>
        <div className='title'>
          {'Available Projects'}
        </div>
        <div className='projects'>
          {
            Object.values(permissions).sort(
              (p1, p2)=>(p1.role[0] + p1.project_name) > (p2.role[0] + p2.project_name)
            ).map( ({project_name,role, privileged}) =>
              <div key={project_name} className='project'>
                <a className='project_name'
                  key={project_name}
                  href={`/${project_name}`}>
                  {project_name}
                </a>
                - <span className='role'>{role}{ privileged ? ', privileged access' : '' }</span>
              </div>
            )
          }
        </div>
      </div>
    );
  }
}

export default connect(
  (state)=>({
    permissions: selectUserPermissions(state)
  })
)(RootView);
