import React, { useState, useEffect, useCallback } from 'react';
import { json_post, checkStatus} from 'etna-js/utils/fetch';
import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import { isAdmin, isSuperuser } from 'etna-js/utils/janus';
import Icon from 'etna-js/components/icon';

const NEWUSER_INPUTS = {
  select: [ 'role' ],
  text: [ 'name', 'email' ]
};

const updateNewUser = function(e) {
  let select = e.target;
  let form = select.closest('.item');

  if (inputsChanged(form, NEWUSER_INPUTS)) showForm(form);
  else hideForm(form);
}

const cancelNewUser = function(e) {
  let select = e.target;
  let form = select.closest('.item');

  resetInputs(form, NEWUSER_INPUTS);
  hideForm(form);
}

const value = (item) => item.value || item.innerHTML;

const filter = (e) => {
  let filter = e.target.value;

  let items = document.querySelectorAll('.items .item');
  items.forEach(item => {
    let columns = [ 'name', 'email', 'role', 'affiliation' ];
    if (columns.some(value_class =>
      value(item.querySelector(`.${value_class}`)).match(new RegExp(filter, 'i'))
      )) {
      show(item, 'flex');
    } else hide(item);
  })
}

const postUpdatePermission = (project_name, email, revision) => json_post(`/update_permission/${project_name}`, {email, ...revision});
const postAddUser = (project_name, {email, name, role}) => json_post(`/add_user/${project_name}`, {email, name, role});

const Permission = ({roles, editable, onSave, permission={}, create}) => {
  let {user_email, user_name, role, affiliation, privileged} = permission;
  let [ revision, updateRevision ] = useState({});

  let update = useCallback( (key, value) => {
    let newRevision = {...revision};
    if (permission[key] == value) {
      delete newRevision[key];
    } else {
      newRevision[key] = value;
    }
    updateRevision(newRevision);
  }, [ permission, revision ]);

  return <div className='item permission'>
    <div className='cell'>
      {
        create
          ? <input type='text' placeholder='Name' name='name' onChange={ e => update('name', e.target.value) } value={ revision.name  || ''} />
          : <span className='name'>{ user_name }</span>
      }
    </div>
    <div className='cell'>
      {
        create
          ? <input type='text' placeholder='Email' name='email' onChange={ e => update('email', e.target.value) } value={ revision.email || '' } />
          : <span className='email'>{ user_email }</span>
      }
    </div>
    <div className='cell'>
      {
        editable
          ? <select className='role' name='role' value={revision.role || role || ''} onChange={ (e) => update('role', e.target.value) }>
             {
               roles.map(
                 r => <option key={r} value={ r }>{ r }</option>
               )
             }
           </select>
         : <span className='role'>{ role }</span>
      }
    </div>
    <div className='cell'>
      {
        editable
          ? <input className='affiliation' type='text' name='affiliation' value={ revision.affiliation || affiliation || ''}
              onChange={ (e) => update( 'affiliation', e.target.value ) } />
          : <span className='affiliation'>{ affiliation }</span>
      }
    </div>
    <div className='cell'>
      {
        editable && !create
          ? <input type='checkbox' name='privileged' checked={ 'privileged' in revision ? revision.privileged : privileged } onChange={ e => update('privileged', e.target.checked) } />
          : (privileged && 'Yes')
      }
    </div>
    <div className='cell submit'>
      {
        Object.keys(revision).length > 0 && <React.Fragment>
          <Icon className='approve' icon='save' onClick={ () => { onSave({revision,user_email}); updateRevision({}); } } />
          <Icon className='cancel' icon='ban' onClick={ () => updateRevision({}) }/>
        </React.Fragment>
      }
    </div>
  </div>
}

const userCount = (count, txt) => count == 0 ? null : `${count} ${txt}${ count === 1 ? '' : 's' }`;
const ProjectView = ({project_name}) => {
  let user = useReduxState( state => selectUser(state) );
  let [ project, setProject ] = useState({});
  let [ filter, setFilter ] = useState(null);

  let retrieveProject = useCallback(
    () => {
      fetch(`/project/${project_name}`).then(checkStatus)
        .then(({project}) => setProject(project))
    }, [project_name, setProject]
  );
  useEffect( () => retrieveProject(), [] );

  let { permissions=[], project_name_full } = project;

  let privileged = permissions.filter(p => p.privileged);

  let editable = isAdmin(user, project_name);
  let roles = isSuperuser(user,project_name)
    ?  [ 'administrator', 'editor', 'viewer', 'disabled' ]
    : isAdmin(user, project_name)
      ? [ 'editor', 'viewer', 'disabled' ]
      : [ ];

  return <div id='project-view'>
    <div className='title'> { project_name_full }</div>
    <div className='item summary'>
      {
        [ 'administrator', 'editor', 'viewer' ].map(
          role => {
            let ps = permissions.filter(p => p.role == role);
            return userCount(ps.length, role)
          }
        ).concat(
          userCount(privileged.length, 'privileged user')
        ).filter(_=>_).join(', ')
      }
    </div>
    <div className='item header'>
      <div className='cell'>Name</div>
      <div className='cell'>Email</div>
      <div className='cell'>Role</div>
      <div className='cell'>Affiliation</div>
      <div className='cell'>Privileged</div>
      <div className='cell submit'></div>
    </div>
    {
      permissions.length > 10 && <div className='item'>
        <input className='filter' type='text' placeholder='Filter rows' name='filter' onChange={ e => setFilter(e.target.value) }/>
      </div>
    }
    {
      permissions.sort((a,b) => (a.role+a.user_email).localeCompare(b.role+b.user_email) ).filter(
        p => !filter || [ 'name', 'user_email', 'role', 'affiliation' ].some( column => column in p && p[column].match(new RegExp(filter, 'i')))
      ).map(
        p => <Permission
          key={p.user_email}
          permission={p}
          editable={editable}
          roles={roles}
          onSave={({revision, user_email}) => postUpdatePermission(project_name, user_email, revision).then(() => retrieveProject())}
        />
      )
    }
    {
      editable && <React.Fragment>
        <div className='new'>New User</div>
        <Permission create={true} editable={true} roles={['Select role', 'viewer', 'editor']} 
          onSave={({revision}) => postAddUser(project_name, revision).then(() => retrieveProject()).catch()} />
      </React.Fragment>

    }
  </div>
}
export default ProjectView;
