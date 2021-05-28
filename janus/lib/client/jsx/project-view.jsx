import React, { useState, useEffect, useCallback } from 'react';
import { json_post, json_get } from 'etna-js/utils/fetch';
import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import { isAdmin, isSuperuser } from 'etna-js/utils/janus';
import Icon from 'etna-js/components/icon';

import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Checkbox from '@material-ui/core/Checkbox';
import TextField from '@material-ui/core/TextField';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';

const useStyles = makeStyles((theme) => ({
  table: {
    minWidth: 650
  },
  role: {
    width: '300px'
  }
}));

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

  const classes = useStyles();

  return <TableRow>
    <TableCell>
      {
        create
          ? <TextField
            placeholder='New User Name'
            value={ revision.name  || ''}
            onChange={ e => update('name', e.target.value) }/>
          : user_name
      }
    </TableCell>
    <TableCell>
      {
        create
          ? <TextField
            placeholder='Email'
            value={ revision.email  || ''}
            onChange={ e => update('email', e.target.value) }/>
          : user_email
      }
    </TableCell>
    <TableCell>
      {
        editable
          ? <Select
            className={classes.role}
            value={revision.role || role || ''}
            onChange={ (e) => update('role', e.target.value) }>
             {
               roles.map(
                 r => <MenuItem key={r} value={ r }>{ r }</MenuItem>
               )
             }
           </Select>
         : <span className='role'>{ role }</span>
      }
    </TableCell>
    <TableCell>
      {
        editable
          ? <TextField
            value={ revision.affiliation || affiliation || ''}
            onChange={ (e) => update( 'affiliation', e.target.value ) }
          />
          : <span className='affiliation'>{ affiliation }</span>
      }
    </TableCell>
    <TableCell padding='checkbox'>
      {
        editable && !create
          ? <Checkbox
            checked={ 'privileged' in revision ? revision.privileged : privileged }
            onChange={ e => update('privileged', e.target.checked) }
            inputProps={{'aria-label': 'privileged user'}}
          />
          : (privileged && 'Yes')
      }
    </TableCell>
    <TableCell>
      {
        Object.keys(revision).length > 0 && <React.Fragment>
          <Icon className='approve' icon='save' onClick={ () => { onSave({revision,user_email}); updateRevision({}); } } />
          <Icon className='cancel' icon='ban' onClick={ () => updateRevision({}) }/>
        </React.Fragment>
      }
    </TableCell>
  </TableRow>
}

const userCount = (count, txt) => count == 0 ? null : `${count} ${txt}${ count === 1 ? '' : 's' }`;

const displayPermissions = (permissions, filter) => (
  permissions.sort(
    (a,b) => (a.role+a.user_email).localeCompare(b.role+b.user_email)
  ).filter(
    p => !filter || [ 'user_name', 'user_email', 'role', 'affiliation' ].some(
      column => column in p && null != p[column] && p[column].match(new RegExp(filter, 'i'))
    )
  )
);

const ProjectView = ({project_name}) => {
  let user = useReduxState( state => selectUser(state) );
  let [ project, setProject ] = useState({});
  let [ filter, setFilter ] = useState(null);

  let retrieveProject = useCallback(
    () => {
      json_get(`/project/${project_name}`)
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
  const classes = useStyles();

  return <div id='project-view'> 
    <Grid container direction='column'>
      <div className='title'> { project_name_full }</div>
      <div className='summary'>
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

      <TableContainer className={classes.table} component={Paper}>
        <Table aria-label='project users'>
          <TableHead>
            <TableRow>
              <TableCell>Name</TableCell>
              <TableCell>Email</TableCell>
              <TableCell>Role</TableCell>
              <TableCell>Affiliation</TableCell>
              <TableCell padding='checkbox'>Privileged</TableCell>
              <TableCell width="80"></TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
          {
            permissions.length > 2 && <TableRow>
              <TableCell colspan='6'>
                <TextField placeholder='Filter rows' onChange={ e => setFilter(e.target.value) }/>
              </TableCell>
            </TableRow>
          }
          {
            displayPermissions(permissions, filter).map(
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
            editable && <Permission
              create={true}
              editable={true}
              roles={['viewer', 'editor']} 
              onSave={({revision}) => postAddUser(project_name, revision).then(() => retrieveProject()).catch()} />
          }
          </TableBody>
        </Table>
      </TableContainer>
      {
        permissions.length > 10 && <div className='item'>
          <input className='filter' type='text' placeholder='Filter rows' name='filter' onChange={ e => setFilter(e.target.value) }/>
        </div>
      }
    </Grid>
  </div>
}
export default ProjectView;
