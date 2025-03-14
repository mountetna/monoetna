import React, {useState, useEffect, useCallback, useMemo} from 'react';
import {json_post, json_get} from 'etna-js/utils/fetch';
import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {isAdmin, isSuperuser} from 'etna-js/utils/janus';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
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

import {basicSetup, EditorView} from 'codemirror';
import {markdown} from '@codemirror/lang-markdown';
import {defaultHighlightStyle, syntaxHighlighting} from '@codemirror/language';
import {EditorState} from '@codemirror/state';
import {gutter, lineNumbers} from '@codemirror/view';
import CodeMirror from 'rodemirror';
import {lintGutter} from '@codemirror/lint';

import ErrorBoundary from 'etna-js/components/error_boundary';

import SaveCancel from './save-cancel';
import {projectTypeFor} from './utils/project';

const useStyles = makeStyles((theme) => ({
  table: {
    minWidth: 650
  },
  title: {
    padding: '10px 15px'
  },
  subheading: {
    padding: '6px 15px',
    background: '#eee',
    borderTop: '1px solid #4431c6'
  },
  summary: {
    fontSize: '0.8em',
    color: '#1d6a86',
    fontStyle: 'italic',
    padding: '10px 15px',
    background: '#eee'
  },
  role: {
    width: '300px'
  },
  editor: {
    border: '1px solid #ccc',
    height: '100px',
    resize: 'vertical',
    overflow: 'hidden'
  }
}));

const postUpdatePermission = (project_name, email, revision) =>
  json_post(`/api/admin/${project_name}/update_permission`, {
    email,
    ...revision
  });
const postUpdateProject = (project_name, revision) =>
  json_post(`/api/admin/${project_name}/update`, revision);
const postAddUser = (project_name, revision) =>
  json_post(`/api/admin/${project_name}/add_user`, revision);

const Permission = ({roles, editable, onSave, permission = {}, create}) => {
  const {user_email, user_name, role, affiliation, privileged} = permission;
  const [revision, updateRevision] = useState({});

  const update = useCallback(
    (key, value) => {
      let newRevision = {...revision};
      if (permission[key] == value) {
        delete newRevision[key];
      } else {
        newRevision[key] = value;
      }
      updateRevision(newRevision);
    },
    [permission, revision]
  );

  const classes = useStyles();

  const handleOnSave = useCallback(() => {
    onSave({revision, user_email});
    updateRevision({});
  }, [revision, user_email]);

  const handleOnCancel = useCallback(() => {
    updateRevision({});
  });

  if (role == 'guest') roles = [...roles, 'guest'];

  return (
    <TableRow>
      <TableCell>
        {create ? (
          <TextField
            placeholder='New User Name'
            value={revision.name || ''}
            onChange={(e) => update('name', e.target.value)}
          />
        ) : (
          user_name
        )}
      </TableCell>
      <TableCell>
        {create ? (
          <TextField
            placeholder='Email'
            value={revision.email || ''}
            onChange={(e) => update('email', e.target.value)}
          />
        ) : (
          user_email
        )}
      </TableCell>
      <TableCell>
        {editable ? (
          <Select
            className={classes.role}
            value={revision.role || role || ''}
            onChange={(e) => update('role', e.target.value)}
          >
            {roles.map((r) => (
              <MenuItem key={r} value={r}>
                {r}
              </MenuItem>
            ))}
          </Select>
        ) : (
          <span className='role'>{role}</span>
        )}
      </TableCell>
      <TableCell>
        {editable ? (
          <TextField
            value={revision.affiliation || affiliation || ''}
            onChange={(e) => update('affiliation', e.target.value)}
          />
        ) : (
          <span className='affiliation'>{affiliation}</span>
        )}
      </TableCell>
      <TableCell padding='checkbox'>
        {editable && !create ? (
          <Checkbox
            checked={
              'privileged' in revision ? revision.privileged : privileged
            }
            onChange={(e) => update('privileged', e.target.checked)}
            inputProps={{'aria-label': 'privileged user'}}
          />
        ) : (
          privileged && 'Yes'
        )}
      </TableCell>
      <TableCell>
        {Object.keys(revision).length > 0 && (
          <SaveCancel onSave={handleOnSave} onCancel={handleOnCancel} />
        )}
      </TableCell>
    </TableRow>
  );
};

const userCount = (count, txt) =>
  count == 0 ? null : `${count} ${txt}${count === 1 ? '' : 's'}`;

const displayPermissions = (permissions, filter) =>
  permissions
    .sort((a, b) =>
      (a.role + a.user_email).localeCompare(b.role + b.user_email)
    )
    .filter(
      (p) =>
        !filter ||
        ['user_name', 'user_email', 'role', 'affiliation'].some(
          (column) =>
            column in p &&
            null != p[column] &&
            p[column].match(new RegExp(filter, 'i'))
        )
    );

const settingStyles = makeStyles((theme) => ({
  setting: {
    padding: '10px 15px'
  }
}));

const Setting = ({title, children}) => {
  const classes = settingStyles();

  return (
    <Grid alignItems='center' className={classes.setting} item container>
      <Grid item xs={2}>
        {title}
      </Grid>
      <Grid item xs={10}>
        {children}
      </Grid>
    </Grid>
  );
};

const ProjectView = ({project_name}) => {
  const classes = useStyles();

  const [projectType, setProjectType] = useState('');
  const [projectCoc, setProjectCoc] = useState('ff');
  const [projectContact, setProjectContact] = useState('');
  const [error, setError] = useState(null);
  const invoke = useActionInvoker();

  let user = useReduxState((state) => selectUser(state));
  let [project, setProject] = useState({});
  let [filter, setFilter] = useState(null);

  const setProjectDetails = (newProject) => {
    setProjectType(projectTypeFor(newProject));
    setProjectCoc(newProject.cc_text);
    setProjectContact(newProject.contact_email);
  };

  let retrieveProject = useCallback(() => {
    json_get(`/api/admin/${project_name}/info`).then(({project}) => {
      setProject(project);
      setProjectDetails(project);
    });
  }, [project_name, setProject]);
  useEffect(() => retrieveProject(), []);

  const updateProject = () => {
    postUpdateProject(project_name, {
      requires_agreement: projectType == 'community',
      resource: projectType == 'community' || projectType == 'resource',
      ...(projectCoc != project.cc_text && {cc_text: projectCoc}),
      ...(projectContact != project.contact_email && {
        contact_email: projectContact
      }),
      ...(projectType != 'community' && {cc_text: '', contact_email: ''})
    })
      .then((project) => {
        setProject(project);
        setProjectDetails(project);
        setError(null);
      })
      .catch((e) => e.then(({error}) => setError(error)));
  };

  const resetProject = () => {
    setProjectDetails(project);
    setError(null);
  };

  let {permissions = [], project_name_full} = project;

  let privileged = permissions.filter((p) => p.privileged);

  let editable = isAdmin(user, project_name);
  let roles = isSuperuser(user, project_name)
    ? ['administrator', 'editor', 'viewer', 'disabled']
    : isAdmin(user, project_name)
    ? ['editor', 'viewer', 'disabled']
    : [];

  const updateProjectType = (type) => {
    if (type == 'community')
      setError(
        "WARNING! Setting this project to 'community' will allow any library user to add themselves as a guest to the project."
      );
    else if (type == 'resource')
      setError(
        "WARNING! Setting this project to 'resource' will allow any library user to view the project data."
      );
    else setError(null);
    setProjectType(type);
  };

  const isChanged =
    projectType != projectTypeFor(project) ||
    projectCoc != project.cc_text ||
    projectContact != project.contact_email;

  const extensions = useMemo(
    () => [
      basicSetup,
      syntaxHighlighting(defaultHighlightStyle, {fallback: true}),
      markdown(),
      EditorView.editable.of(true),
      EditorState.readOnly.of(false),
      EditorView.lineWrapping,
      EditorState.tabSize.of(2),
      gutter({class: 'CodeMirror-lint-markers'}),
      lineNumbers(),
      lintGutter(),
      EditorView.updateListener.of(function (e) {
        if (e.docChanged) {
          setProjectCoc(e.state.doc.toString());
        }
      })
    ],
    []
  );

  return (
    <div id='project-view'>
      <Grid container direction='column'>
        <Grid className={classes.title}> {project_name_full}</Grid>
        <Grid container direction='column'>
          <Grid className={classes.subheading}>
            <Typography>Settings</Typography>
          </Grid>
          <Grid container direction='column'>
            <Setting title='Project Type'>
              <Select
                value={projectType}
                onChange={(e) => updateProjectType(e.target.value)}
              >
                {['team', 'community', 'resource'].map((r) => (
                  <MenuItem key={r} value={r}>
                    {r}
                  </MenuItem>
                ))}
              </Select>
            </Setting>
            {projectType == 'community' && (
              <>
                <Setting title='Code of Conduct Agreement'>
                  <Grid className={classes.editor}>
                    <ErrorBoundary>
                      <CodeMirror extensions={extensions} value={projectCoc} />
                    </ErrorBoundary>
                  </Grid>
                </Setting>
                <Setting title='Contact email'>
                  <TextField
                    placeholder='Email'
                    value={projectContact}
                    onChange={(e) => setProjectContact(e.target.value)}
                  />
                </Setting>
              </>
            )}
            {error && (
              <Setting>
                <Typography color='error'>{error}</Typography>
              </Setting>
            )}
            {isChanged && (
              <Setting>
                <SaveCancel onSave={updateProject} onCancel={resetProject} />
              </Setting>
            )}
          </Grid>
        </Grid>
        <Grid container direction='column'>
          <Grid className={classes.subheading}>
            <Typography>Users</Typography>
          </Grid>
          <Grid className={classes.summary}>
            {['administrator', 'editor', 'viewer', 'guest']
              .map((role) => {
                let ps = permissions.filter((p) => p.role == role);
                return userCount(ps.length, role);
              })
              .concat(userCount(privileged.length, 'privileged user'))
              .filter((_) => _)
              .join(', ')}
          </Grid>

          <TableContainer className={classes.table} component={Paper}>
            <Table aria-label='project users'>
              <TableHead>
                <TableRow>
                  <TableCell>Name</TableCell>
                  <TableCell>Email</TableCell>
                  <TableCell>Role</TableCell>
                  <TableCell>Affiliation</TableCell>
                  <TableCell padding='checkbox'>Privileged</TableCell>
                  <TableCell width='80'></TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {permissions.length > 1 && (
                  <TableRow>
                    <TableCell colspan='6'>
                      <TextField
                        placeholder='Filter rows'
                        onChange={(e) => setFilter(e.target.value)}
                      />
                    </TableCell>
                  </TableRow>
                )}
                {displayPermissions(permissions, filter).map((p) => (
                  <Permission
                    key={p.user_email}
                    permission={p}
                    editable={editable}
                    roles={roles}
                    onSave={({revision, user_email}) =>
                      postUpdatePermission(project_name, user_email, revision)
                        .then(() => retrieveProject())
                        .catch((e) =>
                          e.then(({error}) => invoke(showMessages([error])))
                        )
                    }
                  />
                ))}
                {editable && (
                  <Permission
                    create={true}
                    editable={true}
                    roles={['viewer', 'editor']}
                    onSave={({revision}) =>
                      postAddUser(project_name, revision)
                        .then(() => retrieveProject())
                        .catch((e) =>
                          e.then(({error}) => invoke(showMessages([error])))
                        )
                    }
                  />
                )}
              </TableBody>
            </Table>
          </TableContainer>
        </Grid>
      </Grid>
    </div>
  );
};
export default ProjectView;
