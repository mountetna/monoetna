import React, {useState, useEffect, useCallback, useContext} from 'react';
import {json_post} from 'etna-js/utils/fetch';
import {getModels} from 'etna-js/api/magma_api';
import {MagmaContext} from 'etna-js/contexts/magma-context';

import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import FormControl from '@material-ui/core/FormControl';
import InputLabel from '@material-ui/core/InputLabel';
import Input from '@material-ui/core/Input';
import Breadcrumbs from '@material-ui/core/Breadcrumbs';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import CircularProgress from '@material-ui/core/CircularProgress';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableSortLabel from '@material-ui/core/TableSortLabel';
import TableRow from '@material-ui/core/TableRow';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import Paper from '@material-ui/core/Paper';
import Chip from '@material-ui/core/Chip';
import CodeIcon from '@material-ui/icons/Code';
import IconButton from '@material-ui/core/IconButton';
import Collapse from '@material-ui/core/Collapse';

import {formatTime} from './workflow/run-state';
import {Log} from './polyphemus';
import { getItem } from 'etna-js/utils/cookies';
import { parseToken } from 'etna-js/utils/janus';

declare var CONFIG: { [key: string]: string };


const useStyles = makeStyles((theme) => ({
  title: {
    padding: '10px 0px',
    display: 'flex',
    alignItems: 'center',
    color: 'black'
  },
  breadcrumb: {
    padding: '0px 10px',
    flex: '1 1 auto'
  },
  log_list: {
    border: '1px solid #ccc',
    height: 'calc(100vh - 125px)'
  },
  filters: {
    borderBottom: '1px solid #ccc',
    padding: '10px',
    height: '70px',
    alignItems: 'end'
  },
  filter: {
    marginRight: '10px'
  },
  app_filter: {
    minWidth: '120px'
  },
  janus: {
    background: 'rgba(155,86,255,0.33)'
  },
  magma: {
    background: '#cd4a34'
  },
  timur: {
    background: 'rgba(117,255,117,0.30)'
  },
  metis: {
    background: '#5fe2e4'
  },
  gnomon: {
    background: '#e6e6e6'
  },
  polyphemus: {
    background: '#d0d835'
  },
  vulcan: {
    background: 'rgba(255,8,4,0.49)'
  }
}));

const applications = [
  'janus',
  'magma',
  'timur',
  'metis',
  'gnomon',
  'polyphemus',
  'vulcan'
];

const LogRow = ({log, privileged}:{log: Log, privileged: boolean}) => {
  const [ showPayload, setShowPayload ] = useState(false);
  const [ payload, setPayload ] = useState(false);
  const loadPayload = useCallback( () => {
    if (showPayload) {
      setShowPayload(false);
      return;
    }
    if (!payload) {
      json_post(`/api/log/${log.project_name}/payload/${log.id}`).then(
        ({payload}) => setPayload(payload)
      );
    }
    setShowPayload(true);
  }, [payload, showPayload] );
  const classes: any = useStyles();
  return <>
    <TableRow>
      <TableCell align="right">{formatTime(log.created_at)}</TableCell>
      <TableCell> <Chip label={log.application} className={ classes[log.application]} /> </TableCell>
      <TableCell>{log.project_name}</TableCell>
      <TableCell>{log.user}</TableCell>
      <TableCell>{log.event}</TableCell>
      <TableCell>{log.message}</TableCell>
      { privileged && <TableCell align='right'>{
        log.payload && <IconButton onClick={ loadPayload } size="small"><CodeIcon fontSize="small"/></IconButton>
      }</TableCell> }
    </TableRow>
    { privileged && <TableRow>
        <TableCell style={{ paddingBottom: 0, paddingTop: 0 }} colSpan={7}>
          <Collapse in={showPayload} timeout="auto" unmountOnExit>
          <pre>{ JSON.stringify(payload, null, 2) }</pre>
          </Collapse>
        </TableCell>
      </TableRow>
    }
  </>
}

const PolyphemusLogs = ({project_name}: {project_name: string}) => {
  const classes: any = useStyles();
  const [order, setOrder] = useState< 'desc' | 'asc' | undefined>('asc');
  const [orderBy, setOrderBy] = useState('Job Type');
  const [logs, setLogs] = useState<Log[]>([]);

  const [ filterApplications, setFilterApplications ] = useState<string[]>([]);
  const [ filterUser, setFilterUser ] = useState('');
  const [ filterEvent, setFilterEvent ] = useState('');
  const [ filterMessage, setFilterMessage ] = useState('');
  const [ filterProjects, setFilterProjects ] = useState<string[]>([]);
  const [ filterFrom, setFilterFrom ] = useState('');
  const [ filterTo, setFilterTo ] = useState('');

  const token = parseToken(getItem(CONFIG.token_name) as string);

  const isAdmin = token.permissions['administration'] && (
    token.permissions.administration.role == 'administrator' ||
    token.permissions.administration.role == 'editor');
  const isPrivileged = isAdmin || token.permissions[project_name] && (
    token.permissions[project_name].privileged
  );

  const getLogs = useCallback(() => {
    let params: any = {};

    if (filterUser) params.user = filterUser;
    if (filterEvent) params.event = filterEvent;
    if (filterMessage) params.message = filterMessage;
    if (filterApplications.length) params.application = filterApplications;
    if (filterTo && filterFrom) {
      params.from = filterFrom;
      params.to = filterTo;
    }
    if (filterProjects.length && project_name == 'administration') {
      params.project_names = filterProjects;
    }
    json_post(
      `/api/log/${project_name}/read`,
      params
    ).then(({logs}) => setLogs(logs));
  }, [ filterApplications, filterUser, filterEvent, filterMessage, filterFrom, filterTo, filterProjects ]);

  useEffect(() => {
    getLogs();
  }, []);

  useEffect(() => {
    const timeout = setInterval(getLogs, 10000);

    return () => clearInterval(timeout);
  }, [getLogs]);

  const headCells = [
    {
      id: 'Time',
      align: 'right' as const,
      key: 'created_at'
    },
    {
      id: 'Application',
      align: 'left' as const,
      key: 'application'
    },
    {
      id: 'Project',
      align: 'left' as const,
      key: 'project_name'
    },
    {
      id: 'User',
      align: 'left' as const,
      key: 'user'
    },
    {
      id: 'Event',
      align: 'left' as const,
      key: 'event'
    },
    {
      id: 'Message',
      align: 'left' as const,
      key: 'message'
    }
  ];

  const orderCell = headCells.find( cell => cell.id === orderBy );

  const sortBy = useCallback(
    (id) => {
      const isAsc = orderBy === id && order === 'asc';
      setOrder(isAsc ? 'desc' : 'asc');
      setOrderBy(id);
    }, [ orderBy, order ]
  );

  const compareBy = useCallback(
    (a, b) => {
      const dir = (order == 'asc' ? 1 : -1);
      const comp = (orderCell
        ? (a[ orderCell.key ] || '').localeCompare(b[ orderCell.key ] || '')
        : 0) || -1 * a.created_at.localeCompare(b.created_at);
      return dir * comp;
    }, [ orderBy, order ]
  );

  return <Grid id='polyphemus-main'>
    <Grid item xs={12}>
      <Grid className={classes.title}>
        <Breadcrumbs className={classes.breadcrumb}>
          <Typography>
            {project_name}
          </Typography>
          <Typography>logs</Typography>
        </Breadcrumbs>
      </Grid>
      <TableContainer className={ classes.log_list } >
        <Grid container className={ classes.filters }>
          <TextField className={ classes.filter } placeholder='User' value={ filterUser } onChange={ e => setFilterUser(e.target.value) }/>
          <TextField className={ classes.filter } placeholder='Event' value={ filterEvent } onChange={ e => setFilterEvent(e.target.value) }/>
          <TextField className={ classes.filter } placeholder='Message' value={ filterMessage } onChange={ e => setFilterMessage(e.target.value) }/>
          {
            project_name == 'administration' && <TextField className={ classes.filter } placeholder='Projects (comma-separated)' value={ filterProjects } onChange={ e => setFilterProjects(e.target.value.split(',').filter(_=>_)) }/>
          }
          <TextField className={classes.filter}
            label="From"
            type="date"
            InputLabelProps={{ shrink: true }}
            value={filterFrom}
            onChange={ e => setFilterFrom(e.target.value) }
          />
          <TextField className={classes.filter}
            label="To"
            type="date"
            InputLabelProps={{ shrink: true }}
            value={filterTo}
            onChange={ e => setFilterTo(e.target.value) }
          />
          <FormControl className={classes.app_filter}>
            <InputLabel id="app-label">Applications</InputLabel>
            <Select
              labelId='app-label'
              multiple
              value={filterApplications}
              onChange={ e => setFilterApplications(e.target.value as string[]) }
              input={<Input />}
              renderValue={
                (selected) => <Grid>
                  {
                    (selected as string[]).map(
                      (value: string) => <Chip key={value} label={value} className={classes[value]} />
                    )
                  }
                </Grid>
              }
            >
              {
                applications.map(app => <MenuItem key={app} value={app}>{app}</MenuItem>)
              }
            </Select>
          </FormControl>
        </Grid>
        <Table size="small" stickyHeader>
          <TableHead>
            <TableRow>
              {headCells.map((cell) =>
              <TableCell key={cell.id}
                align={cell.align}
                sortDirection={ orderBy === cell.id ? order : 'asc' }>
                <TableSortLabel
                  active={ orderBy === cell.id }
                  direction={ orderBy === cell.id ? order : 'asc' }
                  onClick={ () => sortBy(cell.id) }>
                  {cell.id}
                </TableSortLabel>
              </TableCell>
              )}
              {
                isPrivileged && <TableCell/>
              }
            </TableRow>
          </TableHead>
          <TableBody>
            {logs
              .sort((a, b) => compareBy(a, b))
              .map((log: Log) => (
                <LogRow
                  key={log.id}
                  log={log}
                  privileged={isPrivileged}
                />
              ))}
          </TableBody>
        </Table>
      </TableContainer>
    </Grid>
  </Grid>
};

export default PolyphemusLogs;
