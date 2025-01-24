import React, {useState, useEffect, useCallback, useContext} from 'react';
import {json_get} from 'etna-js/utils/fetch';
import {getModels} from 'etna-js/api/magma_api';
import {MagmaContext} from 'etna-js/contexts/magma-context';

import Typography from '@material-ui/core/Typography';
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

import { WorkflowConfig } from './workflow/workflow-config';
import WorkflowStatusRow from './workflow/workflow-status';
import WorkflowCreate from './workflow/workflow-create';
import {Workflow, Job, Status} from './polyphemus';
import {runTime} from './workflow/run-state';

const useStyles = makeStyles((theme) => ({
  workflows: {
    minWidth: '800px'
  },
  workflow_list: {
    border: '1px solid #ccc',
    height: 'calc(100vh - 125px)'
  },
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
  link: {
    color: theme.palette.primary.main,
    cursor: 'pointer'
  }
}));

const PolyphemusMain = ({project_name}: {project_name: string}) => {
  const [workflows, setWorkflows] = useState<Workflow[]>([]);
  const [runtimeConfigs, setRuntimeConfigs] = useState<{ [ config_id: string ]: RuntimeConfig }>({});
  const [jobs, setJobs] = useState<Job[] | null>(null);
  const [statuses, setStatuses] = useState<Status[]>([]);
  const [create, setCreate] = useState(false);
  const {models, setModels} = useContext(MagmaContext);

  const addWorkflow = (workflow: Workflow) => {
    let index = workflows.findIndex((e) => e.config_id == workflow.config_id);
    let new_workflows =
      index == -1
        ? workflows.concat(workflow)
        : workflow.run_interval == -2
        ? workflows.filter((e, i) => i != index)
        : workflows.map((e, i) => (i == index ? workflow : e));

    setWorkflows(new_workflows);
  };

  const classes = useStyles();

  const updateStatus = () => json_get(`/api/workflows/${project_name}/status`).then(setStatuses);

  useEffect( () => {
    const interval = setInterval(updateStatus, 30000);

    return () => clearInterval(interval);
  }, []);

  useEffect(() => {
    json_get('/api/workflows').then(setJobs);
    updateStatus();
    getModels(project_name).then(({models}) => setModels(models));
  }, []);

  const [order, setOrder] = useState< 'desc' | 'asc' | undefined>('asc');
  const [orderBy, setOrderBy] = useState('Job Type');
  const [selectedWorkflow, setSelectedWorkflow] = useState<Workflow|null>(null);
  const [ loading, setLoading ] = useState(false);

  const selectWorkflow = (config_id) => {
    setLoading(true);
    json_get(`/api/workflows/${project_name}/configs/${config_id}`).then(
      workflow => {
        setSelectedWorkflow(workflow);
        setLoading(false);
      }
    );
  }

  const headCells = [
    {
      id: 'Job Type',
      align: 'left' as const,
      key: 'workflow_type'
    },
    {
      id: 'Name',
      align: 'left' as const,
      key: 'workflow_name'
    },
    {
      id: 'Last Status',
      align: 'right' as const,
      key: 'pipeline_state'
    },
    {
      id: 'Last Completed',
      align: 'right' as const,
      key: 'pipeline_finished_at'
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
        ? (orderCell.compare
          ? orderCell.compare(a, b)
          : (a[ orderCell.key ] || '').localeCompare(b[ orderCell.key ] || ''))
        : 0) || a.workflow_name.localeCompare(b.workflow_name);
      return dir * comp;
    }, [ orderBy, order ]
  );

  return (
    <Grid id='polyphemus-main'>
      {!jobs ? null : (
        <Grid className={classes.workflows} item xs={12}>
          <Grid className={classes.title}>
            <Breadcrumbs className={classes.breadcrumb}>
              <Typography>
                {project_name}
              </Typography>
              {
                selectedWorkflow
                ? <Typography className={classes.link} onClick={ () => selectWorkflow(null) }>data loaders</Typography>
                : <Typography>data loaders</Typography>
              }
              {
                selectedWorkflow ? <Typography>{ selectedWorkflow?.name }</Typography> : null
              }
            </Breadcrumbs>
            <Button onClick={() => setCreate(true)}>Add Loader</Button>
            <WorkflowCreate
              project_name={project_name}
              open={create}
              onClose={() => setCreate(false)}
              onCreate={addWorkflow}
              jobs={jobs}
            />
          </Grid>
          {
            loading 
            ? <CircularProgress size={24}/>
            : selectedWorkflow
            ? <WorkflowConfig
                workflow={selectedWorkflow}
                status={ statuses.find(status => status.config_id == selectedWorkflow.config_id) }
                onUpdate={addWorkflow}
                job={jobs.find((j) => j.name == selectedWorkflow.workflow_type)}
              />
            : <>
              <TableContainer className={classes.workflow_list}>
                <Table size="small">
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
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {statuses
                      .sort((a, b) => compareBy(a, b))
                      .map((status: Status) => (
                        <WorkflowStatusRow
                          key={status.workflow_name}
                          status={status}
                          onClick={ () => selectWorkflow(status.config_id) }
                        />
                      ))}
                  </TableBody>
                </Table>
              </TableContainer>
            </>
          }
        </Grid>
      )}
    </Grid>
  );
};

export default PolyphemusMain;
