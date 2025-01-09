import React, {useState, useEffect, useCallback, useContext} from 'react';
import {json_get} from 'etna-js/utils/fetch';
import {getDocuments} from 'etna-js/api/magma_api';

import Typography from '@material-ui/core/Typography';
import Breadcrumbs from '@material-ui/core/Breadcrumbs';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import { WorkflowConfigRow, WorkflowConfig } from './workflow/workflow-config';
import WorkflowCreate from './workflow/workflow-create';
import {Workflow, Job} from './polyphemus';

import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableSortLabel from '@material-ui/core/TableSortLabel';
import TableRow from '@material-ui/core/TableRow';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import Paper from '@material-ui/core/Paper';

import {runTime} from './workflow/run-state';
import {MagmaContext} from 'etna-js/contexts/magma-context';

const useStyles = makeStyles((theme) => ({
  workflows: {
    minWidth: '800px'
  },
  workflow_list: {
    border: '1px solid #ccc',
    height: 'calc(100vh - 160px)',
    marginBottom: '5px'
  },
  title: {
    padding: '15px',
    color: 'black'
  },
  link: {
    color: theme.palette.primary.main,
    cursor: 'pointer'
  }
}));

const PolyphemusMain = ({project_name}: {project_name: string}) => {
  const [workflows, setWorkflows] = useState<Workflow[]>([]);
  const [jobs, setJobs] = useState<Job[] | null>(null);
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

  useEffect( () => {
    const updateWorkflows = () => {
      json_get(`/api/workflow/${project_name}/configs`).then(
        new_workflows => {
          const newWorkflows = workflows.map(
            workflow => {
              const new_workflow = new_workflows.find( (e:Workflow) => e.config_id == workflow.config_id );
              if (new_workflow) {
                const { ran_at, run_interval, status } = new_workflow;

                workflow = { ...workflow, ran_at, run_interval, status };
              }

              return workflow
            }
          );

          setWorkflows(newWorkflows);
        }
      );
    };

    const interval = setInterval(updateWorkflows, 30000);

    return () => clearInterval(interval);
  }, [workflows]);

  useEffect(() => {
    json_get('/api/workflow/workflows').then(setJobs);
    json_get(`/api/workflow/${project_name}/configs`).then(setWorkflows);

    getDocuments(
      {
        project_name,
        model_name: 'all',
        record_names: [],
        attribute_names: 'all'
      },
      fetch
    )
      .then(({models}) => setModels(models))
      .catch((e) => console.log({e}));
  }, []);

  const collator = new Intl.Collator(undefined, {
    numeric: true,
    sensitivity: 'base'
  });


  const [order, setOrder] = useState< 'desc' | 'asc' | undefined>('asc');
  const [orderBy, setOrderBy] = useState('Job Type');
  const [selected, setSelected] = useState<number|null>(null);

  const selectedWorkflow = workflows.find( e => e.config_id === selected );

  const headCells = [
    {
      id: 'Job Type',
      align: 'left' as const,
      key: 'workflow'
    },
    {
      id: 'Name',
      align: 'left' as const,
      key: 'name'
    },
    {
      id: 'Last Status',
      align: 'right' as const,
      key: 'status'
    },
    {
      id: 'Last Ran',
      align: 'right' as const,
      key: 'ran_at'
    },
    {
      id: 'Next Run',
      align: 'right' as const,
      compare: (a:Workflow,b:Workflow) => runTime(a.ran_at,a.run_interval).localeCompare(runTime(b.ran_at,b.run_interval))
    },
    {
      id: 'Run State',
      align: 'right' as const,
      compare: (a:Workflow,b:Workflow) => a.run_interval - b.run_interval
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
        : 0) || a.name.localeCompare(b.name);
      return dir * comp;
    }, [ orderBy, order ]
  );

  return (
    <Grid id='polyphemus-main'>
      {!jobs ? null : (
        <Grid className={classes.workflows} item xs={12}>
          <Breadcrumbs className={classes.title}>
            <Typography>
              {project_name}
            </Typography>
            {
              selected
              ? <Typography className={classes.link} onClick={ () => setSelected(null) }>data loaders</Typography>
              : <Typography>data loaders</Typography>
            }
            {
              selected ? <Typography>{ selectedWorkflow?.name }</Typography> : null
            }
          </Breadcrumbs>
          {
            selected && selectedWorkflow
            ? <WorkflowConfig
                {...selectedWorkflow}
                onUpdate={addWorkflow}
                job={jobs.find((j) => j.name == selectedWorkflow?.workflow)}
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
                    {workflows
                      .sort((a, b) => compareBy(a, b))
                      .map((workflow: Workflow) => (
                        <WorkflowConfigRow
                          key={workflow.name}
                          {...workflow}
                          onClick={ () => setSelected(workflow.config_id) }
                          job={jobs.find((j) => j.name == workflow.workflow_type)}
                        />
                      ))}
                  </TableBody>
                </Table>
              </TableContainer>
              <Grid>
                <Button onClick={() => setCreate(true)}>Add Loader</Button>
                <WorkflowCreate
                  project_name={project_name}
                  open={create}
                  onClose={() => setCreate(false)}
                  onCreate={addWorkflow}
                  jobs={jobs}
                />
              </Grid>
            </>
          }
        </Grid>
      )}
    </Grid>
  );
};

export default PolyphemusMain;
