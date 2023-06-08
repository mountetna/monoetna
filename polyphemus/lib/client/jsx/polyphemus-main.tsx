import React, {useState, useEffect, useCallback, useContext} from 'react';
import {json_get} from 'etna-js/utils/fetch';
import {getDocuments} from 'etna-js/api/magma_api';

import Typography from '@material-ui/core/Typography';
import Breadcrumbs from '@material-ui/core/Breadcrumbs';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import { EtlConfigRow, EtlConfig } from './etl/etl-config';
import EtlCreate from './etl/etl-create';
import {Etl, Job} from './polyphemus';

import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableSortLabel from '@material-ui/core/TableSortLabel';
import TableRow from '@material-ui/core/TableRow';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import Paper from '@material-ui/core/Paper';

import {runTime} from './etl/run-state';
import {MagmaContext} from 'etna-js/contexts/magma-context';

const useStyles = makeStyles((theme) => ({
  etls: {
    minWidth: '800px'
  },
  etl_list: {
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
  const [etls, setEtls] = useState<Etl[]>([]);
  const [jobs, setJobs] = useState<Job[] | null>(null);
  const [create, setCreate] = useState(false);
  const {models, setModels} = useContext(MagmaContext);

  const addEtl = (etl: Etl) => {
    let index = etls.findIndex((e) => e.config_id == etl.config_id);
    let new_etls =
      index == -1
        ? etls.concat(etl)
        : etl.run_interval == -2
        ? etls.filter((e, i) => i != index)
        : etls.map((e, i) => (i == index ? etl : e));

    setEtls(new_etls);
  };

  const classes = useStyles();

  useEffect(() => {
    json_get(`/api/etl/${project_name}/configs`).then(setEtls);
    json_get('/api/etl/jobs').then(setJobs);
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


  const [order, setOrder] = React.useState('asc');
  const [orderBy, setOrderBy] = React.useState('Job Type');
  const [selected, setSelected] = React.useState<number|null>(null);

  const selectedEtl = etls.find( e => e.config_id === selected );

  const headCells = [
    {
      id: 'Job Type',
      align: 'left',
      key: 'etl'
    },
    {
      id: 'Name',
      align: 'left',
      key: 'name'
    },
    {
      id: 'Last Status',
      align: 'right',
      key: 'status'
    },
    {
      id: 'Last Ran',
      align: 'right',
      key: 'ran_at'
    },
    {
      id: 'Next Run',
      align: 'right',
      compare: (a,b) => runTime(a.ran_at,a.run_interval).localeCompare(runTime(b.ran_at,b.run_interval))
    },
    {
      id: 'Run State',
      align: 'right',
      compare: (a,b) => a.run_interval - b.run_interval
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
        <Grid className={classes.etls} item xs={12}>
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
              selected ? <Typography>{ selectedEtl.name }</Typography> : null
            }
          </Breadcrumbs>
          {
            selected
            ? <EtlConfig
                {...selectedEtl}
                onUpdate={addEtl}
                job={jobs.find((j) => j.name == selectedEtl.etl)}
              />
            : <>
              <TableContainer className={classes.etl_list}>
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
                    {etls
                      .sort((a, b) => compareBy(a, b))
                      .map((etl: Etl) => (
                        <EtlConfigRow
                          key={etl.name}
                          {...etl}
                          onUpdate={addEtl}
                          onClick={ () => setSelected(etl.config_id) }
                          job={jobs.find((j) => j.name == etl.etl)}
                        />
                      ))}
                  </TableBody>
                </Table>
              </TableContainer>
              <Grid>
                <Button onClick={() => setCreate(true)}>Add Loader</Button>
                <EtlCreate
                  project_name={project_name}
                  open={create}
                  onClose={() => setCreate(false)}
                  onCreate={addEtl}
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
