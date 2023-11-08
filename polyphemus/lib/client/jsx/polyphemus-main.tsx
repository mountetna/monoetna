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
import { ReactElement } from 'react';

const useStyles = makeStyles((theme) => ({
  etls: {
    minWidth: '800px'
  },
  etl_list: {
    border: '1px solid #ccc',
    height: 'calc(100vh - 127px)'
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

  useEffect( () => {
    const updateEtls = () => {
      json_get(`/api/etl/${project_name}/configs`).then(
        new_etls => {
          const newEtls = etls.map(
            etl => {
              const new_etl = new_etls.find( (e:Etl) => e.config_id == etl.config_id );
              if (new_etl) {
                const { ran_at, run_interval, status } = new_etl;

                etl = { ...etl, ran_at, run_interval, status };
              }

              return etl
            }
          );

          setEtls(newEtls);
        }
      );
    };

    const interval = setInterval(updateEtls, 30000);

    return () => clearInterval(interval);
  }, [etls]);

  useEffect(() => {
    json_get('/api/etl/jobs').then(setJobs);
    json_get(`/api/etl/${project_name}/configs`).then(setEtls);

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

  const selectedEtl = etls.find( e => e.config_id === selected );

  const headCells = [
    {
      id: 'Job Type',
      align: 'left' as const,
      key: 'etl'
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
      compare: (a:Etl,b:Etl) => runTime(a.ran_at,a.run_interval).localeCompare(runTime(b.ran_at,b.run_interval))
    },
    {
      id: 'Run State',
      align: 'right' as const,
      compare: (a:Etl,b:Etl) => a.run_interval - b.run_interval
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

  let add_button: ReactElement = null
  let primary_contents: ReactElement = <></>
  if (selected && selectedEtl) {
    primary_contents = <EtlConfig
      {...selectedEtl}
      onUpdate={addEtl}
      job={jobs.find((j) => j.name == selectedEtl?.etl)}
    />
  } else {
    add_button = <Grid item>
      <Button disabled={selected && selectedEtl} onClick={() => setCreate(true)}>Add Loader</Button>
      <EtlCreate
        project_name={project_name}
        open={create}
        onClose={() => setCreate(false)}
        onCreate={addEtl}
        jobs={jobs}
      />
    </Grid>
    primary_contents = <TableContainer className={classes.etl_list}>
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
                onClick={ () => setSelected(etl.config_id) }
                job={jobs.find((j) => j.name == etl.etl)}
              />
            ))}
        </TableBody>
      </Table>
    </TableContainer>
  }

  return (
    <Grid id='polyphemus-main'>
      {!jobs ? null : (
        <Grid className={classes.etls} item xs={12}>
          <Grid item container direction alignItems="center">
            <Grid item xs={5}>
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
                  selected ? <Typography>{ selectedEtl?.name }</Typography> : null
                }
              </Breadcrumbs>
            </Grid>
            {add_button}
          </Grid>
          {primary_contents}
        </Grid>
      )}
    </Grid>
  );
};

export default PolyphemusMain;
