import React, {useContext, useMemo} from 'react';
import 'regenerator-runtime/runtime';

import Chip from '@material-ui/core/Chip';
import {
  DataGrid,
  GridCellParams,
  GridValueFormatterParams
} from '@material-ui/data-grid';
import {makeStyles} from '@material-ui/core/styles';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import {VulcanContext} from '../../contexts/vulcan_context';

const useStyles = makeStyles((theme) => ({
  title: {
    padding: '10px 15px 5px',
    color: '#444'
  },
  workflows: {
    padding: '15px'
  },
  none: {
    color: '#f44'
  }
}));

export default function WorkflowsTable({project_name}: {project_name: string}) {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  const classes = useStyles();

  const projectWorkflows = workflows
    ? workflows.filter(({projects}: {projects?: string[]}) =>
        projects?.includes(project_name)
      )
    : [];

  const columns = [
    {
      field: 'displayName',
      headerName: 'Name',
      flex: 1
    },
    {
      field: 'authors',
      headerName: 'Authors',
      valueFormatter: (params: GridValueFormatterParams) => {
        return ((params.value || []) as string[]).join(',');
      },
      flex: 1
    },
    {
      field: 'description',
      headerName: 'Description',
      flex: 3
    },
    {
      field: 'lastModified',
      headerName: 'Last Modified',
      flex: 1
    },
    {
      field: 'tags',
      headerName: 'Tags',
      renderCell: (params: GridCellParams) => {
        console.log('params', params);
        return (
          <>
            {((params.value || []) as string[]).map((tag: string) => {
              return <Chip label={tag} />;
            })}
          </>
        );
      },
      flex: 1
    }
  ];

  const rows = useMemo(() => {
    return projectWorkflows.sort((a, b) =>
      (a.displayName || a.name).localeCompare(b.displayName || b.name)
    );
  }, [projectWorkflows]);

  return (
    <DataGrid
      autoHeight={true}
      rows={rows}
      columns={columns}
      getRowId={(row) => row.name}
      pageSize={20}
      onRowClick={(params, e) => {}}
    />
  );
}
