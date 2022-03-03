import React, {useCallback, useContext, useMemo, useState} from 'react';
import * as _ from 'lodash';
import 'regenerator-runtime/runtime';

import Chip from '@material-ui/core/Chip';
import Tooltip from '@material-ui/core/Tooltip';
import IconButton from '@material-ui/core/IconButton';
import InsertChartIcon from '@material-ui/icons/InsertChart';
import EditIcon from '@material-ui/icons/Edit';
import {
  DataGrid,
  GridCellParams,
  GridValueFormatterParams,
  GridRowId,
  GridSelectionModel,
  GridColumnHeaderParams,
  GridAlignment
} from '@material-ui/data-grid';
import {makeStyles} from '@material-ui/core/styles';

import {pushLocation} from 'etna-js/actions/location_actions';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import {VulcanContext} from '../../contexts/vulcan_context';
import ImageMemo from './image_memo';

// To get webpack to pick up the files.
require('../../../img/umap.png');
require('../../../img/barplot.png');
require('../../../img/scatter.png');
require('../../../img/yplot.png');
require('../../../img/default.png');
require('../../../img/add_integers.png');

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
  },
  thumbnail: {
    maxWidth: '32px',
    maxHeight: '32px'
  }
}));

export default function WorkflowsTable({project_name}: {project_name: string}) {
  const [selectionModel, setSelectionModel] = useState<GridRowId[]>([]);
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  const classes = useStyles();

  const projectWorkflows = workflows
    ? workflows.filter(({projects}: {projects?: string[]}) =>
        projects?.includes(project_name)
      )
    : [];

  const handleCreateFigure = useCallback(
    (workflowName: string) => {
      invoke(
        pushLocation(
          `/${project_name}/figure/new/${workflowName.replace('.cwl', '')}`
        )
      );
    },
    [invoke, project_name]
  );

  const handleEditWorkflow = useCallback(() => {}, []);

  const columns = [
    {
      field: 'image',
      flex: 1,
      renderCell: (params: GridCellParams) => {
        return (
          <ImageMemo
            className={classes.thumbnail}
            src={`/images/${params.value || 'default.png'}`}
            alt='Workflow image'
          />
        );
      },
      renderHeader: (params: GridColumnHeaderParams) => <></>,
      hideSortIcons: true,
      disableColumnMenu: true,
      resizable: false,
      editable: false,
      align: 'center' as GridAlignment
    },
    {
      field: 'displayName',
      headerName: 'Name',
      flex: 2
    },
    {
      field: 'authors',
      headerName: 'Authors',
      valueFormatter: (params: GridValueFormatterParams) => {
        return ((params.value || []) as string[]).join(', ');
      },
      flex: 2
    },
    {
      field: 'description',
      headerName: 'Description',
      flex: 4
    },
    {
      field: 'lastModified',
      headerName: 'Last Modified',
      flex: 2
    },
    {
      field: 'tags',
      headerName: 'Tags',
      renderCell: (params: GridCellParams) => {
        return (
          <>
            {((params.value || []) as string[]).map(
              (tag: string, index: number) => {
                return <Chip key={index} label={tag} />;
              }
            )}
          </>
        );
      },
      flex: 2
    },
    {
      field: 'actions',
      type: 'actions',
      flex: 1,
      renderHeader: (params: GridColumnHeaderParams) => <></>,
      renderCell: (params: GridCellParams) => (
        <>
          <Tooltip title='Create figure'>
            <IconButton
              aria-label='Create figure'
              onClick={() => handleCreateFigure(params.id as string)}
              color='primary'
            >
              <InsertChartIcon />
            </IconButton>
          </Tooltip>
          {/* <IconButton
            aria-label='Edit workflow'
            disabled
            onClick={handleEditWorkflow}
            color='secondary'
          >
            <EditIcon />
          </IconButton> */}
        </>
      ),
      align: 'center' as GridAlignment,
      hideSortIcons: true,
      disableColumnMenu: true,
      resizable: false,
      editable: false
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
      selectionModel={selectionModel}
      columns={columns}
      getRowId={(row) => row.name}
      pageSize={20}
      onRowClick={(params, e) => {
        // Here, call prop to filter out saved figures?
      }}
      hideFooterSelectedRowCount={true}
      onSelectionModelChange={(model: GridSelectionModel) => {
        setSelectionModel(_.isEqual(selectionModel, model) ? [] : model);
      }}
    />
  );
}
