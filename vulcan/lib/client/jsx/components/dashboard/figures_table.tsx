import React, {
  useCallback,
  useContext,
  useState,
  useEffect,
  useMemo
} from 'react';
import 'regenerator-runtime/runtime';
import moment from 'moment-timezone';

import {
  DataGrid,
  GridCellParams,
  GridColumnHeaderParams,
  GridAlignment,
  GridValueFormatterParams
} from '@material-ui/data-grid';
import {makeStyles} from '@material-ui/core/styles';
import Chip from '@material-ui/core/Chip';

import {pushLocation} from 'etna-js/actions/location_actions';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import FigureMenu from './figure_menu';
import {VulcanContext} from '../../contexts/vulcan_context';
import ImageMemo from './image_memo';
import {VulcanFigureSession} from '../../api_types';
import {workflowByName} from '../../selectors/workflow_selectors';

const figureListStyles = makeStyles((theme) => ({
  title: {
    padding: '10px 0px',
    color: '#444'
  },
  figures: {
    padding: '15px'
  },
  thumbnail: {
    maxWidth: '32px',
    maxHeight: '32px'
  }
}));

export default function FiguresTable({
  project_name,
  workflowName
}: {
  project_name: string;
  workflowName?: string;
}) {
  const {
    showErrors,
    fetchFigures,
    createFigure,
    updateFigure,
    deleteFigure,
    state
  } = useContext(VulcanContext);

  const [allFigureSessions, setAllFigureSessions] = useState<
    VulcanFigureSession[]
  >([]);
  const [filteredFigureSessions, setFilteredFigureSessions] = useState<
    VulcanFigureSession[]
  >([]);

  const invoke = useActionInvoker();
  const classes = figureListStyles();

  useEffect(() => {
    showErrors(
      fetchFigures(
        project_name
      ).then(({figures}: {figures: VulcanFigureSession[]}) =>
        setAllFigureSessions(figures)
      )
    );
  }, [showErrors, fetchFigures, project_name]);

  useEffect(() => {
    if (workflowName) {
      setFilteredFigureSessions(
        allFigureSessions.filter(
          (figure) => figure.workflow_name === workflowName
        )
      );
    } else {
      setFilteredFigureSessions([...allFigureSessions]);
    }
  }, [workflowName, allFigureSessions]);

  const handleOnCopy = useCallback(
    (figure: VulcanFigureSession) => {
      const copy = {
        ...figure,
        figure_id: null,
        title: `${figure.title} - copy`
      };
      showErrors(
        createFigure(project_name, copy).then((newFigure) => {
          setAllFigureSessions([...allFigureSessions].concat([newFigure]));
        })
      );
    },
    [showErrors, createFigure, project_name, allFigureSessions]
  );

  const handleOnRename = useCallback(
    (figure: VulcanFigureSession) => {
      const newTitle = prompt(
        'Please enter a new figure title',
        figure.title || ''
      );

      if (!newTitle) return;

      showErrors(
        updateFigure(project_name, {
          ...figure,
          title: newTitle
        }).then((updatedFigure) => {
          const updated = allFigureSessions.map((oldFigure) => {
            if (oldFigure.figure_id === updatedFigure.figure_id) {
              return updatedFigure;
            }

            return oldFigure;
          });

          setAllFigureSessions(updated);
        })
      );
    },
    [showErrors, updateFigure, project_name, allFigureSessions]
  );

  const handleOnRemove = useCallback(
    (figure: VulcanFigureSession) => {
      if (!figure.figure_id) return;
      showErrors(
        deleteFigure(project_name, figure.figure_id).then(() => {
          const updated = allFigureSessions.filter((oldFigure) => {
            return oldFigure.figure_id !== figure.figure_id;
          });
          setAllFigureSessions(updated);
        })
      );
    },
    [showErrors, deleteFigure, project_name, allFigureSessions]
  );

  const visitFigure = useCallback(
    (figureId: number) => {
      invoke(pushLocation(`/${project_name}/figure/${figureId}`));
    },
    [invoke]
  );

  const columns = [
    {
      field: 'workflow',
      flex: 1,
      renderCell: (params: GridCellParams) => {
        const workflow = workflowByName(params.value as string, state);
        return (
          <ImageMemo
            className={classes.thumbnail}
            src={`/images/${workflow?.image || 'default.png'}`}
            alt='Figure thumbnail'
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
      field: 'title',
      headerName: 'Title',
      flex: 2
    },
    {
      field: 'author',
      headerName: 'Author',
      flex: 2
    },
    {
      field: 'workflow_name',
      headerName: 'Workflow',
      flex: 2,
      valueFormatter: (params: GridValueFormatterParams) => {
        const workflow = workflowByName(params.value as string, state);
        return workflow?.displayName || workflow?.name;
      }
    },
    {
      field: 'updated_at',
      headerName: 'Last Modified',
      flex: 2,
      valueFormatter: (params: GridValueFormatterParams) => {
        return moment
          .utc(params.value as string)
          .tz(moment.tz.guess())
          .format('MMM D, YYYY');
      }
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
        <FigureMenu
          figureId={params.id as number}
          onCopy={() => handleOnCopy(params.row as VulcanFigureSession)}
          onRename={() => handleOnRename(params.row as VulcanFigureSession)}
          onRemove={() => handleOnRemove(params.row as VulcanFigureSession)}
        />
      ),
      align: 'center' as GridAlignment,
      hideSortIcons: true,
      disableColumnMenu: true,
      resizable: false,
      editable: false
    }
  ];

  const rows = useMemo(() => {
    return filteredFigureSessions.sort((a, b) => {
      if (null == a.figure_id) return -1;
      if (null == b.figure_id) return 1;
      return a.figure_id - b.figure_id;
    });
  }, [filteredFigureSessions]);

  return (
    <DataGrid
      autoHeight={true}
      rows={rows}
      columns={columns}
      getRowId={(row) => row.figure_id}
      pageSize={20}
      onRowClick={(params, e) => {
        visitFigure(params.id as number);
      }}
      hideFooterSelectedRowCount={true}
    />
  );
}
