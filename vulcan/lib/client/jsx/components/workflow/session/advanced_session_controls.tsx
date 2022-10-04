import React, {useState, useContext, useCallback, useMemo} from 'react';

import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import IconButton from '@material-ui/core/IconButton';
import MoreVertIcon from '@material-ui/icons/MoreVert';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {VulcanFigure, VulcanSession} from '../../../api_types';
import {
  setSessionAndFigure,
  setWorkflow
} from '../../../actions/vulcan_actions';

export default function AdvancedSessionControls({
  session,
  figure
}: {
  figure: VulcanFigure;
  session: VulcanSession;
}) {
  const [menuAnchor, setMenuAnchor] = useState(
    null as HTMLButtonElement | null
  );
  const {dispatch, updateFigureDependencies} = useContext(VulcanContext);

  const handleClose = () => {
    setMenuAnchor(null);
  };

  const updateDependencies = useCallback(() => {
    // Post an update request for the figure,
    //   then set local workflow / session with
    //   the new reference_figure_id and snapshot.
    if (!figure.figure_id) return;

    updateFigureDependencies(session.project_name, figure.figure_id).then(
      (updatedFigure) => {
        dispatch(setSessionAndFigure(updatedFigure));
        if (updatedFigure.workflow_snapshot)
          {dispatch(
            setWorkflow(updatedFigure.workflow_snapshot, session.project_name)
          );}
        handleClose();
      }
    );
  }, [session, figure, updateFigureDependencies, dispatch]);

  const viewingRevision = useMemo(() => {
    return figure.id !== session.reference_figure_id;
  }, [figure, session]);

  return (
    <>
      <IconButton
        onClick={(e) => {
          e.stopPropagation();
          setMenuAnchor(e.currentTarget);
        }}
      >
        <MoreVertIcon />
      </IconButton>
      <Menu
        id={`advanced-session-menu-${session.key}`}
        open={Boolean(menuAnchor)}
        anchorEl={menuAnchor}
        onClose={handleClose}
      >
        <MenuItem onClick={updateDependencies} disabled={viewingRevision}>
          Update dependencies
        </MenuItem>
      </Menu>
    </>
  );
}
