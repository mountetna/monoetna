import React, {useContext, useMemo} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';

import OutputUI from './drawers/output_user_input';
import { outputUIsWithInputsReady } from '../../selectors/workflow_selectors';
import {useWorkspace} from '../../contexts/workspace_context';
import LoadingIcon from '../dashboard/loading_icon';

export default function OutputFeed() {
  // Shows stream of Output, Plots, etc.,
  //   as the session object updates.
  const {state} = useContext(VulcanContext);
  const {workspace} = useWorkspace();
  const {status, update_files} = state;

  const outputFeed = useMemo(() => {
    return state.update_files ?
      <div>
        <LoadingIcon/>
        Refreshing Files
      </div> :
      outputUIsWithInputsReady(workspace, status).map((s, index) => (
        <OutputUI key={index} step={s}/>
      ));
  }, [workspace, status, update_files]);

  return (
    <div className="session-output-feed">
      {outputFeed}
    </div>
  );
}
