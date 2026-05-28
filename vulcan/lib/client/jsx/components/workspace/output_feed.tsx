import React, {useContext, useMemo} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';

import OutputUI from './drawers/output_user_input';
import { outputUIsWithInputsReady } from '../../selectors/workflow_selectors';
import {WorkspaceContext} from '../../contexts/workspace_context';
import { LoadingIconWithText } from '../dashboard/loading_icon';

export default function OutputFeed() {
  const {state: { workspaceDetails, file_contents } } = useContext(WorkspaceContext);

  const outputFeed = useMemo(() => {
    return outputUIsWithInputsReady(workspaceDetails, file_contents).map((s, index) => (
        <OutputUI key={index} step={s}/>
      ));
  }, [workspaceDetails, file_contents]);

  return (
    <div className="session-output-feed">
      {outputFeed}
    </div>
  );
}
