import React, {useContext} from 'react';
import {VulcanContext} from '../../../../contexts/vulcan_context';
import Link from 'etna-js/components/link';

export default function LinkOutput({data}) {
  let {state, vulcanPath} = useContext(VulcanContext);
  let {workspace, projectName} = state;
  if (!workspace) {
    return
  }
  return <React.Fragment>
    {Object.keys(data).map(file_name =>
        <React.Fragment key={file_name}>
          <Link link={vulcanPath(`/api/v2/${projectName}/workspace/${workspace.workspace_id}/file/download/${file_name}`)}>{file_name}</Link>
        </React.Fragment>
    )}
  </React.Fragment>;
}
