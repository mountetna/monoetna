import React, {useContext} from 'react';
import {VulcanContext} from '../../../../contexts/vulcan_context';
import Link from 'etna-js/components/link';

export function LinkOutput({data}) {
  let {state, vulcanPath} = useContext(VulcanContext);
  let {workspace, projectName} = state;
  if (!workspace) {
    return
  }
  const target = `${projectName}/workspace/${workspace.workspace_id}/file/stream-download`
  return <React.Fragment>
    {Object.keys(data).map(file_name =>
        <React.Fragment key={file_name}>
          <Link link={vulcanPath(`/api/v2/${target}/output/${file_name}`)}>{file_name}</Link>
        </React.Fragment>
    )}
  </React.Fragment>;
}
