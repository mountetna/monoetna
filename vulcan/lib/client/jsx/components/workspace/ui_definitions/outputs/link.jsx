import React, {useContext} from 'react';
import {VulcanContext} from '../../../../contexts/vulcan_context';
import Link from 'etna-js/components/link';

function _LinkOutput({data, stream}) {
  let {state, vulcanPath} = useContext(VulcanContext);
  let {workspace, projectName} = state;
  if (!workspace) {
    return
  }
  const target = `${projectName}/workspace/${workspace.workspace_id}/file/${stream ? 'stream-download' : 'download'}`
  return <React.Fragment>
    {Object.keys(data).map(file_name =>
        <React.Fragment key={file_name}>
          <Link link={vulcanPath(`/api/v2/${target}/${file_name}`)}>{file_name}</Link>
        </React.Fragment>
    )}
  </React.Fragment>;
}

export function LinkOutput({data}) {
  return _LinkOutput({data, stream: false})
}

export function LinkLargeOutput({data}) {
  return _LinkOutput({data, stream: true})
}
