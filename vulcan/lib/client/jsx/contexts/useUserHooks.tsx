import React, {useMemo, useState, useCallback} from 'react';
import _ from 'lodash';

import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {VulcanStorage, WorkspaceMinimal} from '../api_types';
import {isGuest} from 'etna-js/utils/janus';

export default function useUserHooks() {
  const user = useReduxState((state: any) => selectUser(state));

  const canEdit = useCallback(
    (workspace: WorkspaceMinimal | VulcanStorage) => {
      if (Object.keys(workspace).includes("workspace")) workspace = (workspace as VulcanStorage).workspace;
      // Todo before release: Remove null check once author implemented in API!
      return user.email === (workspace as WorkspaceMinimal).user_email
    }, [user]
  );

  const guest = useCallback(
    (projectName: string) => isGuest(user, projectName),
    [user]
  );

  return {
    canEdit,
    guest
  };
}
