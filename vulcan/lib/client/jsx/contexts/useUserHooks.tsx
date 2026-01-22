import React, {useMemo, useState, useCallback} from 'react';
import _ from 'lodash';

import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {VulcanStorage, Workspace, WorkspaceMinimal} from '../api_types';
import {isGuest, isSuperuser} from 'etna-js/utils/janus';

export default function useUserHooks() {
  const user = useReduxState((state: any) => selectUser(state));
  
  const canEdit = useCallback(
    (workspaceOrStorage: Workspace | WorkspaceMinimal | VulcanStorage) => {
      const workspace = ('workspace' in workspaceOrStorage) ? workspaceOrStorage.workspace : workspaceOrStorage
      return user.email === workspace.user_email
    }, [user]
  );

  const guest = useCallback(
    (projectName: string) => isGuest(user, projectName) as boolean,
    [user]
  );

  const superuser = useMemo(
    () => isSuperuser(user) as boolean,
    [user]
  )

  return {
    canEdit,
    guest,
    superuser,
  };
}
