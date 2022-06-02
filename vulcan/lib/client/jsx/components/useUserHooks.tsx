import React, {useMemo, useState, useCallback} from 'react';
import _ from 'lodash';

import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {VulcanFigure, VulcanFigureSession} from '../api_types';
import {isGuest} from 'etna-js/utils/janus';

export default function useUserHooks() {
  const user = useReduxState((state: any) => selectUser(state));

  console.log('user');
  const canEdit = useCallback(
    (figureSession: VulcanFigureSession | VulcanFigure) =>
      user.name === figureSession.author,
    [user]
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
