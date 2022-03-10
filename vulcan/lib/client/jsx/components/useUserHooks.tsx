import React, {useMemo, useState, useCallback} from 'react';
import _ from 'lodash';

import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {VulcanFigure, VulcanFigureSession} from '../api_types';

export default function useUserHooks() {
  const user = useReduxState((state: any) => selectUser(state));

  const canEdit = useCallback(
    (figureSession: VulcanFigureSession | VulcanFigure) =>
      user.name === figureSession.author,
    [user]
  );

  return {
    canEdit
  };
}
