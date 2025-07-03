import React, {useState, useContext, useEffect, createContext, useCallback} from 'react';
import { Models } from 'etna-js/models/magma-model';
import { MagmaContext } from 'etna-js/contexts/magma-context';
import {getModels} from 'etna-js/api/magma_api';

const defaultDashboardState = {
  models: null as Models | null,
};

export type DashboardState = Readonly<typeof defaultDashboardState>;

export const defaultDashboardContext = {
  state: defaultDashboardState as DashboardState,
};

export type DashboardContextData = typeof defaultDashboardContext;

export const DashboardContext = createContext(defaultDashboardContext);
export type DashboardContext = typeof DashboardContext;
export type ProviderProps = {params?: {}; children: any};

export const DashboardProvider = (
  props: ProviderProps & Partial<DashboardContextData>
) => {
  const [state, setState] = useState(props.state || defaultDashboardContext.state);

  const { models, setModels } = useContext(MagmaContext);

  useEffect( () => {
    if (Object.keys(models).length == 0) {
      getModels(project_name).then(({models}) => setModels(models));
    } else {
      setState({
        ...state,
        models: new Models(models)
      });
    }
  }, [models] );

  useEffect( () => {
    if (!('projectInfo' in state)) {
      json_get(
        `${CONFIG.janus_host}/api/admin/${CONFIG.project_name}/info`
      ).then(
        ({project}) => setState({
          ...state,
          projectInfo: project
        })
      )
    }
  }, [] );


  return (
    <DashboardContext.Provider
      value={ state }
    >
      {props.children}
    </DashboardContext.Provider>
  );
};
