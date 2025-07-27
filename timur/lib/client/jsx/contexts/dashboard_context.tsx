import React, {useState, useContext, useEffect, createContext, useCallback} from 'react';
import { Models } from 'etna-js/models/magma-model';
import { MagmaContext } from 'etna-js/contexts/magma-context';
import { getModels, magmaPath } from 'etna-js/api/magma_api';
import {json_get, json_delete} from 'etna-js/utils/fetch';

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
  const { models, setModels } = useContext(MagmaContext);
  useEffect( () => {
    if (Object.keys(models).length == 0) {
      getModels(CONFIG.project_name).then(({models}) => setModels(models));
    }
  }, [models] );

  const [ rules, setRules ] = useState(null);
  useEffect( () => {
    json_get(magmaPath(`gnomon/${CONFIG.project_name}/`)).then(
      ({config: { rules } }) => setRules(rules)
    );
  }, [] );

  const [ projectInfo, setProjectInfo ] = useState(null);
  useEffect( () => {
    json_get(`${CONFIG.janus_host}/api/admin/${CONFIG.project_name}/info`).then(
      ({project}) => setProjectInfo(project)
    );
  }, [] );

  const [ buckets, setBuckets ] = useState(null);
  useEffect( () => {
    json_get(`${CONFIG.metis_host}/${CONFIG.project_name}/list`).then(
      ({buckets}) => setBuckets(buckets)
    );
  }, [] );

  const [ files, setFiles ] = useState(null);
  useEffect( () => {
    json_get(`${CONFIG.metis_host}/api/stats/files/${CONFIG.project_name}`).then(
      ({[CONFIG.project_name]: files}) => setFiles(files)
    );
  }, [] );

  const [ bytes, setBytes ] = useState(null);
  useEffect( () => {
    json_get(`${CONFIG.metis_host}/api/stats/bytes/${CONFIG.project_name}`).then(
      ({[CONFIG.project_name]: bytes}) => setBytes(bytes)
    );
  }, [] );

  const [ records, setRecords ] = useState(null);
  useEffect( () => {
  }, [] );

  return (
    <DashboardContext.Provider
      value={ {
        models, rules, bytes, files, projectInfo, buckets
      }}
    >
      {props.children}
    </DashboardContext.Provider>
  );
};
