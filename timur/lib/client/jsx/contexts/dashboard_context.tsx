import React, {useState, useContext, useEffect, createContext, useCallback} from 'react';
import { Models, ModelsObject } from 'etna-js/models/magma-model';
import { CountsObject, MagmaContext } from 'etna-js/contexts/magma-context';
import { getModels, getAnswer, getIdentifiers, magmaPath } from 'etna-js/api/magma_api';
import {json_get, json_delete} from 'etna-js/utils/fetch';

type Workspace = {
  updated_at: string,
};

type Loader = {
  pipeline_finished_at: string,
}

const defaultDashboardState = {
  models: null as ModelsObject | null,
  rules: null as any | null,
  bytes: null as number | null,
  files: null as number | null,
  projectInfo: null as any | null,
  buckets: null as any | null,
  counts: null as CountsObject | null,
  loaders: null as Loader[] | null,
  workflows: null as any | null,
  workspaces: null as Workspace[] | null,
};

export type DashboardState = Readonly<typeof defaultDashboardState>;

export const defaultDashboardContext = {
  ...defaultDashboardState as DashboardState,
};

export type DashboardContextData = typeof defaultDashboardContext;

export const DashboardContext = createContext(defaultDashboardContext);
export type DashboardContext = typeof DashboardContext;
export type ProviderProps = {params?: {}; children: any};

export const DashboardProvider = (
  props: ProviderProps & Partial<DashboardContextData>
) => {
  const { models, setModels, mergeModels, counts, setCounts } = useContext(MagmaContext);
  useEffect( () => {
    const getMagmaData = async () => {
      const { models: m1 } = await getModels(CONFIG.project_name);
      const { models: m2 } = await getIdentifiers(CONFIG.project_name);
      const mergedModels = mergeModels(m1,m2);
      setModels(mergedModels);
    }

    if (Object.keys(models).length != 0) return;

    getMagmaData();
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

  useEffect( () => {
    const getRecords = async () => {
      let count = 0;
      const modelNames = Object.keys(models).filter( n => n != 'project' );

      let newCounts: CountsObject = {};

      for (const modelName of modelNames) {
        const { answer } = await getAnswer({query: [ modelName, '::count' ], project_name: CONFIG.project_name});
        newCounts[ modelName ] = answer;
      }
      setCounts({ ...counts, ...newCounts });
    }

    if (Object.keys(models).length == 0) return;

    getRecords();
  }, [models] );

  const [ loaders, setLoaders ] = useState<Loader[] | null>(null);
  useEffect( () => {
    json_get(`${CONFIG.polyphemus_host}/api/workflows/${CONFIG.project_name}/status`).then(
      loaders => setLoaders(loaders)
    );
  }, [] );

  const [ workflows, setWorkflows ] = useState(null);
  useEffect( () => {
    json_get(`${CONFIG.vulcan_host}/api/v2/${CONFIG.project_name}/workflows`).then(
      ({workflows}) => setWorkflows(workflows.length)
    );
  }, [] );

  const [ workspaces, setWorkspaces ] = useState(null);
  useEffect( () => {
    json_get(`${CONFIG.vulcan_host}/api/v2/${CONFIG.project_name}/workspace`).then(
      ({workspaces}) => {
        setWorkspaces(workspaces)
      }
    );
  }, [] );

  return (
    <DashboardContext.Provider
      value={ {
        models, rules, bytes, files, projectInfo, buckets, counts, loaders, workflows, workspaces
      }}
    >
      {props.children}
    </DashboardContext.Provider>
  );
};
