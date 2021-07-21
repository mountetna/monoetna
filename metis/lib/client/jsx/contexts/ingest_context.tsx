import React, {useState, createContext, useCallback} from 'react';

export type Host = {
  alias: string;
  host: string;
  directories: string[];
};

const defaultState = {
  hosts: {} as {[key: string]: Host}
};

export const defaultContext = {
  state: defaultState as IngestState,
  setHosts: (hosts: {[key: string]: Host}) => {}
};

export type IngestState = Readonly<typeof defaultState>;
export type IngestContextData = typeof defaultContext;

export const IngestContext = createContext(defaultContext);
export type IngestContext = typeof IngestContext;
export type ProviderProps = {params?: {}; children: any};

export const IngestProvider = (
  props: ProviderProps & Partial<IngestContextData>
) => {
  const [state, setState] = useState(props.state || defaultContext.state);

  const setHosts = useCallback(
    (hosts: {[key: string]: Host}) => {
      setState({
        ...state,
        hosts: {...hosts}
      });
    },
    [state]
  );

  return (
    <IngestContext.Provider
      value={{
        state,
        setHosts
      }}
    >
      {props.children}
    </IngestContext.Provider>
  );
};
