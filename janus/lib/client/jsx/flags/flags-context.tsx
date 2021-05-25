import React, {useState, useEffect, createContext} from 'react';

import {UserFlagsInterface, Project} from '../types/janus_types';

const defaultFlagsState = {
  users: [] as UserFlagsInterface[],
  projects: [] as Project[]
};

export type FlagsState = Readonly<typeof defaultFlagsState>;

export const defaultContext = {
  state: defaultFlagsState as FlagsState,
  setUsers: (users: UserFlagsInterface[]) => {},
  setProjects: (projects: Project[]) => {}
};

export type FlagsContextData = typeof defaultContext;
export const FlagsContext = createContext(defaultContext);
export type FlagsContext = typeof FlagsContext;
export type ProviderProps = {
  params?: {};
  children: any;
};

export const FlagsProvider = (
  props: ProviderProps & Partial<FlagsContextData>
) => {
  const [state, setState] = useState(props.state || defaultFlagsState);
  const [projects, setProjects] = useState([] as Project[]);
  const [users, setUsers] = useState([] as UserFlagsInterface[]);

  useEffect(() => {
    setState({...state, users, projects});
  }, [users, projects]);

  return (
    <FlagsContext.Provider
      value={{
        state,
        setUsers,
        setProjects
      }}
    >
      {props.children}
    </FlagsContext.Provider>
  );
};
