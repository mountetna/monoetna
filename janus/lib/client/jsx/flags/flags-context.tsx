import React, {useState, createContext} from 'react';

import {UserFlagsInterface, Project} from '../types/janus_types';

const defaultFlagsState = {
  users: [] as UserFlagsInterface[],
  projects: [] as Project[]
};

type FlagsState = Readonly<typeof defaultFlagsState>;

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
  const [state, setState] = useState({} as FlagsState);

  function setUsers(users: UserFlagsInterface[]) {
    setState({...state, users});
  }

  function setProjects(projects: Project[]) {
    setState({...state, projects});
  }

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
