'use client'

import * as React from 'react'

const defaultProjectExplorerState = {
  visibleColumns: [ 'Name', 'Project ID', 'Data types', 'Theme' ],
  projectData: []
}

export const ProjectExplorerContext = React.createContext(defaultProjectExplorerState)

export function ProjectExplorerContextProvider({projectData, children}:{
    children: React.ReactNode,
}) {
  const [ state, setState ] = React.useState({
    ...defaultProjectExplorerState,
    projectData
  });

  const toggleColumnVisibility = React.useCallback(
    columnName => setState(
      {
        ...state,
        visibleColumns: state.visibleColumns.includes(columnName)
          ? state.visibleColumns.filter(c => columnName != c)
          : [ ...state.visibleColumns, columnName ]
      }
    )
  );

  return (
    <ProjectExplorerContext.Provider value={{
      state,
      toggleColumnVisibility
    }}>
      {children}
    </ProjectExplorerContext.Provider>
  )
}
