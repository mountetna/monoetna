'use client'

import * as React from 'react'

const defaultProjectExplorerState = {
  visibleColumns: [ 'Name', 'Project ID', 'Theme', 'Investigators' ],
  projectData: [],
  filters: {},
  filterItemSet: {},
  matchAllFilters: false
}

export const ProjectExplorerContext = React.createContext(defaultProjectExplorerState)

export function ProjectExplorerContextProvider({projectData, children}:{
    children: React.ReactNode,
}) {
  const [ state, setState ] = React.useState({
    ...defaultProjectExplorerState,
    projectData
  });

  const setMatchAllFilters = React.useCallback(
    (matchAllFilters) => setState({
      ...state,
      matchAllFilters
    }), [ state ]
  );

  const updateFilter = React.useCallback(
    (filterName, filter, filterItems) => {
      let newState = {
        ...state,
        filters: {
          ...state.filters,
          [filterName]: filter
        }
      };

      if (filterItems) {
        newState.filterItemSet = {
          ...state.filterItemSet,
          [ filterName ]: filterItems
        }
      }

      if (!filter) {
        delete newState.filters[filterName];
        delete newState.filterItemSet[filterName];
      }

      setState(newState);
    }, [ state ]
  );

  const filteredProjectData = projectData.filter(
    project => !Object.keys(state.filters).length ? true :
      Object.values(state.filters)[
        state.matchAllFilters ? 'every' : 'some'
      ](
        filter => filter(project, state.matchAllFilters)
      )
  )

  console.log({state});

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
      filteredProjectData,
      toggleColumnVisibility,
      setMatchAllFilters,
      updateFilter
    }}>
      {children}
    </ProjectExplorerContext.Provider>
  )
}
