'use client'

import * as React from 'react'
import { projectDataTypes } from '@/lib/utils/filters';

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

  const createFilter = React.useCallback(
    (title, filter, collect, render) => {
      console.log("Creating "+title);
      setState({
        ...state,
        filters: {
          ...state.filters,
          [title]: { filter, collect, render }
        }
      })
    }, [state]
  )

  const updateFilterItems = React.useCallback(
    (filterName, filterItems) => {
      let newState = { ...state };
      newState.filterItemSet = {
        ...state.filterItemSet,
        [ filterName ]: filterItems
      }

      if (!filterItems || filterItems.length == 0) {
        delete newState.filterItemSet[filterName];
      }
      console.log({oldState: state, newState});
      setState(newState);
    }, [ state ]
  );

  const filteredProjectData = projectData.filter(
    project => !Object.keys(state.filterItemSet).length ? true :
      Object.keys(state.filterItemSet)[
        state.matchAllFilters ? 'every' : 'some'
      ](
        filterName => state.filterItemSet[filterName][
          state.matchAllFilters ? 'every' : 'some'
        ](
          filterItem => state.filters[filterName].filter(filterItem, project, state.matchAllFilters)
        )
      )
  )

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
      createFilter,
      updateFilterItems
    }}>
      {children}
    </ProjectExplorerContext.Provider>
  )
}
