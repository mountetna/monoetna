'use client'

import * as React from 'react'
import { projectDataTypes } from '@/lib/utils/filters';
import { FilterItem, Project } from './models';

interface FilterSet {
  [filterName: string]: any
}

const defaultProjectExplorerState = {
  visibleColumns: [ 'Name', 'Project ID', 'Theme', 'Investigators' ],
  projectData: [] as Project[],
  filters: { } as FilterSet,
  filterItemSet: {} as FilterSet,
  matchAllFilters: false
}

interface ProjectExplorerContextValues {
  state: typeof defaultProjectExplorerState;
  filteredProjectData: Project[];
  updateFilterItems: (title: string, filterItems: any[] | null) => void;
  searchOptions: FilterItem[];
  setMatchAllFilters: (m:boolean) => void;
  createFilter: (filterName: string, filter: Function, items: Function | null) => void;
  toggleColumnVisibility: (columnName:string) => void;
}

export const ProjectExplorerContext = React.createContext<ProjectExplorerContextValues>({} as ProjectExplorerContextValues)

export function ProjectExplorerContextProvider({projectData, children}:{
    projectData: Project[],
    children: React.ReactNode,
}) {
  const [ state, setState ] = React.useState({
    ...defaultProjectExplorerState,
    projectData
  });

  const { filters } = state;

  const setMatchAllFilters = React.useCallback(
    (matchAllFilters:boolean) => setState({
      ...state,
      matchAllFilters
    }), [ state ]
  );

  const searchOptions = React.useMemo(() => {
    let options: FilterItem[] = [];

    Object.keys(filters).forEach(
      type => {
        if (!filters[type].items) return;
        let items: FilterSet = {};
        projectData.forEach((project:Project) => {
          items = { ...items, ...filters[type].items(project) }
        });
        Object.keys(items).forEach(
          itemName => options.push({
            type,
            value: items[itemName],
            label: typeof(items[itemName]) == 'object' ? items[itemName].name : items[itemName],
            key: type + '.' + itemName
          })
        );
      }
    );
    return options;
  }, [ filters ] );

  const freeFilter = React.useCallback(
    (filterItem:FilterItem['value'], project:Project, matchAllFilters:boolean) => {
      const labels: string[] = Object.values(filters).map(
        filter => filter.items ? Object.keys(filter.items(project)) : []
      ).flat().concat([
        project.name,
        project.fullName
      ]);
      return labels.some( label => label.toLowerCase().includes((filterItem as string).toLowerCase()) )
    }, [ filters ]
  );

  const createFilter = React.useCallback(
    (title:string, filter: Function, items: Function|null) => {
      setState( prevState => ({
        ...prevState,
        filters: {
          ...prevState.filters,
          [title]: { filter, items }
        }
      }))
    }, []
  )

  const updateFilterItems = React.useCallback(
    (filterName: string, filterItems: FilterItem['value'][]|null) => {
      let newState = { ...state };
      newState.filterItemSet = {
        ...state.filterItemSet,
        [ filterName ]: filterItems
      }

      if (!filterItems || filterItems.length == 0) {
        delete newState.filterItemSet[filterName];
      }
      setState(newState);
    }, [ state ]
  );

  const filteredProjectData = projectData.filter(
    (project:Project) => !Object.keys(state.filterItemSet).length ? true :
      Object.keys(state.filterItemSet)[
        state.matchAllFilters ? 'every' : 'some'
      ](
        filterName => state.filterItemSet[filterName][
          state.matchAllFilters ? 'every' : 'some'
        ](
          (filterItem:FilterItem['value']) => (
            filterName == 'free' ? freeFilter : state.filters[filterName].filter
          )(filterItem, project, state.matchAllFilters)
        )
      )
  )

  const toggleColumnVisibility = React.useCallback(
    (columnName:string) => setState( prevState => 
      ({
        ...prevState,
        visibleColumns: prevState.visibleColumns.includes(columnName)
          ? prevState.visibleColumns.filter(c => columnName != c)
          : [ ...prevState.visibleColumns, columnName ]
      })
    ), [state]
  );

  return (
    <ProjectExplorerContext.Provider value={{
      state,
      searchOptions,
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
