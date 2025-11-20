'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Grid from '@mui/system/Grid'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import Fade from '@mui/material/Fade';
import _ from 'lodash'
import { usePathname, useRouter, useSearchParams } from 'next/navigation';

import { ExternalProjectStatus, getExternalProjectStatus, PrincipalInvestigator, Project, ProjectHeadingInfoSet, PROJECTS_SEARCH_PARAMS_KEY, ProjectsSearchParamsControls, ProjectsSearchParamsState } from './models';
import ProjectListing from './project-listing';
import ProjectTable from './project-table';
import { ProjectExplorerContext } from './context';
import DrawerButton from '@/components/searchable-list/controls/drawer/button';
import Drawer from '@/components/searchable-list/controls/drawer/drawer';
import Autocomplete from '@/components/searchable-list/controls/autocomplete';
import { ThemeData } from '../themes/models';
import ProjectPI from './project-pi';
import FilterPill from '../searchable-list/filter-pill';
import FilterList from './filter-list';
import { ValueOf } from '@/lib/utils/types';
import { FILE_EXPORT_STATUS, handleExportFile, MIME_FILE_FORMATS } from '@/lib/utils/file-export';
import { parseSearchParams, toSearchParamsString } from '@/lib/utils/uri';
import { scrollTo } from '@/lib/utils/scroll';
import { useBreakpoint } from '@/lib/utils/responsive';
import { defaultDict, flattenObject } from '@/lib/utils/object';
import { FilterMethod } from '../searchable-list/models';

import filterLightIcon from '/public/images/icons/filter-light.svg'
import filterDarkIcon from '/public/images/icons/filter-dark.svg'
import searchDarkIcon from '/public/images/icons/search.svg'


interface FilterItem {
    // TODO: define value and type as separate interfaces? using Omit is a bit confusing
    value: (ValueOf<Omit<SearchableProjectData, 'principalInvestigators' | 'dataTypes' | 'theme'>>) | PrincipalInvestigator | ThemeData
    type: (keyof Omit<SearchableProjectData, 'principalInvestigators' | 'dataTypes' | 'fullName'>) | 'principalInvestigator' | 'dataType' | 'name' | 'hasClinicalData'
    label: string;
    key: string;
    projectKey: keyof SearchableProjectData;
}

const drawerFilterItemTypes: FilterItem['type'][] = ['theme', 'status', 'dataType']

const viewSets: FilterItem[] = Object.entries(ProjectHeadingInfoSet).map(([key, label]) => ({
    label,
    key,
    // Placeholder data just to satifsy the shape
    value: key,
    projectKey: 'fullName',
    type: 'name',
}))

const filterMethods: FilterItem[] = Object.entries(FilterMethod).map(([key, label]) => ({
    label,
    key,
    // Placeholder data just to satifsy the shape
    value: key,
    projectKey: 'fullName',
    type: 'name',
}))

const PIExportAttrs: (keyof PrincipalInvestigator)[] = ['name', 'title']
const ThemeExportAttrs: (keyof ThemeData)[] = ['name', 'description', 'projectsLink']


// TODO: use separate list component in searchable-list module
function _ProjectExplorer({ }) {
    // Manage search params sync
    const { state: { projectData } } = React.useContext(ProjectExplorerContext);
    const router = useRouter()
    const pathname = usePathname()
    const searchParams = useSearchParams()
    const filterItems = [];

    React.useEffect(() => {
        const parsedSearchParams = parseSearchParams(searchParams)

        if (PROJECTS_SEARCH_PARAMS_KEY in parsedSearchParams) {
            const state = parsedSearchParams[PROJECTS_SEARCH_PARAMS_KEY]
            let updatedState = false

            const stateFilterItems = parseSearchOptionsFromState(state, searchOptions)

            if (
                stateFilterItems !== undefined &&
                // Prevent item shuffling if params match state but different order
                !_.isEqual(
                    [...stateFilterItems].sort((a, b) => a.key.localeCompare(b.key)),
                    [...filterItems].sort((a, b) => a.key.localeCompare(b.key)),
                )
            ) {
                setFilterItems(stateFilterItems)
                updatedState = true
            }

            const stateViewSet = parseViewSetFromState(state, viewSets)
            if (stateViewSet !== undefined && viewSet.key !== stateViewSet.key) {
                setViewSet(stateViewSet)
                updatedState = true
            }

            const stateFilterMethod = parseFilterMethodFromState(state, filterMethods)
            if (stateFilterMethod !== undefined && viewSet.key !== stateFilterMethod.key) {
                setFilterMethod(stateFilterMethod)
                updatedState = true
            }

            const stateCurrentPage = parseCurrentPageFromState(state)
            if (currentPage !== stateCurrentPage) {
                setCurrentPage(stateCurrentPage)
                updatedState = true
            }

            if (updatedState) {
                closeAllProjects()
            }
        }
    }, [searchParams])

    const theme = useTheme()
    const breakpoint = useBreakpoint()
    const isMobile = breakpoint === 'mobile'

    const countsByFilterItemType: Record<string, Record<string, number>> = defaultDict(_ => defaultDict(_ => 0))

    const searchPlaceholder = 'Search e.g. "Fibrosis"'

    const [drawerOpen, setDrawerOpen] = React.useState(true)

    const [viewSet, setViewSet] = React.useState<(typeof viewSets)[number]>(viewSets[0])

    // Sync url with filter+control states
    const updateUrl = (filterItems: FilterItem[], viewSet: FilterItem, filterMethod: FilterItem, currentPage: number) => {
        // get filter states
        const filters: Record<string, string[]> = {}
        filterItems.forEach(item => {
            let val: string
            switch (item.type) {
                case 'theme':
                case 'type':
                case 'status':
                case 'dataType':
                case 'name':
                case 'hasClinicalData':
                    val = (item.value as string)
                    break
                case 'principalInvestigator':
                    val = (item.value as PrincipalInvestigator).name
                    break
            }
            if (!(item.type in filters)) {
                filters[item.type] = []
            }
            filters[item.type].push(val)
        })

        // get control states
        const controls: ProjectsSearchParamsControls = {
            viewSet: viewSet.key,
            filterMethod: filterMethod.key,
            page: currentPage,
        }

        const projectsState: ProjectsSearchParamsState = {
            filters,
            controls,
        }

        // push to router
        router.push(pathname + '?' + toSearchParamsString({ [PROJECTS_SEARCH_PARAMS_KEY]: projectsState }) + window.location.hash, { scroll: false })
    }


    const handleSetCurrentPage = (page: number) => {
        closeAllProjects()
        setCurrentPage(page)
        updateUrl(filterItems, viewSet, filterMethod, page)
    }

    const handleChangeFilterItems = (filterItems: FilterItem[]) => {
        setFilterItems(filterItems)
        closeAllProjects()
        setCurrentPage(0)
        updateUrl(filterItems, viewSet, filterMethod, 0)
    }

    const handleClickRemoveFilterItem = (filterItem: FilterItem) => {
        const newFilterItems = filterItems.filter((item) => item.key !== filterItem.key)
        setFilterItems(newFilterItems)
        closeAllProjects()
        setCurrentPage(0)
        updateUrl(newFilterItems, viewSet, filterMethod, 0)
    }

    const handleChangeViewSet = (viewSet: FilterItem) => {
        setViewSet(viewSet)
        updateUrl(filterItems, viewSet, filterMethod, currentPage)
    }

    const handleChangeFilterMethod = (filterMethod: FilterItem) => {
        setFilterMethod(filterMethod)
        updateUrl(filterItems, viewSet, filterMethod, currentPage)
    }

    // Manage file export
    const [fileExportStatus, setFileExportStatus] = React.useState(FILE_EXPORT_STATUS.idle)

    const handleClickExportButton = async () => {
        await handleExportFile(
            filteredProjectData.map(proj => {
                const _proj = _.cloneDeep(proj)
                _proj.theme = _.pick(_proj.theme, ThemeExportAttrs)
                _proj.principalInvestigators = _proj.principalInvestigators.map(pi => _.pick(pi, PIExportAttrs))

                return flattenObject(_proj, '.', '', true)
            }),
            MIME_FILE_FORMATS.csv,
            setFileExportStatus,
            'ucsf-data-library-projects',
        )
    }

    return (
        <Container
            sx={{
                display: 'flex',
                flexDirection: 'column',
                gap: '24px',
                [theme.breakpoints.up('tablet')]: {
                    gap: '38px',
                },
                mb: '50px'
            }}
        >
            <Typography
                variant='h2'
            >
                Project explorer
            </Typography>

            <Box
                sx={{
                    display: 'flex',
                    flexDirection: 'column',
                    gap: '16px',
                    [theme.breakpoints.up('tablet')]: {
                        gap: '8px',
                    },
                }}
            >
                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
                        [theme.breakpoints.up('tablet')]: {
                            display: 'flex',
                            flexDirection: 'row',
                            justifyContent: 'space-between',
                            width: '100%',
                        },
                        paddingBottom: '16px',
                        borderBottom: `1px solid ${theme.palette.ground.grade50}`,
                        borderEndEndRadius: '2px'
                    }}
                >
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            columnGap: '8px',
                            pb: '16px',
                            [theme.breakpoints.up('tablet')]: {
                                width: '100%',
                                overflow: 'visible',
                                pb: '0',
                                rowGap: '0',
                            },
                        }}
                    >
                        <Box
                            sx={{
                                display: 'flex',
                                flexDirection: 'row',
                                columnGap: '8px',
                                height: '52px',
                                [theme.breakpoints.up('tablet')]: {
                                    rowGap: '0',
                                    columnGap: '10px',
                                    justifyContent: 'space-between',
                                    width: '100%',
                                },
                            }}
                        >
                            <Box
                                sx={{
                                    display: 'flex',
                                    flexDirection: 'row',
                                    columnGap: '8px',
                                    flexGrow: 1,
                                    width: '100%',
                                    [theme.breakpoints.up('tablet')]: {
                                        flexGrow: 'unset',
                                        columnGap: '10px',
                                    },
                                }}
                            >
                                <DrawerButton
                                    label='Filters'
                                    iconLight={filterLightIcon}
                                    iconDark={filterDarkIcon}
                                    onClick={() => setDrawerOpen(!drawerOpen)}
                                    activated={filterItems.find(item => drawerFilterItemTypes.includes(item.type)) !== undefined}
                                    open={drawerOpen}
                                />

                                <Autocomplete
                                    multiple
                                    filterSelectedOptions
                                    icon={searchDarkIcon}
                                    placeholder={searchPlaceholder}
                                    options={[]}
                                    showResultsOnEmpty={false}
                                    // @ts-ignore
                                    onChange={(_, value: FilterItem[], reason) => {
                                        // disables removing options with backspace
                                        if (reason === 'removeOption') return

                                        handleChangeFilterItems(value)
                                    }}
                                    value={filterItems}
                                    // @ts-ignore
                                    getOptionKey={(option: FilterItem) => option.key}
                                    isOptionEqualToValue={(option, value) => option.key === value.key}
                                    renderOption={(params) => (
                                        <SearchOption
                                            key={params.option.key}
                                            option={params.option}
                                        />
                                    )}
                                    groupBy={option => option.type}
                                    renderGroup={(params) => (
                                        <Box
                                            key={params.key}
                                            role='listitem'
                                            sx={{
                                                '&:not(:last-child)': {
                                                    mb: '10px',
                                                },
                                            }}
                                        >
                                            <Typography
                                                variant='pMediumMediumWt'
                                                component='div'
                                                sx={{
                                                    color: 'ground.grade10',
                                                    opacity: 0.4,
                                                    p: '4.11px',
                                                    pl: '0px',
                                                    mb: '8px',
                                                }}
                                            >
                                                {_.startCase(params.group)}
                                            </Typography>
                                            <Box
                                                role='list'
                                                sx={{
                                                    display: 'flex',
                                                    flexDirection: 'column',
                                                    gap: '8px',
                                                }}
                                            >
                                                {params.options}
                                            </Box>
                                        </Box>
                                    )}
                                    renderNoResults={() => (
                                        <Typography variant='pMedium'>
                                            No results
                                        </Typography>
                                    )}
                                    sx={{
                                        flexGrow: 1,
                                        maxHeight: '100%',
                                        [theme.breakpoints.up('tablet')]: {
                                            flexGrow: 'unset',
                                        },
                                        '& .MuiBox-root': {
                                            width: '100%',
                                            [theme.breakpoints.up('tablet')]: {
                                                width: 'auto',
                                                '& input': {
                                                    flexGrow: 'unset',
                                                    width: `${searchPlaceholder.length * 0.45}em`,
                                                },
                                            },
                                        },
                                    }}
                                />

                            </Box>
                        </Box>

                          {filterItems.length > 0 && <Box
                              sx={{
                                  display: 'flex',
                                  flex: '1 1 auto',
                                  flexDirection: 'row',
                                  flexWrap: 'wrap',
                                  columnGap: '10px',
                                  rowGap: '16px',
                                  pt: '16px',
                                  [theme.breakpoints.up('tablet')]: {
                                  },
                              }}
                          >
                              {filterItems.map((item => (
                                  <FilterPill
                                      key={item.key}
                                      label={item.label}
                                      removeable
                                      onClickRemove={() => handleClickRemoveFilterItem(item)}
                                  />
                              )))}
                          </Box>}
                    </Box>
                </Box>

                <Box sx={{
                  display: 'flex'
                }}>
                  {drawerOpen && <Box
                    sx={ theme => ({
                      flex: '0 0 405px',
                      pt: '9px',
                      pr: '5px',
                      borderRight: `1px solid ${theme.palette.ground.grade50}`,
                      transition: theme.transitions.create(
                        ['flex'],
                        {
                           easing: theme.transitions.easing.ease,
                           duration: theme.transitions.duration.ease,
                        }
                      ),
                      borderEndEndRadius: '2px'
                    })}>
                     <FilterList/>
                    </Box>
                  }
                  <Box
                      role='list'
                      sx={{ flex: '1 1 auto' }}
                  >
                    <Box>
                      <ProjectTable/>
                    </Box>
                  </Box>
              </Box>
            </Box>
        </Container >
    )
}

// Needed for https://nextjs.org/docs/messages/missing-suspense-with-csr-bailout
export default function ProjectExplorer({ }) {
    return (
        <React.Suspense fallback={null}>
            <_ProjectExplorer/>
        </React.Suspense>
    )
}


function SearchOption({
    option,
}: {
    option: FilterItem,
}) {
    return (
        <Box
            sx={{
                p: '8px',
            }}
        >
            {option.type === 'principalInvestigator' ?
                (
                    <ProjectPI
                        data={option.value as PrincipalInvestigator}
                        variant='filled'
                    />
                ) : (
                    <Typography variant='pMedium'>
                        {/* @ts-ignore */}
                        {option.value}
                    </Typography>
                )
            }

        </Box>
    )
}

function getFilterItemLabel(value: FilterItem['value'], type: FilterItem['type']): string {
    switch (type) {
        case 'name':
        case 'type':
        case 'status':
        case 'dataType':
        case 'theme':
        case 'hasClinicalData':
            // @ts-ignore
            return value
        case 'principalInvestigator':
            return (value as PrincipalInvestigator).name
    }
}

function getFilterItemKey(value: FilterItem['value'], type: FilterItem['type']): string {
    return `${type}_${getFilterItemLabel(value, type)}`
}

function parseSearchOptionsFromState(state: ProjectsSearchParamsState, searchOptions: FilterItem[]): FilterItem[] | undefined {
    if (!(state.filters)) {
        return undefined
    }
    return searchOptions.filter(searchOption => {
        if (state.filters && searchOption.type in state.filters) {
            // @ts-ignore
            const vals: string[] = state.filters[searchOption.type]

            for (const val of vals) {
                let optionVal: string
                switch (searchOption.type) {
                    case 'theme':
                    case 'type':
                    case 'status':
                    case 'dataType':
                    case 'name':
                    case 'hasClinicalData':
                        // @ts-ignore
                        optionVal = searchOption.value
                        break
                    case 'principalInvestigator':
                        optionVal = (searchOption.value as PrincipalInvestigator).name
                }
                if (val === optionVal) {
                    return true
                }
            }
            return false
        }
        return false
    })
}

function parseViewSetFromState(state: ProjectsSearchParamsState, viewSets: FilterItem[]): FilterItem | undefined {
    const viewSetValue = state.controls?.viewSet
    if (viewSetValue !== undefined) {
        return viewSets.find(item => item.key === viewSetValue)
    }
}

function parseFilterMethodFromState(state: ProjectsSearchParamsState, filterMethods: FilterItem[]): FilterItem | undefined {
    const filterMethodValue = state.controls?.filterMethod
    if (filterMethodValue !== undefined) {
        return filterMethods.find(item => item.key === filterMethodValue)
    }
}

function parseCurrentPageFromState(state: ProjectsSearchParamsState): number {
    const paramsPage = state.controls?.page
    return paramsPage !== undefined ? paramsPage : 0
}
