'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import Fade from '@mui/material/Fade';
import _ from 'lodash'
import { usePathname, useRouter, useSearchParams } from 'next/navigation';

import { ExternalProjectStatus, getExternalProjectStatus, PrincipalInvestigator, Project, ProjectHeadingInfoSet, PROJECTS_SEARCH_PARAMS_KEY, ProjectsSearchParamsControls, ProjectsSearchParamsState } from './models';
import ProjectListing from './project-listing';
import Pagination, { PaginationClasses } from '@/components/searchable-list/controls/pagination'
import DrawerButton from '@/components/searchable-list/controls/drawer/button';
import Drawer from '@/components/searchable-list/controls/drawer/drawer';
import Autocomplete from '@/components/searchable-list/controls/autocomplete';
import { ThemeData } from '../themes/models';
import ProjectPI from './project-pi';
import FilterPill from '../searchable-list/filter-pill';
import { ValueOf } from '@/lib/utils/types';
import { FILE_EXPORT_STATUS, handleExportFile, MIME_FILE_FORMATS } from '@/lib/utils/file-export';
import { DrawerSectionClasses, DrawerClasses } from '../searchable-list/controls/drawer/models';
import { parseSearchParams, toSearchParamsString } from '@/lib/utils/uri';
import { scrollTo } from '@/lib/utils/scroll';
import { useBreakpoint } from '@/lib/utils/responsive';
import { defaultDict, flattenObject } from '@/lib/utils/object';
import { FilterMethod } from '../searchable-list/models';

import filterLightIcon from '/public/images/icons/filter-light.svg'
import filterDarkIcon from '/public/images/icons/filter-dark.svg'
import searchDarkIcon from '/public/images/icons/search.svg'


const searchableProjectKeys: (keyof SearchableProjectData)[] = ['fullName', 'principalInvestigators', 'type', 'dataTypes', 'theme', 'status', 'hasClinicalData']

interface SearchableProjectData extends Pick<Project, 'fullName' | 'principalInvestigators' | 'type' | 'dataTypes' | 'theme' | 'hasClinicalData'> {
    status: ExternalProjectStatus
}

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
function _ProjectListings({
    projectData,
}: {
    projectData: Project[],
}) {
    // Manage search params sync
    const router = useRouter()
    const pathname = usePathname()
    const searchParams = useSearchParams()

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

    const [currentPage, setCurrentPage] = React.useState(0)
    const pageSize = 12

    const countsByFilterItemType: Record<string, Record<string, number>> = defaultDict(_ => defaultDict(_ => 0))

    const [searchOptions] = React.useState(() => {

        const _searchOptions: Set<FilterItem> = new Set()

        const existingOptions = new Set<string>()

        function addOption(value: FilterItem['value'], type: FilterItem['type'], projectKey: keyof SearchableProjectData) {
            const key = getFilterItemKey(value, type)

            countsByFilterItemType[type][key] += 1

            if (!existingOptions.has(key)) {
                existingOptions.add(key)
                _searchOptions.add({
                    value,
                    type,
                    label: getFilterItemLabel(value, type),
                    key,
                    projectKey,
                })
            }
        }

        for (const project of projectData) {
            for (const projectKey of searchableProjectKeys) {

                let value = project[projectKey]
                let type: FilterItem['type']

                switch (projectKey) {
                    case 'type':
                    case 'hasClinicalData':
                        type = projectKey
                        break
                    case 'status':
                        value = getExternalProjectStatus(project)
                        type = projectKey
                        break
                    case 'theme':
                        value = (value as ThemeData).name
                        type = projectKey
                        break
                    case 'fullName':
                        type = 'name'
                        break
                    case 'principalInvestigators':
                        type = 'principalInvestigator'
                        break
                    case 'dataTypes':
                        type = 'dataType'
                        break
                }

                if (Array.isArray(value)) {
                    for (const item of value) {
                        addOption(
                            item,
                            type,
                            projectKey
                        )
                    }
                } else {
                    addOption(
                        value,
                        type,
                        projectKey,
                    )
                }
            }
        }
        return [..._searchOptions]
    })

    const searchPlaceholder = 'Search e.g. "Fibrosis"'

    const [drawerItems] = React.useState(() => {
        const filtered = searchOptions
            .filter((option) => {
                return drawerFilterItemTypes.includes(option.type)
            })

        const itemsByType: Record<string, FilterItem[]> = {}
        for (const item of filtered) {
            if (!(item.type in itemsByType)) {
                itemsByType[item.type] = []
            }
            itemsByType[item.type].push(item)
        }

        const sorted: FilterItem[] = []
        for (const [type, items] of Object.entries(itemsByType)) {
            const counts = countsByFilterItemType[type]
            items.sort((a, b) => counts[b.key] - counts[a.key])
            sorted.push(...items)
        }
        return sorted
    })

    const [drawerOpen, setDrawerOpen] = React.useState(false)

    const [filterItems, setFilterItems] = React.useState<FilterItem[]>([])

    const [filterMethod, setFilterMethod] = React.useState<(typeof filterMethods)[number]>(filterMethods[0])

    const filteredProjectData = projectData
        .filter((project) => {
            if (filterItems.length === 0) return true;

            function doesFilterItemMatch(filterItem: FilterItem): boolean {
                const projectKey = filterItem.projectKey
                let projectValue

                switch (projectKey) {
                    case 'fullName':
                    case 'type':
                    case 'hasClinicalData':
                        projectValue = project[projectKey]
                        break
                    case 'dataTypes':
                        for (const dt of project.dataTypes) {
                            if (dt === filterItem.value) {
                                return true
                            }
                        }
                        return false
                    case 'theme':
                        projectValue = project.theme.name
                        break
                    case 'status':
                        projectValue = getExternalProjectStatus(project)
                        break
                    case 'principalInvestigators':
                        for (const pi of project.principalInvestigators) {
                            if (pi.name === (filterItem.value as PrincipalInvestigator).name) {
                                return true
                            }
                        }
                        return false
                }

                return projectValue === filterItem.value
            }

            switch (filterMethod.label) {
                case FilterMethod.all:
                    return filterItems.every(item => doesFilterItemMatch(item))
                case FilterMethod.any:
                    return filterItems.some(item => doesFilterItemMatch(item))
            }
        })

    const paginatedFilteredProjectData = filteredProjectData
        .slice(currentPage * pageSize, currentPage * pageSize + pageSize)
        .map(proj => ({
            project: proj,
            nodeRef: React.createRef<HTMLElement>(),
        }))


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


    // Handle book open state to coordinate only having one open at a time
    const [projectOpens, setProjectOpens] = React.useState(projectData.map(_ => false))

    const handleSetProjectOpen = (newOpenState: boolean, projectIdx: number) => {
        setProjectOpens(opens => opens.map((_, idx) => {
            if (idx === projectIdx) {
                return newOpenState
            }
            return false
        }))

    }

    const closeAllProjects = () => {
        setProjectOpens(opens => opens.map(_ => false))
    }

    const handleSetCurrentPage = (page: number) => {
        closeAllProjects()
        setCurrentPage(page)
        updateUrl(filterItems, viewSet, filterMethod, page)
    }

    const scrollToProject = (nodeRef: React.RefObject<HTMLElement>) => {
        const projectEl = nodeRef.current
        if (projectEl === null) {
            return
        }

        scrollTo({ top: projectEl.offsetTop }, breakpoint)
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

    const paginationEl = (
        <Pagination
            currentPage={currentPage}
            pageSize={pageSize}
            listSize={filteredProjectData.length}
            onClickPrev={() => handleSetCurrentPage(currentPage - 1)}
            onClickNext={() => handleSetCurrentPage(currentPage + 1)}
            listItemLabel='projects'
            sx={{
                [theme.breakpoints.up('tablet')]: {
                    height: '100%',
                    [`& .${PaginationClasses.pageInfo}`]: {
                        height: '100%',
                    },
                },
            }}
        />
    )

    return (
        <Container
            sx={{
                display: 'flex',
                flexDirection: 'column',
                gap: '24px',
                [theme.breakpoints.up('tablet')]: {
                    gap: '38px',
                },
            }}
        >
            <Typography
                variant='h1'
            >
                Projects
            </Typography>

            <Box
                sx={{
                    display: 'flex',
                    flexDirection: 'column',
                    gap: '16px',
                    [theme.breakpoints.up('tablet')]: {
                        gap: '38px',
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
                        }
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
                                    options={[...searchOptions].sort((a, b) => a.type.localeCompare(b.type))}
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

                            {!isMobile && paginationEl}
                        </Box>

                        <Box
                            sx={{
                                [`& .${DrawerClasses.root}`]: {
                                    pt: '16px',
                                    [theme.breakpoints.up('tablet')]: {
                                        pt: '13px',
                                    },
                                    [theme.breakpoints.up('desktop')]: {
                                    },
                                },
                                [`& .${DrawerSectionClasses.base}`]: {
                                    [theme.breakpoints.up('tablet')]: {
                                        width: '33%',
                                    },
                                    [theme.breakpoints.up('desktop')]: {
                                        minWidth: 'unset',
                                        width: '33%',
                                    },
                                },
                            }}
                        >
                            <Drawer
                                items={drawerItems}
                                viewSetItems={viewSets}
                                filterMethodItems={filterMethods}
                                getItemLabel={item => item.label}
                                getItemKey={item => item.key}
                                getItemType={item => item.type}
                                activeItems={filterItems}
                                onChange={handleChangeFilterItems}
                                activeViewSetItem={viewSet}
                                onChangeViewSet={handleChangeViewSet}
                                showViewSets={true}
                                activeFilterMethodItem={filterMethod}
                                onChangeFilterMethod={handleChangeFilterMethod}
                                showFilterMethods={true}
                                renderFilterMethod={(toggleEl) => (
                                    <React.Fragment>
                                        <Typography variant='pBody'>
                                            Show projects that match
                                        </Typography>

                                        {toggleEl}

                                        <Typography variant='pBody'>
                                            filters
                                        </Typography>
                                    </React.Fragment>
                                )}
                                open={drawerOpen}
                                showButton={true}
                                buttonLabel={'Export to CSV'}
                                onClickButton={handleClickExportButton}
                                buttonDisabled={fileExportStatus === FILE_EXPORT_STATUS.inProgress}
                                displayStyle={isMobile ? 'expandable' : 'default'}
                            />
                        </Box>

                        {filterItems.length > 0 && <Box
                            sx={{
                                display: 'flex',
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

                    {isMobile && paginationEl}
                </Box>


                <Box
                    role='list'
                    sx={{
                        bgcolor: 'utilityWhite.main',
                        borderRadius: '20px',
                        overflow: 'hidden',
                    }}
                >
                    <Fade
                        key={paginatedFilteredProjectData.reduce((prev, curr) => prev + curr.project.fullName, '')}
                        in={true}
                        easing={theme.transitions.easing.ease}
                        timeout={theme.transitions.duration.ease}
                        exit={false}
                    >
                        <Box>
                            {paginatedFilteredProjectData.map(({ project, nodeRef }, i) => (
                                <Box
                                    key={project.fullName}
                                    ref={nodeRef}
                                    role='listitem'
                                >
                                    <ProjectListing
                                        data={project}
                                        open={projectOpens[i]}
                                        headingInfoSet={(viewSet.label as ProjectHeadingInfoSet)}
                                        onSetOpen={open => handleSetProjectOpen(open, i)}
                                        onFinishOpen={() => scrollToProject(nodeRef)}
                                    />
                                </Box>
                            ))}
                        </Box>
                    </Fade>
                </Box>
            </Box>
        </Container >
    )
}

// Needed for https://nextjs.org/docs/messages/missing-suspense-with-csr-bailout
export default function ProjectListings({
    projectData,
}: {
    projectData: Project[],
}) {
    return (
        <React.Suspense fallback={null}>
            <_ProjectListings
                projectData={projectData}
            />
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