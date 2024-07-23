'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useMediaQuery, useTheme } from '@mui/material';
import _ from 'lodash'

import { DataType, ExternalProjectStatus, getExternalProjectStatus, PrincipalInvestigator, Project, ProjectHeadingInfoSet } from './models';
import ProjectListing from './project-listing';
import Pagination from '@/components/searchable-list/controls/pagination'
import DrawerButton from '@/components/searchable-list/controls/drawer/button';
import Drawer from '@/components/searchable-list/controls/drawer/drawer';
import Autocomplete from '@/components/searchable-list/controls/autocomplete';
import ToggleGroup from '../inputs/toggle-group';
import { ThemeData } from '../themes/models';
import filterDarkIcon from '/public/images/icons/filter-dark.svg'
import searchDarkIcon from '/public/images/icons/search.svg'
import ProjectPI from './project-pi';
import FilterPill from '../searchable-list/filter-pill';
import { ValueOf } from '@/lib/utils/types';
import { FILE_EXPORT_STATUS, handleExportFile, MIME_FILE_FORMATS } from '@/lib/utils/file-export';


// TODO: clean this up
const searchableProjectKeys: (keyof SearchableProjectData)[] = ['fullName', 'principalInvestigators', 'type', 'dataTypes', 'theme', 'status']

interface SearchableProjectData extends Pick<Project, 'fullName' | 'principalInvestigators' | 'type' | 'dataTypes' | 'theme'> {
    status: ExternalProjectStatus
}

interface FilterItem {
    value: (ValueOf<Omit<SearchableProjectData, 'principalInvestigators' | 'dataTypes' | 'theme'>>) | PrincipalInvestigator | DataType
    type: (keyof Omit<SearchableProjectData, 'principalInvestigators' | 'dataTypes' | 'fullName'>) | 'principalInvestigator' | 'dataType' | 'name'
    label: string;
    key: string;
    projectKey: keyof SearchableProjectData;
}

const drawerFilterItemTypes: FilterItem['type'][] = ['theme', 'status', 'dataType']

const viewSets: FilterItem[] = Object.entries(ProjectHeadingInfoSet).map(([key, label]) => ({
    label: label,
    key: key,
    // Placeholder data just to satifsy the shape
    value: '',
    projectKey: 'fullName',
    type: 'name',
}))


// TODO: use separate list component in searchable-list module
export default function ProjectListings({
    projectData,
}: {
    projectData: Project[],
}) {
    const theme = useTheme()
    const isMobile = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.mobile,
        theme.breakpoints.values.tablet,
    ))

    const [currentPage, setCurrentPage] = React.useState(0)
    const pageSize = 12
    const [searchOptions] = React.useState(() => {

        const _searchOptions: Set<FilterItem> = new Set()

        const existingOptions = new Set<string>()

        function addOption(value: FilterItem['value'], type: FilterItem['type'], projectKey: keyof SearchableProjectData) {
            const key = getFilterItemKey(value, type)

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

    const [drawerItems] = React.useState(() => {
        return searchOptions
            .filter((option) => {
                return drawerFilterItemTypes.includes(option.type)
            })
    })

    const [drawerOpen, setDrawerOpen] = React.useState(false)

    const [filterItems, setFilterItems] = React.useState<FilterItem[]>([])

    const filteredProjectData = projectData
        .filter((project) => {
            if (filterItems.length === 0) return true;

            for (const filterItem of filterItems) {

                const projectKey = filterItem.projectKey
                let projectValue
                switch (projectKey) {
                    case 'fullName':
                    case 'type':
                        projectValue = project[projectKey]
                        break
                    case 'dataTypes':
                        for (const dt of project.dataTypes) {
                            if (dt === filterItem.value) {
                                return true
                            }
                        }
                        continue
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
                        continue
                }

                if (projectValue === filterItem.value) {
                    return true
                }
            }

            return false
        })

    const paginatedFilteredProjectData = filteredProjectData
        .slice(currentPage * pageSize, currentPage * pageSize + pageSize)

    const projectListRef = React.createRef<HTMLDivElement>()

    const [viewSet, setViewSet] = React.useState<(typeof viewSets)[number]>(viewSets[0])

    // handle book open state to coordinate only having one open at a time
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
    }

    const scrollToProject = (bookIdx: number) => {
        const projectEl = projectListRef.current?.children[bookIdx] as HTMLElement | null
        if (projectEl === null) {
            return
        }

        window.scroll({
            top: projectEl.offsetTop,
            behavior: 'smooth',
        })
    }

    const handleChangeFilterItems = (filterItems: FilterItem[]) => {
        setFilterItems(filterItems)
        closeAllProjects()
        setCurrentPage(0)
    }

    const handleClickRemoveFilterItem = (filterItem: FilterItem) => {
        const newFilterItems = filterItems.filter((item) => item.key !== filterItem.key)
        setFilterItems(newFilterItems)
        closeAllProjects()
        setCurrentPage(0)
    }

    // Manage file export
    const [fileExportStatus, setFileExportStatus] = React.useState(FILE_EXPORT_STATUS.idle)

    const handleClickExportButton = async () => {
        await handleExportFile(
            projectData,
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
                        rowGap: '16px',
                        [theme.breakpoints.up('tablet')]: {
                            display: 'flex',
                            flexDirection: 'row',
                            justifyContent: 'space-between',
                            rowGap: '38px',
                        }
                    }}
                >
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            columnGap: '8px',
                            rowGap: '16px',
                            [theme.breakpoints.up('tablet')]: {

                            },
                        }}
                    >
                        <Box
                            sx={{
                                display: 'flex',
                                flexDirection: 'row',
                                columnGap: '8px',
                                rowGap: '16px',
                                height: '52px',
                            }}
                        >
                            <DrawerButton
                                label='Filters'
                                icon={filterDarkIcon}
                                onClick={() => setDrawerOpen(!drawerOpen)}
                                activated={filterItems.find(item => drawerFilterItemTypes.includes(item.type)) !== undefined}
                                open={drawerOpen}
                            />

                            <Autocomplete
                                multiple
                                filterSelectedOptions
                                icon={searchDarkIcon}
                                placeholder='Search e.g. "Fibrosis"'
                                options={[...searchOptions].sort((a, b) => a.type.localeCompare(b.type))}
                                showResultsOnEmpty={false}
                                // @ts-ignore
                                onChange={(_, value: FilterItem[]) => handleChangeFilterItems(value)}
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
                                }}
                            />
                        </Box>

                        <Drawer
                            items={drawerItems}
                            viewSetItems={viewSets}
                            getItemLabel={item => item.label}
                            getItemKey={item => item.key}
                            getItemType={item => item.type}
                            activeItems={filterItems}
                            onChange={handleChangeFilterItems}
                            activeViewSetItem={viewSet}
                            onChangeViewSet={setViewSet}
                            showViewSets={isMobile}
                            open={drawerOpen}
                            showButton={isMobile}
                            buttonLabel='Export to CSV'
                            onClickButton={handleClickExportButton}
                            buttonDisabled={fileExportStatus === FILE_EXPORT_STATUS.inProgress}
                            displayStyle={isMobile ? 'expandable' : 'default'}
                        />

                        {filterItems.length > 0 && <Box
                            sx={{
                                display: 'flex',
                                flexDirection: 'row',
                                flexWrap: 'wrap',
                                columnGap: '10px',
                                rowGap: '16px',
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

                        <Box
                            sx={{
                                display: 'none',
                                [theme.breakpoints.up('tablet')]: {
                                    display: 'block',
                                },
                            }}
                        >
                            <ToggleGroup
                                valueIdx={viewSets.indexOf(viewSet)}
                                values={viewSets.map(val => val.label)}
                                onChange={(idx) => setViewSet(viewSets[idx])}
                            />
                        </Box>
                    </Box>

                    <Pagination
                        currentPage={currentPage}
                        pageSize={pageSize}
                        listSize={filteredProjectData.length}
                        onClickPrev={() => handleSetCurrentPage(currentPage - 1)}
                        onClickNext={() => handleSetCurrentPage(currentPage + 1)}
                        listItemLabel='projects'
                    />
                </Box>

                <Box
                    ref={projectListRef}
                    role='list'
                    sx={{
                        bgcolor: 'utilityWhite.main',
                        borderRadius: '20px',
                        overflow: 'hidden',
                    }}
                >
                    {paginatedFilteredProjectData.map((project, i) => (
                        <Box
                            key={project.fullName}
                            role='listitem'
                        >
                            <ProjectListing
                                data={project}
                                open={projectOpens[i]}
                                headingInfoSet={(viewSet.label as ProjectHeadingInfoSet)}
                                onSetOpen={open => handleSetProjectOpen(open, i)}
                                onFinishOpen={() => scrollToProject(i)}
                            />
                        </Box>
                    ))}
                </Box>
            </Box>
        </Container>
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
                <ProjectPI
                    data={option.value as PrincipalInvestigator}
                /> :
                <Typography variant='pMedium'>
                    {/* @ts-ignore */}
                    {option.value}
                </Typography>
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
            // @ts-ignore
            return value
        case 'principalInvestigator':
            return (value as PrincipalInvestigator).name
    }
}

function getFilterItemKey(value: FilterItem['value'], type: FilterItem['type']): string {
    return `${type}_${getFilterItemLabel(value, type)}`
}