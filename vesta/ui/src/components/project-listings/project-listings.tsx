'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';

import { DataType, ExternalProjectStatus, getExternalProjectStatus, PrincipalInvestigator, Project } from './models';
import ProjectListing from './project-listing';
import Pagination from '@/components/searchable-list/controls/pagination'
import DrawerButton from '@/components/searchable-list/controls/drawer-button';
import Autocomplete from '@/components/searchable-list/controls/autocomplete';
import Tabs from '../inputs/tabs';
import { ThemeData } from '../themes/models';
import { ValueOf } from '@/lib/utils/types';

import filterDarkIcon from '/public/images/icons/filter-dark.svg'
import searchDarkIcon from '/public/images/icons/search.svg'


const collapsedInfoSets = [
    'Default',
    'Data Types',
    'Principal Investigators',
]

// TODO: clean this up
const searchableProjectKeys: (keyof SearchableProjectData)[] = ['fullName', 'principalInvestigators', 'type', 'dataTypes', 'theme', 'status']

interface SearchableProjectData extends Pick<Project, 'fullName' | 'principalInvestigators' | 'type' | 'dataTypes' | 'theme'> {
    status: ExternalProjectStatus
}

interface SearchOption {
    value: (ValueOf<Omit<SearchableProjectData, 'principalInvestigators' | 'dataTypes'>>) | PrincipalInvestigator | DataType
    type: (keyof Omit<SearchableProjectData, 'principalInvestigators' | 'dataTypes'>) | 'principalInvestigator' | 'dataType'
}


// TODO: use separate list component in searchable-list module
export default function ProjectListings({
    projectData,
}: {
    projectData: Project[],
}) {
    const theme = useTheme()

    const [currentPage, setCurrentPage] = React.useState(0)
    const pageSize = 12
    const [searchOptions] = React.useState(() => {

        const _searchOptions: Set<SearchOption> = new Set()

        const searchOptionsByType = {} as Record<keyof SearchableProjectData, Set<string>>
        for (const key of searchableProjectKeys) {
            searchOptionsByType[key] = new Set()
        }

        function addOption(option: SearchOption, projectKey: keyof SearchableProjectData) {
            const label = getOptionLabel(option)
            const existingOptions = searchOptionsByType[projectKey]

            if (!existingOptions.has(label)) {
                existingOptions.add(label)
                _searchOptions.add(option)
            }
        }

        for (const project of projectData) {
            for (const projectKey of searchableProjectKeys) {

                let value = project[projectKey]
                let type: SearchOption['type']

                switch (projectKey) {
                    case 'fullName':
                    case 'type':
                    case 'status':
                        type = projectKey
                        break
                    case 'theme':
                        value = (value as ThemeData).name
                        type = projectKey
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
                        addOption({
                            value: item,
                            type,
                        }, projectKey)
                    }
                } else {
                    addOption({
                        value: value,
                        type,
                    }, projectKey)
                }
            }
        }
        // Must sort by type in order for grouping to work
        return [..._searchOptions].sort((a, b) => a.type.localeCompare(b.type))
    })

    const [searchOptionValue, setSearchOptionValue] = React.useState<SearchOption[]>([])

    const filteredProjectData = projectData
        .slice(currentPage * pageSize, currentPage * pageSize + pageSize)

    const projectListRef = React.createRef<HTMLDivElement>()

    const [infoSetsIdx, setInfoSetsIdx] = React.useState<string | number | null>(0)

    // handle book open state to coordinate only having one open at a time
    // eslint-disable-next-line
    const projectOpens = projectData.map(_ => React.useState(false))

    const handleSetProjectOpen = (newOpenState: boolean, projectIdx: number) => {
        for (const [idx, [open, setOpen]] of projectOpens.entries()) {
            if (idx === projectIdx) {
                setOpen(newOpenState)
                continue
            }
            if (open) {
                setOpen(false)
            }
        }
    }

    const handleSetCurrentPage = (page: number) => {
        // close all open projects
        for (const [_, setOpen] of projectOpens) {
            setOpen(false)
        }
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
                            gap: '8px',
                            [theme.breakpoints.up('tablet')]: {

                            },
                        }}
                    >
                        <Box
                            sx={{
                                display: 'flex',
                                flexDirection: 'row',
                                gap: '8px',
                            }}
                        >
                            <DrawerButton
                                label='Filters'
                                icon={filterDarkIcon}
                                onClick={() => null}
                                activated={false}
                                open={false}
                            />

                            <Autocomplete
                                icon={searchDarkIcon}
                                placeholder='Search e.g. "Fibrosis"'
                                options={searchOptions}
                                multiple
                                freeSolo
                                onChange={(_, value) => console.log(value)}
                                value={searchOptionValue}
                                getOptionLabel={getOptionLabel}
                                renderOption={(params) => (
                                    <SearchOption
                                        option={params.option}
                                        index={params.index}
                                    />
                                )}
                                groupBy={option => option.type}
                                renderGroup={(params) => (
                                    <li key={params.key}>
                                        <h1>{params.group}</h1>
                                        <ul>{params.options}</ul>
                                    </li>
                                )}
                            />
                        </Box>

                        <Box
                            sx={{
                                display: 'none',
                                [theme.breakpoints.up('desktop')]: {
                                    display: 'block',
                                },
                            }}
                        >
                            <Tabs
                                valueIdx={infoSetsIdx}
                                values={collapsedInfoSets}
                                onChange={setInfoSetsIdx}
                            />
                        </Box>
                    </Box>

                    <Pagination
                        currentPage={currentPage}
                        pageSize={pageSize}
                        listSize={projectData.length}
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
                    {filteredProjectData.map((project, i) => (
                        <Box
                            key={project.fullName}
                            role='listitem'
                        >
                            <ProjectListing
                                data={project}
                                open={projectOpens[i][0]}
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
    index,
}: {
    option: SearchOption,
    index: number,
}) {
    return (
        <Box>
            type: {option.type}
            index: {index}
        </Box>
    )
}

function getOptionLabel(option: SearchOption | string): string {
    if (typeof option === 'string') {
        return option
    } else {
        let label: string

        switch (option.type) {
            case 'fullName':
            case 'type':
            case 'status':
            case 'dataType':
            case 'theme':
                // @ts-ignore
                label = option.value
                break
            case 'principalInvestigator':
                label = (option.value as PrincipalInvestigator).name
                break
        }

        return `${option.type}_${label}`
    }
}