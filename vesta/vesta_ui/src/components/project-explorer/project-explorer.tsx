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

import { ExternalProjectStatus, getExternalProjectStatus, FilterItem, PrincipalInvestigator, Project, ProjectHeadingInfoSet, PROJECTS_SEARCH_PARAMS_KEY, ProjectsSearchParamsControls, ProjectsSearchParamsState } from './models';
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

const drawerFilterItemTypes: FilterItem['type'][] = ['theme', 'status', 'dataType']

const PIExportAttrs: (keyof PrincipalInvestigator)[] = ['name', 'title']
const ThemeExportAttrs: (keyof ThemeData)[] = ['name', 'description', 'projectsLink']


function _ProjectExplorer({ }) {
    // Manage search params sync
    const { state: { projectData, filters, filterItemSet }, searchOptions, updateFilterItems } = React.useContext(ProjectExplorerContext);
    const router = useRouter()
    const pathname = usePathname()
    const searchParams = useSearchParams()

    const theme = useTheme()
    const breakpoint = useBreakpoint()
    const isMobile = breakpoint === 'mobile'

    const searchPlaceholder = 'Search e.g. "Fibrosis"'

    const [drawerOpen, setDrawerOpen] = React.useState(true);
    const [inputValue, setInputValue] = React.useState('');
    const [currentPage, setCurrentPage] = React.useState(0);

    const handleChangeFilterItems = React.useCallback(
      ([ filterItem, ...others ]: FilterItem[]) => {
        if (!(filterItem.type in filterItemSet) || !filterItemSet[filterItem.type].includes(filterItem.value))
        updateFilterItems(filterItem.type, (filterItemSet[filterItem.type] || []).concat(filterItem.value))
        setCurrentPage(0)
      }, [ filterItemSet, updateFilterItems ]
    );

    return (
        <Container
            id='project-explorer'
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
                                    activated={ Object.keys(filterItemSet).length > 0 }
                                    open={drawerOpen}
                                />

                                <Autocomplete
                                    multiple
                                    icon={searchDarkIcon}
                                    placeholder={searchPlaceholder}
                                    options={[...searchOptions].sort((a, b) => a.type.localeCompare(b.type))}
                                    showResultsOnEmpty={false}
                                    // @ts-ignore
                                    inputValue={inputValue}
                                    onInputChange={ (_, value: string, reason) => setInputValue(value) }
                                    onKeyDown={
                                      (e) => {
                                        if (e.key === 'Enter') {
                                          if (inputValue) {
                                            if ('free' in filterItemSet && filterItemSet['free'].includes(inputValue)) return;
                                            updateFilterItems( 'free', (filterItemSet['free'] || []).concat(inputValue) );
                                          }
                                          setInputValue('');
                                        }
                                      }
                                    }
                                    // @ts-ignore
                                    onChange={(_, value: FilterItem[], reason) => {
                                        // disables removing options with backspace
                                        if (reason === 'removeOption') return

                                        handleChangeFilterItems(value)
                                    }}
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
                                                {params.group}
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
                    </Box>
                </Box>

                <Box sx={{
                  display: 'flex'
                }}>
                  {drawerOpen && <Box
                    sx={
                      theme => ({
                        flex: '0 0 405px',
                        height: '640px',
                        overflowY: 'scroll',
                        pt: '9px',
                        pr: '5px',
                        borderRight: `1px solid ${theme.palette.ground.grade50}`,
                        // @ts-ignore
                        transition: theme.transitions ? theme.transitions.create(
                          ['flex'],
                          {
                        // @ts-ignore
                             easing: theme.transitions.easing.ease,
                        // @ts-ignore
                             duration: theme.transitions.duration.ease,
                          }
                        ) : 'unset',
                        borderEndEndRadius: '2px'
                      })
                     }>
                     <FilterList/>
                    </Box>
                  }
                  <Box
                      role='list'
                      sx={{ flex: '1 1 auto', overflowX: 'hidden' }}
                  >
                    <Box>
                      <ProjectTable currentPage={currentPage} setCurrentPage={setCurrentPage}/>
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
            {option.type === 'Investigators' ?
                (
                    <ProjectPI
                        data={option.value as PrincipalInvestigator}
                        variant='filled'
                    />
                ) : (
                    <Typography variant='pMedium'>
                        {/* @ts-ignore */}
                        {typeof(option.value) === 'object' ? option.value.name : option.value}
                    </Typography>
                )
            }

        </Box>
    )
}
