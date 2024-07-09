'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';

import { Project } from './models';
import ProjectListing from './project-listing';
import Pagination from '@/components/searchable-list/controls/pagination'


export default function ProjectListings({
    projectData,
}: {
    projectData: Project[],
}) {
    const theme = useTheme()

    const [currentPage, setCurrentPage] = React.useState(0)
    const pageSize = 12

    const filteredProjectData = projectData
        .slice(currentPage * pageSize, currentPage * pageSize + pageSize)

    const projectListRef = React.createRef<HTMLDivElement>()

    // handle book open state to coordinate only having one open at a time
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

    const handleSetCurrentPage = (page: number) => {
        // close all open projects
        handleSetProjectOpen(false, 0)
        setCurrentPage(page)
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
                        [theme.breakpoints.up('tablet')]: {
                            display: 'flex',
                            flexDirection: 'row',
                            justifyContent: 'space-between',
                        }
                    }}
                >
                    <Box>Search and filter placeholder</Box>

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