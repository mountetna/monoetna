'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';

import { Project } from './models';
import ProjectListingItem from './project-listing-item';


export default function ProjectListings({
    projectData,
}: {
    projectData: Project[],
}) {
    const theme = useTheme()

    const projectListRef = React.createRef<HTMLDivElement>()

    // handle book open state to coordinate only having one open at a time
    const projectOpens = projectData.map((_, i) => React.useState(false))

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
                <Box>
                    controls placeholder
                </Box>

                <Box
                    ref={projectListRef}
                    role='list'
                    sx={{
                        bgcolor: 'utilityWhite.main',
                        borderRadius: '20px',
                    }}
                >
                    {projectData.map((project, i) => (
                        <Box
                            key={project.fullName}
                            role='listitem'
                        >
                            <ProjectListingItem
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