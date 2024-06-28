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
                    role='list'
                >
                    {projectData.map(project => (
                        <Box
                            key={project.fullName}
                            role='listitem'
                        >
                            <ProjectListingItem
                                data={project}
                            />
                        </Box>
                    ))}
                </Box>
            </Box>
        </Container>
    )
}