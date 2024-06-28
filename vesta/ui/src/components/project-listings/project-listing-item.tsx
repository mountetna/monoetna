'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import Image from 'next/image';

import { Project } from "./models";


export default function ProjectListingItem({
    data,
}: {
    data: Project,
}) {
    const theme = useTheme()

    return (
        <Box
            className='project-listing-item'
        >
            {/* HEADING */}
            <Box
                sx={{
                    display: 'flex',
                    flexDirection: 'row',
                    gap: '10px'
                }}
            >
                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'row',
                        gap: '6px',
                    }}
                >
                    <Image
                        src={data.theme.icon}
                        alt={`Abstract icon for "${data.theme.name}" theme`}
                        style={{
                            borderRadius: '50%',
                            backgroundColor: data.theme.color,
                        }}
                    />

                    <Typography variant='h6'>
                        {data.fullName}
                    </Typography>
                </Box>

                <Typography
                variant='h6'
                sx={{
                    borderRadius: '40px',
                    bgcolor: data.theme.color,
                    px: '14px',
                    py: '3px',
                }}
                >
                    {data.theme.name}
                </Typography>
            </Box>
        </Box>
    )
}