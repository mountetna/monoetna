import * as React from 'react'
import Image, { StaticImageData } from "next/image"
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { SxProps, Theme } from '@mui/material';
import { Variant } from '@mui/material/styles/createTypography';


export interface ThemeData {
    name: string
    description: string
    project_count: number
    projects_link: string
    color: string
    image: StaticImageData
}


export default function ThemeBook({
    data,
}: {
    data: ThemeData,
}) {
    return (
        <Box
            role='listitem'
            sx={(theme) => ({
                display: 'flex',
                flexDirection: 'column',
                rowGap: '24px',
                p: '16px',
                bgcolor: 'utilityWhite.main',
                borderRadius: '20px',
                cursor: 'pointer',
            })}
        >
            {/* HEADING */}
            <Box
                sx={(theme) => ({
                    display: 'flex',
                    justifyContent: 'space-between',
                })}
            >
                <Typography
                    variant='pSubtitleMonoCaps'
                    component='span'
                >
                    {data.name.slice(0, 4)}
                </Typography>
                <Typography
                    variant='pSubtitleMonoCaps'
                    component='div'
                    sx={(theme) => ({
                        display: 'none',
                        [theme.breakpoints.up('tablet')]: {
                            display: 'inline-block',
                        },
                    })}
                >
                    <ProjectCount
                        count={data.project_count}

                    />
                </Typography>
                <Typography
                    variant='h5'
                    component='span'
                >
                    {data.name}
                </Typography>
            </Box>

            {/* MAIN CONTENT */}
            <Box>
                <Box
                    sx={(theme) => ({
                        p: '16px',
                        borderRadius: '20px',
                        bgcolor: data.color,
                    })}
                >
                    <Typography
                        variant='pMediumMonoCaps'
                        component='div'
                        sx={(theme) => ({
                            [theme.breakpoints.up('tablet')]: {
                                display: 'none',
                            },
                        })}
                    >
                        <ProjectCount
                            count={data.project_count}
                        />
                    </Typography>
                </Box>
            </Box>
        </Box>
    )
}


function ProjectCount({
    count,
}: {
    count: number,
}) {
    return (
        <Box
            component='span'
            sx={{
                display: 'inline-flex',
                columnGap: '10px',
                py: '10px',
                borderRadius: '30px',
                color: 'utilityWhite.main',
                bgcolor: 'ground.grade10',
                writingMode: 'vertical-lr',
                textOrientation: 'mixed',
                rotate: '180deg',
            }}
        >
            <Box component='span'>
                {count}
            </Box>
            <Box component='span'>
                {`PROJECT${count > 1 ? 'S' : ''}`}
            </Box>
        </Box>
    )
}