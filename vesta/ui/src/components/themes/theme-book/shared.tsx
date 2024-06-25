import * as React from 'react'
import Box from '@mui/system/Box'
import { StaticImageData } from "next/image"


export interface ThemeData {
    name: string
    description: string
    project_count: number
    projects_link: string
    color: string
    textColor: 'light' | 'dark'
    image: StaticImageData
}


export interface ThemeBookProps {
    data: ThemeData,
    onClickSeeProjects: (event: React.MouseEvent<HTMLAnchorElement>, href: string) => void,
    open: boolean,
    onSetOpen: (newState: boolean) => void,
    onFinishOpen: () => void,
}


export function ProjectCount({
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