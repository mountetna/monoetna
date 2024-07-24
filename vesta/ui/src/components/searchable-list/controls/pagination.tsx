import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import PaginationArrows from '@/components/inputs/pagination-arrows';
import { SxProps } from '@mui/material';


export enum PaginationClasses {
    pageInfo = 'autocomplete-page-info',
    arrows = 'autocomplete-arrows',
}



export default function Pagination({
    currentPage,
    pageSize,
    listSize,
    onClickPrev,
    onClickNext,
    listItemLabel,
    showArrows = true,
    sx,
}: {
    currentPage: number,
    pageSize: number,
    listSize: number,
    onClickPrev: () => void,
    onClickNext: () => void,
    listItemLabel: string,
    showArrows?: boolean
    sx?: SxProps,
}) {
    const infoElements = [
        { text: 'Showing', type: 'basic' },
        { text: Math.min(currentPage * pageSize + 1, listSize), type: 'highlight' },
        { text: 'to', type: 'basic' },
        { text: Math.min(currentPage * pageSize + pageSize, listSize), type: 'highlight' },
        { text: 'of', type: 'basic' },
        { text: listSize, type: 'basic' },
        { text: listItemLabel, type: 'basic' },
    ]

    return (
        <Box
            sx={{
                display: 'flex',
                flexDirection: 'row',
                gap: '10px',
                justifyContent: 'space-between',
                alignItems: 'center',
                ...(sx || {})
            }}
        >
            {/* Page info */}
            <Box
                className='autocomplete autocomplete-page-info'
                sx={{
                    display: 'flex',
                    flexDirection: 'row',
                    justifyContent: 'center',
                    alignItems: 'center',
                    gap: '3px',
                    p: '7px 21px',
                    bgcolor: 'utilityWhite.main',
                    borderRadius: '38px',
                    '& .info': {
                        display: 'inline-block',
                        color: 'ground.grade10',
                    },
                    '& .highlight': {
                        color: 'blue.grade50',
                    },
                }}
            >
                {infoElements.map((el, i) => (
                    <Typography
                        key={i}
                        variant={el.type === 'highlight' ? 'pBodyBoldWt' : 'pBody'}
                        className={`info ${el.type}`}
                    >
                        {el.text}
                    </Typography>
                ))}
            </Box>

            {showArrows && <PaginationArrows
                className='autocomplete autocomplete-arrows'
                onClickPrev={onClickPrev}
                onClickNext={onClickNext}
                prevDisabled={currentPage === 0}
                nextDisabled={currentPage * pageSize + pageSize >= listSize}
            />}
        </Box>
    )
}