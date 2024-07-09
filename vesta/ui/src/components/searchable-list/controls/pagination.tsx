import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import PaginationArrows from '@/components/inputs/pagination-arrows';


export default function Pagination({
    currentPage,
    pageSize,
    listSize,
    onClickPrev,
    onClickNext,
    listItemLabel,
    showArrows = true,
}: {
    currentPage: number,
    pageSize: number,
    listSize: number,
    onClickPrev: () => void,
    onClickNext: () => void,
    listItemLabel: string,
    showArrows?: boolean
}) {
    const infoElements = [
        { text: 'Showing', type: 'basic' },
        { text: currentPage * pageSize + 1, type: 'highlight' },
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
            alignItems: 'center',
        }}
        >
            {/* Page info */}
            <Box
                sx={{
                    display: 'flex',
                    flexDirection: 'row',
                    gap: '3px',
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
                onClickPrev={onClickPrev}
                onClickNext={onClickNext}
                prevDisabled={currentPage === 0}
                nextDisabled={currentPage * pageSize + pageSize >= listSize}
            />}
        </Box>
    )
}