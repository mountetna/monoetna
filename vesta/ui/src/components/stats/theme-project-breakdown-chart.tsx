import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useParentSize } from '@visx/responsive'
import Pie from '@visx/shape/lib/shapes/Pie';
import { Group } from '@visx/group';
import { scaleOrdinal } from '@visx/scale';


export interface Project {
    name: string
    theme: string
}


export default function ThemeProjectBreakdownChart({
    themeColors,
    projectsByTheme,
}: {
    themeColors: Record<string, string>,
    projectsByTheme: Project[],
}) {
    const {
        parentRef: chartContainerRef,
        width: chartContainerWidth,
        height: chartContainerHeight,
    } = useParentSize({ debounceTime: 100, })

    const centerSize = 20
    const chartOuterRadius = Math.min(chartContainerWidth, chartContainerHeight) * 0.5 - (centerSize + 10)

    return (
        <Box
            className='theme-project-breakdown-chart'
            sx={(theme) => ({
                display: 'grid',
                p: '8px 8px 12px 8px',
                bgcolor: 'ground.grade10',
                borderRadius: '30px',
                [theme.breakpoints.up('tablet')]: {
                    gridTemplateColumns: 'repeat(12, 1fr)',
                },
                [theme.breakpoints.up('desktop')]: {
                    gridTemplateColumns: 'repeat(6, 1fr)',
                },
            })}
        >
            <Typography
                variant='h5'
                color='utilityWhite.main'
                sx={(theme) => ({
                    p: '8px',
                    [theme.breakpoints.up('tablet')]: {
                        gridColumn: 'span 12',
                    },
                    [theme.breakpoints.up('desktop')]: {
                        gridColumn: 'span 6',
                    },
                })}
            >
                Projects by Theme
            </Typography>
            <Box
                sx={(theme) => ({
                    display: 'flex',
                    justifyContent: 'center',
                    [theme.breakpoints.up('tablet')]: {
                        gridColumn: 'span 5',
                    },
                    [theme.breakpoints.up('desktop')]: {
                        gridColumn: 'span 3',
                    },
                })}
            >
                <Box
                    ref={chartContainerRef}
                    sx={(theme) => ({
                        width: '264px',
                        height: '264px',
                        [theme.breakpoints.up('tablet')]: {

                        },
                        [theme.breakpoints.up('desktop')]: {
                            width: '331px',
                            height: '331px',
                        },
                    })}
                >
                    <svg
                        width={chartContainerWidth}
                        height={chartContainerHeight}
                    >
                    </svg>
                </Box>
            </Box>
            <Box
                sx={(theme) => ({
                    display: 'none',
                    [theme.breakpoints.up('tablet')]: {
                        display: 'block',
                        gridColumn: 'span 7',
                    },
                    [theme.breakpoints.up('desktop')]: {
                        gridColumn: 'span 3',
                    },
                })}
            >

            </Box>
            <Box
                sx={(theme) => ({
                    display: 'block',
                    [theme.breakpoints.up('tablet')]: {
                        display: 'none',
                    },
                })}
            >
                
            </Box>
        </Box>
    )
}