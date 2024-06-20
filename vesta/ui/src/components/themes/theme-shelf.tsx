'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';

import ThemeBook, { ThemeData } from './theme-book'



export default function ThemeShelf({
    themeData,
}: {
    themeData: ThemeData[],
}) {
    return (
        <Container>
            <Typography
                variant='h1'
                sx={{
                    mb: '16px',
                }}
            >
                Our Themes
            </Typography>

            <Typography
                variant='pLarge'
                component='div'
                sx={(theme) => ({
                    mb: '38px',
                    [theme.breakpoints.up('tablet')]: {
                        maxWidth: '659px',
                    },
                    [theme.breakpoints.up('desktop')]: {
                    },
                })}
            >
                Our research themes are focus areas driven to advance our understanding of challenges that face the human race.
            </Typography>

            <Box
                role='list'
                sx={(theme) => ({
                    display: 'flex',
                    flexDirection: 'column',
                    rowGap: '10px',
                })}
            >
                {themeData.map(theme => (
                    <ThemeBook
                        key={theme.name}
                        data={theme}
                    />
                ))}
            </Box>
        </Container>
    )
}