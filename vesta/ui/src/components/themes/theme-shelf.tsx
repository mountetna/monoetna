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
    const themeBookContainerRef = React.createRef<HTMLDivElement>()

    // handle book open state to coordinate only having one open at a time
    const bookOpens = themeData.map(() => React.useState(false))

    const handleSetBookOpen = (newOpenState: boolean, bookIdx: number) => {
        for (const [idx, [open, setOpen]] of bookOpens.entries()) {
            if (idx === bookIdx) {
                setOpen(newOpenState)
                continue
            }
            if (open) {
                setOpen(false)
            }
        }
    }

    const scrollToBook = (bookIdx: number) => {
        const containerEl = themeBookContainerRef.current
        if (containerEl === null) {
            return
        }

        const bookEl = containerEl.children[bookIdx]

        window.scrollTo({
            top: bookEl?.offsetTop,
            behavior: 'smooth',
        })
    }

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
                ref={themeBookContainerRef}
                role='list'
                sx={(theme) => ({
                    display: 'flex',
                    flexDirection: 'column',
                    rowGap: '10px',
                })}
            >
                {themeData.map((theme, i) => (
                    <ThemeBook
                        key={theme.name}
                        data={theme}
                        onClickSeeProjects={handleClickSeeProjects}
                        open={bookOpens[i][0]}
                        onSetOpen={(open: boolean) => handleSetBookOpen(open, i)}
                        onFinishOpen={() => scrollToBook(i)}
                    />
                ))}
            </Box>
        </Container>
    )
}

const handleClickSeeProjects = (event: React.MouseEvent<HTMLAnchorElement>, href: string) => {
    event.preventDefault()
    console.log(href)
}