'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useMediaQuery, useTheme } from '@mui/system';

import ThemeBookVertical from './theme-book/theme-book-vertical'
import ThemeBookHorizontal from './theme-book/theme-book-horizontal'
import { ThemeData } from './theme-book/shared';



export default function ThemeShelf({
    themeData,
}: {
    themeData: ThemeData[],
}) {
    const theme = useTheme()
    const isMobile = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.mobile,
        theme.breakpoints.values.tablet,
    ))

    const mobileThemeBookContainerRef = React.createRef<HTMLDivElement>()
    const tabletDesktopThemeBookContainerRef = React.createRef<HTMLDivElement>()

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
        const containerEl = isMobile ? mobileThemeBookContainerRef.current : tabletDesktopThemeBookContainerRef.current
        if (containerEl === null) {
            return
        }

        const bookEl = containerEl.children[bookIdx] as HTMLElement

        bookEl.scrollIntoView({
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
                ref={mobileThemeBookContainerRef}
                role='list'
                sx={(theme) => ({
                    display: 'flex',
                    flexDirection: 'column',
                    rowGap: '10px',
                    [theme.breakpoints.up('tablet')]: {
                        display: 'none',
                    },
                })}
            >
                {themeData.map((theme, i) => (
                    <ThemeBookVertical
                        key={theme.name}
                        data={theme}
                        onClickSeeProjects={handleClickSeeProjects}
                        open={bookOpens[i][0]}
                        onSetOpen={(open: boolean) => handleSetBookOpen(open, i)}
                        onFinishOpen={() => scrollToBook(i)}
                    />
                ))}
            </Box>

            <Box>
                <Box
                    ref={tabletDesktopThemeBookContainerRef}
                    role='list'
                    sx={(theme) => ({
                        display: 'none',
                        [theme.breakpoints.up('tablet')]: {
                            display: 'flex',
                            flexDirection: 'row',
                            columnGap: '10px',
                            width: '100%',
                            overflow: 'scroll',
                        },
                    })}
                >
                    {themeData.map((theme, i) => (
                        <ThemeBookHorizontal
                            key={theme.name}
                            data={theme}
                            onClickSeeProjects={handleClickSeeProjects}
                            open={bookOpens[i][0]}
                            onSetOpen={(open: boolean) => handleSetBookOpen(open, i)}
                            onFinishOpen={() => scrollToBook(i)}
                        />
                    ))}
                </Box>
            </Box>
        </Container>
    )
}

const handleClickSeeProjects = (event: React.MouseEvent<HTMLAnchorElement>, href: string) => {
    event.preventDefault()
    console.log(href)
    // TODO: scroll to filtered projects
}