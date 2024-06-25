'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useMediaQuery, useTheme } from '@mui/system';

import ThemeBookVertical from './theme-book/theme-book-vertical'
import ThemeBookHorizontal from './theme-book/theme-book-horizontal'
import { ThemeData } from './theme-book/shared';
import PaginationArrows from '../inputs/pagination-arrows';
import { useWindowDimensions } from '@/lib/utils/responsive';
import { containerPadding } from '@/theme';



export default function ThemeShelf({
    themeData,
}: {
    themeData: ThemeData[],
}) {
    const _themeData = [...themeData]
    _themeData.push(..._themeData.slice(0, 8))

    const theme = useTheme()
    const isMobile = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.mobile,
        theme.breakpoints.values.tablet,
    ))

    const mobileThemeBookContainerRef = React.createRef<HTMLDivElement>()
    const tabletDesktopThemeBookContainerRef = React.createRef<HTMLDivElement>()
    const tabletDesktopThemeBookSpacerRef = React.createRef<HTMLDivElement>()
    const tabletDesktopThemeBookListRef = React.createRef<HTMLDivElement>()

    // handle book open state to coordinate only having one open at a time
    const bookOpens = _themeData.map((_, i) => React.useState(i === 0))
    const [openBookIdx, setOpenBookIdx] = React.useState<number | null>(0)

    const handleSetBookOpen = (newOpenState: boolean, bookIdx: number) => {
        for (const [idx, [open, setOpen]] of bookOpens.entries()) {
            if (idx === bookIdx) {
                setOpen(newOpenState)
                setOpenBookIdx(newOpenState === true ? idx : null)
                continue
            }
            if (open) {
                setOpen(false)
            }
        }
    }

    const scrollToBook = (bookIdx: number) => {
        const containerEl = isMobile ? mobileThemeBookContainerRef.current : tabletDesktopThemeBookListRef.current
        if (containerEl === null) {
            return
        }

        const bookEl = containerEl.children[bookIdx] as HTMLElement

        containerEl.scroll({
            left: bookEl.offsetLeft,
            top: bookEl.offsetTop,
            behavior: 'smooth',
        })
    }

    const tabletDesktopThemeBookGapPx = 10

    // handle positioning ThemeBook for max-width desktop
    const {
        dimensions: windowDimensions,
        isResizing: isWindowResizing,
    } = useWindowDimensions()

    React.useEffect(() => {
        const themeBookSpacerRef = tabletDesktopThemeBookSpacerRef.current
        if (themeBookSpacerRef === null || windowDimensions === null) return

        const windowWidth = windowDimensions[0]
        const containerMaxWidth = theme.breakpoints.values.desktopLg

        if (windowWidth > containerMaxWidth) {
            themeBookSpacerRef.style.display = 'block'
            themeBookSpacerRef.style.minWidth = `${(windowWidth - containerMaxWidth) / 2 - tabletDesktopThemeBookGapPx}px`
        } else {
            themeBookSpacerRef.style.display = 'none'
        }
    }, [windowDimensions])

    return (
        <Box>
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
                            mb: '18px',
                        },
                        [theme.breakpoints.up('desktop')]: {
                        },
                    })}
                >
                    Our research themes are focus areas driven to advance our understanding of challenges that face the human race.
                </Typography>
            </Container>

            <Container
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
                {_themeData.map((theme, i) => (
                    <ThemeBookVertical
                        key={theme.name}
                        data={theme}
                        onClickSeeProjects={handleClickSeeProjects}
                        open={bookOpens[i][0]}
                        onSetOpen={(open: boolean) => handleSetBookOpen(open, i)}
                        onFinishOpen={() => scrollToBook(i)}
                    />
                ))}
            </Container>

            <Box
                sx={(theme) => ({
                    display: 'none',
                    [theme.breakpoints.up('tablet')]: {
                        display: 'flex',
                        flexDirection: 'column',
                        gap: '21px',
                    },
                })}
            >
                <Box
                    ref={tabletDesktopThemeBookContainerRef}
                >
                    <Box
                        ref={tabletDesktopThemeBookListRef}
                        role='list'
                        sx={{
                            display: 'flex',
                            flexDirection: 'row',
                            columnGap: `${tabletDesktopThemeBookGapPx}px`,
                            width: '100%',
                            overflowX: 'scroll',
                            // p: '10px',
                            pt: '30px',
                        }}
                    >
                        <Box
                            ref={tabletDesktopThemeBookSpacerRef}
                            sx={{
                                minHeight: '100%',
                            }}
                        />
                        <Typography
                            variant='h4SmallCaps'
                            component='div'
                            sx={{
                                color: 'ground.grade25',
                                writingMode: 'vertical-lr',
                                textOrientation: 'mixed',
                                rotate: '180deg',
                                pl: '10px',
                                borderLeft: '1px solid black',
                            }}
                        >
                            THEMES
                        </Typography>

                        {_themeData.map((theme, i) => (
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

                <Container>
                    <Box
                        sx={{
                            display: 'flex',
                            gap: '16px',
                            alignItems: 'center',
                            justifyContent: 'flex-end',
                        }}
                    >
                        <Typography
                            variant='h4SmallCaps'
                            color='ground.grade25'
                            component='div'
                        >
                            EXPLORE OUR THEMES
                        </Typography>
                        <PaginationArrows
                            onClickPrev={() => openBookIdx !== null && handleSetBookOpen(true, openBookIdx - 1)}
                            prevDisabled={openBookIdx === null || openBookIdx === 0}
                            onClickNext={() => handleSetBookOpen(true, openBookIdx === null ? 0 : openBookIdx + 1)}
                            nextDisabled={openBookIdx === _themeData.length - 1}
                        />
                    </Box>
                </Container>
            </Box>
        </Box>
    )
}

const handleClickSeeProjects = (event: React.MouseEvent<HTMLAnchorElement>, href: string) => {
    event.preventDefault()
    console.log(href)
    // TODO: scroll to filtered projects
}