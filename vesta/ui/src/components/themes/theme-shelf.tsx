'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useMediaQuery, useTheme } from '@mui/system';

import ThemeBookVertical from './theme-book/theme-book-vertical'
import ThemeBookHorizontal from './theme-book/theme-book-horizontal'
import { ThemeData } from './models';
import PaginationArrows from '../inputs/pagination-arrows';
import { useWindowDimensions } from '@/lib/utils/responsive';
import { containerPadding } from '@/theme';


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
    const tabletDesktopThemeBookListRef = React.createRef<HTMLDivElement>()

    // handle book open state to coordinate only having one open at a time
    const bookOpens = themeData.map((_, i) => React.useState(i === 0))
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

    const tabletDesktopThemeBookGapPx = 10
    const tabletThemeBookFadeMaskWidthPx = 47
    const desktopThemeBookFadeMaskWidthPx = 220

    const scrollToBook = (bookIdx: number) => {
        const containerEl = isMobile ? mobileThemeBookContainerRef.current : tabletDesktopThemeBookListRef.current
        if (containerEl === null) {
            return
        }

        const bookEl = containerEl.children[bookIdx] as HTMLElement

        if (isMobile) {
            window.scroll({
                top: bookEl.offsetTop,
                behavior: 'smooth',
            })
        } else {
            containerEl.scroll({
                // TODO: include navbar offset
                left: bookEl.offsetLeft,
                behavior: 'smooth',
            })
        }
    }

    // handle positioning ThemeBook for max-width desktop
    const { dimensions: windowDimensions } = useWindowDimensions()
    const [spacerWidthPx, setSpacerWidthPx] = React.useState<number>(0)

    React.useEffect(() => {
        if (windowDimensions === null) return

        const windowWidth = windowDimensions[0]
        const containerMaxWidth = theme.breakpoints.values.desktopLg

        setSpacerWidthPx(windowWidth > containerMaxWidth ?
            (windowWidth - containerMaxWidth) / 2 - tabletDesktopThemeBookGapPx
            : 0
        )
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
                {themeData.map((theme, i) => (
                    <ThemeBookVertical
                        key={theme.name}
                        data={theme}
                        onClickSeeProjects={handleClickSeeProjects}
                        open={bookOpens[i][0]}
                        onSetOpen={open => handleSetBookOpen(open, i)}
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
                    sx={(theme) => ({
                        position: 'relative',
                        '& .fade-mask': {
                            display: 'none',
                            [theme.breakpoints.up('tablet')]: {
                                '&.fade-mask-right': {
                                    display: 'block',
                                },
                                position: 'absolute',
                                zIndex: 1,
                                top: 0,
                                width: `${tabletThemeBookFadeMaskWidthPx}px`,
                                height: '100%',
                            },
                            [theme.breakpoints.up(theme.breakpoints.values.desktopLg + desktopThemeBookFadeMaskWidthPx + 60)]: {
                                '&.fade-mask-left': {
                                    display: 'block',
                                },
                                width: `${desktopThemeBookFadeMaskWidthPx}px`
                            },
                        },
                    })}
                >
                    <Box
                        className='fade-mask fade-mask-left'
                        sx={(theme) => ({
                            left: 0,
                            background: `linear-gradient(to left, transparent, ${theme.palette.utilityHighlight.main})`,
                        })}
                    />

                    <Box
                        ref={tabletDesktopThemeBookListRef}
                        role='list'
                        sx={{
                            display: 'flex',
                            flexDirection: 'row',
                            columnGap: `${tabletDesktopThemeBookGapPx}px`,
                            // justifyContent: 'space-between',
                            width: '100%',
                            overflowX: 'scroll',
                            ...containerPadding,
                            pt: '30px !important',
                            '& .spacer': {
                                minWidth: `${spacerWidthPx}px`,
                                minHeight: '100%',
                            },
                        }}
                    >
                        <Box className='spacer' />

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

                        <Box className='spacer' />
                    </Box>

                    <Box
                        className='fade-mask fade-mask-right'
                        sx={(theme) => ({
                            right: 0,
                            background: `linear-gradient(to right, transparent, ${theme.palette.utilityHighlight.main})`
                        })}
                    />
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
                            nextDisabled={openBookIdx === themeData.length - 1}
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