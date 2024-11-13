'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useMediaQuery, useTheme } from '@mui/system';
import { useRouter } from 'next/navigation';

import ThemeBookVertical from './theme-book/theme-book-vertical'
import ThemeBookHorizontal from './theme-book/theme-book-horizontal'
import { ThemeData } from './models';
import PaginationArrows from '../inputs/pagination-arrows';
import { useBreakpoint, useWindowDimensions } from '@/lib/utils/responsive';
import { containerPadding } from '@/theme';
import { scrollTo } from '@/lib/utils/scroll';


export default function ThemeShelf({
    themeData,
}: {
    themeData: ThemeData[],
}) {
    const theme = useTheme()
    const breakpoint = useBreakpoint()
    const isMobile = breakpoint === 'mobile'

    const mobileThemeBookContainerRef = React.createRef<HTMLDivElement>()
    const tabletDesktopThemeBookContainerRef = React.createRef<HTMLDivElement>()
    const tabletDesktopThemeBookListRef = React.createRef<HTMLDivElement>()

    // handle book open state to coordinate only having one open at a time
    const [bookOpens, setBookOpens] = React.useState(themeData.map((_, i) => i === 0))
    const [openBookIdx, setOpenBookIdx] = React.useState<number | null>(0)

    const handleSetBookOpen = (newOpenState: boolean, bookIdx: number) => {
        setBookOpens(opens => opens.map((_, idx) => {
            if (idx === bookIdx) {
                setOpenBookIdx(newOpenState === true ? idx : null)
                return newOpenState
            }
            return false
        }))
    }

    const tabletDesktopThemeBookGapPx = 10
    const tabletThemeBookFadeMaskWidthPx = 47
    const desktopThemeBookFadeMaskWidthPx = 170

    const scrollToBook = (bookIdx: number) => {
        const containerEl = isMobile ? mobileThemeBookContainerRef.current : tabletDesktopThemeBookListRef.current
        if (containerEl === null) {
            return
        }

        const bookEl = containerEl.children[bookIdx] as HTMLElement

        if (isMobile) {
            scrollTo({ top: bookEl.offsetTop }, breakpoint)
        } else {
            containerEl.scroll({
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

        setSpacerWidthPx(
            windowWidth > containerMaxWidth ?
                (windowWidth - containerMaxWidth) / 2 - tabletDesktopThemeBookGapPx
                : 0
        )
    }, [windowDimensions])


    const router = useRouter()

    const handleClickSeeProjects = (event: React.MouseEvent<HTMLAnchorElement>, href: string) => {
        event.preventDefault()
        const el = document.getElementById('projects')

        router.push(href, { scroll: false })
        if (el) {
            scrollTo({ top: el.offsetTop }, breakpoint)
        }
    }

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
                    Our research themes serve as broad categories that encompass the foundational biology of our projects.
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
                        open={bookOpens[i]}
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
                                open={bookOpens[i]}
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