import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useParentSize } from '@visx/responsive'
import Pie, { ProvidedProps, PieArcDatum } from '@visx/shape/lib/shapes/Pie';
import { Group } from '@visx/group';
import { scaleOrdinal } from '@visx/scale';
import { animated, useTransition, to, useSpring, UseTransitionProps } from '@react-spring/web';
import { alpha, useTheme } from '@mui/material';
import Image from 'next/image';

import { ThemeData } from '../themes/models';

import plusLight from '/public/images/icons/plus-light.svg'
import minusLight from '/public/images/icons/minus-light.svg'



/*
TODOS:
1. fixed height tablet on svg + legend container—not just svg
*/

export type ThemeProjectBreakdownData = Pick<ThemeData, 'name' | 'color' | 'projectCount'>


export default function ThemeProjectBreakdownChart({
    data,
}: {
    data: ThemeProjectBreakdownData[],
}) {
    const [selectedTheme, setSelectedTheme] = React.useState<string | null>(null)

    const [breakdownOpen, setBreakdownOpen] = React.useState(false)
    const breakdownContainerRef = React.useRef<HTMLElement>()
    const [breakdownStyle, animateBreakdownApi] = useSpring(
        () => ({ height: '0px', opacity: 0 }), []
    )

    const theme = useTheme()

    const toggleBreakdown = () => {
        setBreakdownOpen(!breakdownOpen)
        animateBreakdownApi.start({
            height: `${breakdownOpen ? 0 : breakdownContainerRef.current?.offsetHeight}px`,
            opacity: breakdownOpen ? 0 : 1,
            config: {
                easing: theme.transitions.easing.quintFn,
                duration: theme.transitions.duration.quint,
            },
        })
    }

    const {
        parentRef: chartContainerRef,
        width: chartContainerWidth,
        height: chartContainerHeight,
    } = useParentSize({ debounceTime: 500, })

    let totalProjectCount = 0
    for (const theme of data) {
        totalProjectCount += theme.projectCount
    }

    const themes = data.map(d => d.name)

    const getThemeColor = scaleOrdinal({
        domain: themes,
        range: data.map(d => d.color),
    })

    const radius = Math.min(chartContainerWidth, chartContainerHeight) / 2
    const centerY = chartContainerHeight / 2
    const centerX = chartContainerWidth / 2
    const donutThickness = 10

    const transitionProps: UseTransitionProps<any> = {
        from: { opacity: 0 },
        enter: { opacity: 1 },
        leave: { opacity: 0 },
    }
    const donutLabelCount = selectedTheme ? data.find(val => val.name === selectedTheme)?.projectCount : totalProjectCount
    const donutLabelCountTransitions = useTransition(donutLabelCount, transitionProps)
    const donutLabelText = selectedTheme ?? 'Projects'
    const donutLabelTextTransitions = useTransition(donutLabelText, transitionProps)

    const handleClickThemeItem = (themeName: string) => {
        if (themeName === selectedTheme) return setSelectedTheme(null)
        setSelectedTheme(themeName)
    }

    return (
        <Box
            className='theme-project-breakdown-chart'
            sx={(theme) => ({
                display: 'grid',
                alignItems: 'center',
                p: '8px 8px 12px 8px',
                bgcolor: 'ground.grade10',
                borderRadius: '30px',
                [theme.breakpoints.up('tablet')]: {
                    gridTemplateColumns: 'repeat(12, 1fr)',
                },
                [theme.breakpoints.up('desktop')]: {
                    gridTemplateColumns: 'repeat(6, 1fr)',
                    height: '100%',
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
                        position: 'relative',
                        width: '264px',
                        height: '264px',
                        my: '48px',
                        [theme.breakpoints.up('tablet')]: {
                            mx: '32px',
                            my: '32px',
                        },
                        [theme.breakpoints.up('desktop')]: {
                            width: '331px',
                            height: '331px',
                            mx: '32px',
                            my: '32px',
                        },
                    })}
                >
                    <svg
                        width={chartContainerWidth}
                        height={chartContainerHeight}
                        style={{
                            position: 'absolute',
                            top: 0,
                            left: 0,
                            zIndex: 1,
                        }}
                    >
                        <Group top={centerY} left={centerX}>
                            <Pie
                                data={data}
                                pieValue={d => d.projectCount}
                                outerRadius={radius}
                                innerRadius={radius - donutThickness}
                                cornerRadius={50}
                                padAngle={0.025}
                            >
                                {(pie) => (
                                    <AnimatedPie<ThemeData>
                                        {...pie}
                                        animate={true}
                                        getKey={(arc) => arc.data.name}
                                        onClickDatum={({ data: { name } }) =>
                                            setSelectedTheme(selectedTheme && selectedTheme === name ? null : name)
                                        }
                                        getColor={(arc) => {
                                            const baseColor = getThemeColor(arc.data.name)

                                            if (selectedTheme === null) return baseColor
                                            return selectedTheme === arc.data.name ? baseColor : alpha(baseColor, 0.25)
                                        }}
                                    />
                                )}
                            </Pie>
                        </Group>
                    </svg>
                    <Box
                        sx={{
                            width: '100%',
                            height: '100%',
                            display: 'flex',
                            flexDirection: 'column',
                            justifyContent: 'center',
                        }}
                    >
                        <Typography
                            variant='h1'
                            component='div'
                            sx={{
                                position: 'relative',
                                color: 'utilityWhite.main',
                                textAlign: 'center',
                                lineHeight: '1em',
                            }}
                        >
                            <Box
                                component='span'
                                sx={{
                                    minWidth: '7em',
                                    maxWidth: '10em',
                                    // maxWidth: '80%',
                                    width: '100%',
                                    visibility: 'hidden',
                                }}
                            >
                                {donutLabelCount}
                            </Box>
                            {donutLabelCountTransitions((style, item) => (
                                <animated.span style={{
                                    ...style,
                                    position: 'absolute',
                                    top: 0,
                                    left: '50%',
                                    transform: 'translateX(-50%)',
                                }}>
                                    {item}
                                </animated.span>
                            ))}
                        </Typography>
                        <Typography
                            variant='h5BoldWt'
                            component='div'
                            sx={{
                                position: 'relative',
                                color: 'utilityWhite.main',
                                textAlign: 'center',
                            }}
                        >
                            <Box
                                component='span'
                                sx={{
                                    minWidth: '7em',
                                    maxWidth: '10em',
                                    width: '100%',
                                    visibility: 'hidden',
                                }}
                            >
                                {donutLabelText}
                            </Box>
                            {donutLabelTextTransitions((style, item) => (
                                <animated.span style={{
                                    ...style,
                                    position: 'absolute',
                                    top: 0,
                                    left: '50%',
                                    transform: 'translateX(-50%)',
                                    width: '100%',
                                    minWidth: '7em',
                                    maxWidth: '10em',
                                    overflow: 'hidden',
                                    textOverflow: 'ellipsis',
                                }}>
                                    {item}
                                </animated.span>
                            ))}
                        </Typography>
                    </Box>
                </Box>
            </Box>
            <Box
                sx={(theme) => ({
                    display: 'none',
                    [theme.breakpoints.up('tablet')]: {
                        gridColumn: 'span 7',
                        display: 'flex',
                        flexWrap: 'wrap',
                        gap: '16px',
                        // alignSelf: 'start',
                    },
                    [theme.breakpoints.up('desktop')]: {
                        gridColumn: 'span 3',
                        flexDirection: 'column',
                    },
                    '& .theme-item': {
                        display: 'inline-flex',
                        width: 'fit-content',
                        '&.selected': {
                            bgcolor: 'ground.grade25',
                        },
                    },
                })}
            >
                {data.map(d => (
                    <ThemeItem
                        key={d.name}
                        name={d.name}
                        color={d.color}
                        count={d.projectCount}
                        selected={selectedTheme === d.name}
                        onClick={handleClickThemeItem}
                    />
                ))}
            </Box>
            <Box
                sx={(theme) => ({
                    display: 'block',
                    p: '16px',
                    borderRadius: '16px',
                    bgcolor: 'ground.grade25',
                    [theme.breakpoints.up('tablet')]: {
                        display: 'none',
                    },
                })}
            >
                <Typography
                    variant='h6BoldWt'
                    component='div'
                    color='ground.grade75'
                    sx={{
                        display: 'flex',
                        justifyContent: 'space-between',
                        alignItems: 'center',
                        px: '8px',
                        py: '6px',
                        cursor: 'pointer',
                    }}
                    onClick={toggleBreakdown}
                >
                    <span>{breakdownOpen ? 'Close' : 'Open'} breakdown</span>
                    <Box
                        sx={{
                            display: 'flex',
                            p: '6px',
                            borderRadius: '50%',
                            bgcolor: 'ground.grade10',
                        }}
                    >
                        <Image
                            src={breakdownOpen ? minusLight : plusLight}
                            alt={breakdownOpen ? 'Minus sign' : 'Plus sign'}
                            width={27}
                            height={27}
                        />
                    </Box>
                </Typography>
                <animated.div
                    style={{
                        overflow: 'hidden',
                        ...breakdownStyle,
                    }}
                >
                    <Box
                        ref={breakdownContainerRef}
                        sx={{
                            overflow: 'hidden',
                            pt: '16px',
                            '& > *:not(:last-child)': {
                                mb: '16px',
                            },
                        }}
                    >
                        {data.map(d => (
                            <ThemeItem
                                key={d.name}
                                name={d.name}
                                color={d.color}
                                count={d.projectCount}
                                selected={selectedTheme === d.name}
                                onClick={handleClickThemeItem}
                            />
                        ))}
                    </Box>
                </animated.div>
            </Box>
        </Box>
    )
}


function ThemeItem({
    name,
    color,
    count,
    selected,
    onClick,
}: {
    name: string,
    color: string,
    count: number,
    selected: boolean,
    onClick: (name: string) => void,
}) {
    const theme = useTheme()

    return (
        <Box
            className={`theme-item${selected ? ' selected' : ''}`}
            sx={{
                display: 'flex',
                alignItems: 'center',
                px: '8px',
                py: '6px',
                borderRadius: '8px',
                color: 'ground.grade75',
                cursor: 'pointer',
                transition: theme.transitions.create(
                    ['background-color'],
                    {
                        duration: theme.transitions.duration.ease,
                        easing: theme.transitions.easing.ease,
                    },
                ),
                '&.selected': {
                    bgcolor: 'ground.grade10',
                },
                '& > *:not(:last-child)': {
                    mr: '11px',
                },
            }}
            onClick={_ => onClick(name)}
        >
            <Box
                sx={{
                    width: '34px',
                    height: '17px',
                    bgcolor: color,
                    borderRadius: '4px',
                }}
            />
            <Typography
                variant='pLargeBoldWt'
                component='div'
            >
                {name}
            </Typography>
            <Typography
                variant='pLarge'
                component='div'
            >
                {count}
            </Typography>
        </Box>
    )
}


type AnimatedStyles = { startAngle: number; endAngle: number; opacity: number };

const fromLeaveTransition = ({ endAngle }: PieArcDatum<any>) => ({
    // enter from 360° if end angle is > 180°
    startAngle: endAngle > Math.PI ? 2 * Math.PI : 0,
    endAngle: endAngle > Math.PI ? 2 * Math.PI : 0,
    opacity: 0,
});
const enterUpdateTransition = ({ startAngle, endAngle }: PieArcDatum<any>) => ({
    startAngle,
    endAngle,
    opacity: 1,
});

type AnimatedPieProps<Datum> = ProvidedProps<Datum> & {
    animate?: boolean;
    getKey: (d: PieArcDatum<Datum>) => string;
    getColor: (d: PieArcDatum<Datum>) => string;
    onClickDatum: (d: PieArcDatum<Datum>) => void;
    delay?: number;
};

function AnimatedPie<Datum>({
    animate,
    arcs,
    path,
    getKey,
    getColor,
    onClickDatum,
}: AnimatedPieProps<Datum>) {
    const transitions = useTransition<PieArcDatum<Datum>, AnimatedStyles>(arcs, {
        from: animate ? fromLeaveTransition : enterUpdateTransition,
        enter: enterUpdateTransition,
        update: enterUpdateTransition,
        leave: animate ? fromLeaveTransition : enterUpdateTransition,
        keys: getKey,
    });

    return transitions((props, arc, { key }) => {
        const animatedFill = useSpring({
            fill: getColor(arc)
        })

        return (
            <g key={key} style={{ cursor: 'pointer' }}>
                <animated.path
                    // compute interpolated path d attribute from intermediate angle values
                    d={to([props.startAngle, props.endAngle], (startAngle, endAngle) =>
                        path({
                            ...arc,
                            startAngle,
                            endAngle,
                        }),
                    )}
                    style={animatedFill}
                    onClick={() => onClickDatum(arc)}
                    onTouchStart={() => onClickDatum(arc)}
                />
            </g>
        );
    });
}