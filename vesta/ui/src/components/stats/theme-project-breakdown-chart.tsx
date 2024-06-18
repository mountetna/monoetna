import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useParentSize } from '@visx/responsive'
import Pie, { ProvidedProps, PieArcDatum } from '@visx/shape/lib/shapes/Pie';
import { Group } from '@visx/group';
import { scaleOrdinal } from '@visx/scale';
import { animated, useTransition, to, useSpring, UseTransitionProps } from '@react-spring/web';
import { alpha } from '@mui/material';


export interface ThemeData {
    name: string
    color: string
    project_count: number
}


export default function ThemeProjectBreakdownChart({
    data,
}: {
    data: ThemeData[],
}) {
    const [selectedTheme, setSelectedTheme] = React.useState<string | null>(null)

    const {
        parentRef: chartContainerRef,
        width: chartContainerWidth,
        height: chartContainerHeight,
    } = useParentSize({ debounceTime: 100, })

    let totalProjectCount = 0
    for (const theme of data) {
        totalProjectCount += theme.project_count
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
        trail: 100,
    }
    const donutLabelCount = selectedTheme ? data.find(val => val.name === selectedTheme)?.project_count : totalProjectCount
    const donutLabelCountTransitions = useTransition(donutLabelCount, transitionProps)
    const donutLabelText = selectedTheme ?? 'Projects'
    const donutLabelTextTransitions = useTransition(donutLabelText, transitionProps)

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
                        position: 'relative',
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
                                pieValue={d => d.project_count}
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