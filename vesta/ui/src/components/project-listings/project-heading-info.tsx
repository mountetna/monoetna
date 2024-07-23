'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useMediaQuery, useTheme } from '@mui/material';

import { ProjectHeadingInfoSet, Project } from './models'
import Pill from '../pill/pill';
import ProjectPI from './project-pi';


export default function ProjectHeadingInfo({
    projectData,
    infoSet,
    variant = 'default',
}: {
    projectData: Project,
    infoSet: ProjectHeadingInfoSet,
    variant: 'default' | 'small',
}) {
    const theme = useTheme()
    const isMobile = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.mobile,
        theme.breakpoints.values.tablet,
    ))

    switch (infoSet) {
        case ProjectHeadingInfoSet.default:
            const itemGap = '6.16px'

            return (
                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'row',
                        gap: itemGap,
                        overflow: 'hidden',
                    }}
                >
                    <Typography
                        variant='h6'
                        sx={{
                            borderRadius: '40px',
                            bgcolor: projectData.theme.color,
                            color: projectData.theme.textColor === 'light' ? 'utilityHighlight.main' : 'ground.grade10',
                            px: '14px',
                            py: '3px',
                            overflow: 'hidden',
                            textOverflow: 'ellipsis',
                            textAlign: 'center',
                            textWrap: 'nowrap',
                        }}
                    >
                        {projectData.theme.name}
                    </Typography>

                    {variant === 'default' && (
                        <Box
                            sx={{
                                display: 'flex',
                                flexDirection: 'row',
                                gap: itemGap,
                                '& > *:not(:last-child)': {
                                    borderRight: '2px solid #D0D0D0',
                                },
                            }}
                        >
                            <Typography
                                variant='pTitleRegularWt'
                                sx={{
                                    width: '173px',
                                    pr: '2px',
                                    color: 'ground.grade10',
                                    textAlign: 'center',
                                }}
                            >
                                {projectData.type}
                            </Typography>

                            <Typography
                                variant='pTitleMonoCaps'
                                sx={{
                                    width: '144px',
                                    color: 'ground.grade10',
                                    textAlign: 'center',
                                }}
                            >
                                {projectData.name}
                            </Typography>
                        </Box>
                    )}
                </Box>
            )

        case ProjectHeadingInfoSet.dataTypes:
            const maxTypesShown = variant === 'default' ? 3 : 1
            const typesOverflowCount = Math.max(projectData.dataTypes.length - maxTypesShown, 0)

            return (
                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'row',
                        gap: '12.32px',
                        overflow: 'hidden',
                        '& *': {
                            overflow: 'hidden',
                            textOverflow: 'ellipsis',
                            textWrap: 'nowrap',
                        },
                    }}
                >
                    {projectData.dataTypes.slice(0, maxTypesShown).map(dt => (
                        <Pill
                            key={dt}
                            label={dt}
                            typographyVariant='pMediumMono'
                            variant='stroked'
                        />
                    ))}

                    {typesOverflowCount > 0 && (
                        <Box
                            key={`+${typesOverflowCount}dt`}
                            sx={{
                                overflow: 'visible',
                            }}
                        >
                            <Pill
                                key={`+${typesOverflowCount}`}
                                label={`+${typesOverflowCount}`}
                                typographyVariant='pMediumMono'
                                variant='stroked'
                            />
                        </Box>
                    )}
                </Box>
            )

        case ProjectHeadingInfoSet.pis:
            const maxPIsShown = variant === 'default' ? 3 : 1
            const pisOverflowCount = Math.max(projectData.principalInvestigators.length - maxPIsShown, 0)
            const piItemOverlap = '16px'

            return (
                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'row',
                        ml: piItemOverlap,
                        '& > *': {
                            ml: `-${piItemOverlap} !important`,
                            zIndex: 1,
                        },
                    }}
                >
                    {projectData.principalInvestigators.slice(0, maxPIsShown).map(pi => (
                        <ProjectPI
                            key={pi.name}
                            data={pi}
                            showAvatar={true}
                            showNameAndTitle={!isMobile && projectData.principalInvestigators.length === 1}
                        />
                    ))}

                    {pisOverflowCount > 0 && (
                        <Typography
                            key={`+${pisOverflowCount}pi`}
                            variant='p2XSBoldWt'
                            component='div'
                            sx={{
                                width: '41px',
                                height: '41px',
                                display: 'flex',
                                justifyContent: 'center',
                                alignItems: 'center',
                                bgcolor: 'magenta.grade50',
                                color: 'utilityWhite.main',
                                border: `1px solid ${theme.palette.utilityWhite.main}`,
                                borderRadius: '50%',
                            }}
                        >
                            <Box
                                sx={{
                                    marginLeft: '-2px',
                                }}
                            >
                                +{pisOverflowCount}
                            </Box>
                        </Typography>
                    )}
                </Box>
            )
    }
}