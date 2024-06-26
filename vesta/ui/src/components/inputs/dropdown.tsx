'use client'

import * as React from 'react';
import Box from '@mui/system/Box'
import { Select } from '@mui/base/Select';
import { Option as BaseOption, optionClasses } from '@mui/base/Option';
import Typography, { TypographyOwnProps } from '@mui/material/Typography';
import { useParentSize } from '@visx/responsive'
import _ from 'lodash'
import Image from 'next/image';

import indicatorArrowDark from '/public/images/icons/indicator-arrow-dark.svg'
import { styled } from '@mui/material';


interface OptionType {
    value: any
    label: string
}

const IndicatorArrowDark = (
    <Box
        sx={{
            display: 'flex',
        }}
        component='span'
    >
        <Image
            src={indicatorArrowDark}
            alt='Indicator arrow pointing down'
        />
    </Box>
)

export default function Dropdown({
    options,
    value,
    onChange,
    icon = IndicatorArrowDark,
    listboxStyles = {},
    optionStyles = {},
    optionTypographyVariant = 'pBodyMediumWt',
    selectedOptionTypographyVariant = 'pBody',
}: {
    options: OptionType[],
    value: any,
    onChange: (newValue: any) => void,
    icon?: React.ReactNode,
    listboxStyles?: React.CSSProperties,
    optionStyles?: React.CSSProperties,
    optionTypographyVariant?: TypographyOwnProps['variant'],
    selectedOptionTypographyVariant?: TypographyOwnProps['variant'],
}) {
    const {
        parentRef: selectContainerRef,
        width: selectContainerWidth,
    } = useParentSize({ debounceTime: 100, })

    const Listbox = styled('ul')(
        ({ theme }) => `
        box-sizing: border-box;
        width: ${selectContainerWidth}px;
        padding: 10px;
        margin: 10px 0;
        background-color: ${theme.palette.utilityWhite.main};
        border: 1px solid ${theme.palette.ground.grade75};
        border-radius: 16px;
        overflow: auto;
        list-style: none;
        
        ${Object.entries(listboxStyles)
                .map(([k, v]) => `${_.kebabCase(k)}: ${v};`)
                .join('\n')
            }
        `,
    )

    const Option = styled(BaseOption)(
        ({ theme }) => `
        list-style: none;
        padding: 8px;
        cursor: default;
        margin-bottom: 8px;
        border-radius: 8px;
        color: ${theme.palette.ground.grade10};
      
        &:last-of-type {
          margin-bottom: 0;
        }
      
        &.${optionClasses.selected} {
          background-color: ${theme.palette.blue.grade50};
          color: ${theme.palette.utilityWhite.main};
        }
      
        &.${optionClasses.highlighted} {
          background-color: ${theme.palette.blue.grade50};
          color: ${theme.palette.utilityWhite.main};
        }
      
        &.${optionClasses.highlighted}.${optionClasses.selected} {
          background-color: ${theme.palette.blue.grade50};
          color: ${theme.palette.utilityWhite.main};
        }
      
        &:focus-visible {
          outline: 3px solid ${theme.palette.blue.grade50};
        }
      
        &.${optionClasses.disabled} {
          color: unset;
        }

        &:hover:not(.${optionClasses.disabled}) {
          cursor: pointer;
        }
      
        &:hover:not(.${optionClasses.disabled}, .${optionClasses.selected}, ${optionClasses.highlighted}) {
          background-color: ${theme.palette.utilityHighlight.main};
        }

        ${Object.entries(optionStyles)
                .map(([k, v]) => `${_.kebabCase(k)}: ${v};`)
                .join('\n')
            }
        `,
    )

    return (
        <Box
            className='dropdown'
            sx={{
                bgcolor: 'utilityHighlight.main',
                borderRadius: '10px',
                '& > button:hover': {
                    cursor: 'pointer',
                },
            }}
            ref={selectContainerRef}
        >
            <Select
                className='dropdown-select'
                slots={{ popup: Popup, listbox: Listbox }}
                value={value}
                onChange={(_, newValue) => onChange(newValue)}
                style={{
                    width: '100%',
                    border: 'none',
                    borderRadius: '10px',
                }}
                renderValue={(option) => (
                    <Box
                        className='dropdown-button-content'
                        sx={{
                            display: 'flex',
                            justifyContent: 'space-between',
                            alignItems: 'center',
                            px: '16px',
                            py: '10px',
                        }}
                    >
                        <Typography
                            variant={selectedOptionTypographyVariant}
                            component='span'
                            sx={{
                                mr: '10px',
                                color: 'utilityGround.grade10',
                            }}
                        >
                            {option?.label}
                        </Typography>
                        {icon}
                    </Box>
                )}
            >
                {options.map(option => (

                    <Option
                        key={option.value}
                        value={option.value}
                    >
                        <Typography
                            variant={optionTypographyVariant}
                        >
                            {option.label}
                        </Typography>
                    </Option>
                ))}
            </Select>
        </Box>
    )
}

const Popup = styled(Box)`
  z-index: 1;
`;