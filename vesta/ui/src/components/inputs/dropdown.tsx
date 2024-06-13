import * as React from 'react';
import Box from '@mui/system/Box'
import { Select } from '@mui/base/Select';
import { Option } from '@mui/base/Option';
import { Button } from '@mui/base/Button';
import _ from 'lodash'
import Image from 'next/image';

import indicatorArrowDark from '/public/images/icons/indicator-arrow-dark.svg'


interface OptionType {
    value: any
    label: string
}

const IndicatorArrowDark = (
    <Box
        sx={{
            display: 'inline-block',
            transform: 'rotate(180deg)',
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
}: {
    options: OptionType[],
    value: any,
    onChange: (newValue: any) => void,
    icon?: React.ReactNode,
}) {
    const listboxId = _.uniqueId('dropdown-')

    return (
        <Box
            className='dropdown'
            sx={{
                bgcolor: 'utilityHighlight.main',
                borderRadius: '10px',
            }}
        >
            <Select
                value={value}
                onChange={(_, newValue) => onChange(newValue)}
                style={{
                    width: '100%',
                    border: 'none',
                    borderRadius: '10px',
                }}
                listboxId={listboxId}
                renderValue={(option) => (
                    <Box
                        className='dropdown-button-content'
                        sx={{
                            display: 'flex',
                            justifyContent: 'space-between',
                            alignItems: 'center',
                        }}
                    >
                        <Box component='span'>{option?.label}</Box>
                        {icon}
                    </Box>
                )}
            >
                {options.map(option => (
                    <Option
                        key={option.value}
                        value={option.value}
                        style={{

                        }}
                    >
                        {option.label}
                    </Option>
                ))}
            </Select>
        </Box>
    )
}