import * as React from 'react';
import { styled } from '@mui/system';
import { Tabs as BaseTabs } from '@mui/base/Tabs';
import { TabsList as BaseTabsList } from '@mui/base/TabsList';
import { TabPanel as BaseTabPanel } from '@mui/base/TabPanel';
import { buttonClasses } from '@mui/base/Button';
import { Tab as BaseTab, tabClasses } from '@mui/base/Tab';


type Value = string | number | null


export default function Tabs({
    valueIdx,
    values,
    onChange,
}: {
    valueIdx: Value,
    values: string[],
    onChange: (newValIdx: Value) => void,
}) {
    return (
        <BaseTabs
            value={valueIdx}
            onChange={(_, newIdx) => onChange(newIdx)}
        >
            <TabsList>
                {values.map(val => (
                    <Tab key={val}>{val}</Tab>
                ))}
            </TabsList>
        </BaseTabs>
    );
}

const Tab = styled(BaseTab)`

  &:hover {
    // background-color: black;
  }

  &:focus {
    // color: #fff;
    // outline: 3px solid black;
  }

  &.${tabClasses.selected} {
    // background-color: #fff;
    // color: black;
  }

  &.${buttonClasses.disabled} {
    // opacity: 0.5;
    // cursor: not-allowed;
  }
`;

const TabPanel = styled(BaseTabPanel)`
  width: 100%;
//   font-family: 'IBM Plex Sans', sans-serif;
//   font-size: 0.875rem;
`;

const TabsList = styled(BaseTabsList)(
    ({ theme }) => `
//   min-width: 400px;
//   background-color: black;
//   border-radius: 12px;
//   margin-bottom: 16px;
//   display: flex;
//   align-items: center;
//   justify-content: center;
//   align-content: space-between;
//   box-shadow: 0px 4px 6px ${theme.palette.mode === 'dark' ? 'rgba(0,0,0, 0.4)' : 'rgba(0,0,0, 0.2)'
        };
  `,
);