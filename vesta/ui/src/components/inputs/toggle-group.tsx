import * as React from 'react';
import { styled } from '@mui/system';
import { Tabs as BaseTabs } from '@mui/base/Tabs';
import { TabsList as BaseTabsList } from '@mui/base/TabsList';
import { buttonClasses } from '@mui/base/Button';
import { Tab as BaseTab, tabClasses } from '@mui/base/Tab';
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';


export default function ToggleGroup({
  valueIdx,
  values,
  onChange,
}: {
  valueIdx: number,
  values: string[],
  onChange: (newValIdx: number) => void,
}) {
  const theme = useTheme()

  return (
    <Tabs
      value={valueIdx}
      // @ts-ignore
      onChange={(_, newIdx: number) => onChange(newIdx)}
    >
      <TabsList>
        {values.map(val => (
          <Tab key={val}>
            <Typography
              variant='pMediumMediumWt'
              sx={{
                color: theme.palette.ground.grade10,
              }}
            >
              {val}
            </Typography>
          </Tab>
        ))}
      </TabsList>
    </Tabs>
  );
}

const Tab = styled(BaseTab)(
  ({ theme }) => `
  display: flex;
  justify-content: center;
  align-items: center;
  text-wrap: nowrap;
  padding: 5px 12px;
  outline: none;
  border: none;
  border-radius: 40px;
  background: transparent;
  cursor: pointer;
  transition: ${theme.transitions.create('all', { easing: theme.transitions.easing.ease, duration: theme.transitions.duration.ease, })};

  &:hover {
    // background-color: black;
  }

  &:focus {
    // outline: 1px solid black;
  }

  &.${tabClasses.selected} {
    background: ${theme.palette.teal.grade100};
  }

  &.${buttonClasses.disabled} {
    opacity: 0.5;
    cursor: not-allowed;
  }
`
);

const TabsList = styled(BaseTabsList)(
  ({ theme }) => `
  display: flex;
  gap: 8px;
  `,
);

const Tabs = styled(BaseTabs)(
  ({ theme }) => `
  background-color: ${theme.palette.utilityWhite.main};
  border-radius: 40px;
  display: flex;
  align-items: center;
  justify-content: center;
  max-width: fit-content;
  padding: 8px;
  `
);