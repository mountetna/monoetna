'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import {
  useAutocomplete,
  UseAutocompleteProps,
  UseAutocompleteRenderedOption,
  UseAutocompleteReturnValue,
  AutocompleteGroupedOption,
} from '@mui/base/useAutocomplete';
import { Popper } from '@mui/base/Popper';
import { styled } from '@mui/material';
import useForkRef from '@mui/utils/useForkRef';
import { useTheme, SxProps } from '@mui/material';
import Image, { StaticImageData } from 'next/image';
import _ from 'lodash'

import { forwardRef } from '@/lib/utils/types';


interface _AutocompleteGroupedOption<Value = string> extends Omit<AutocompleteGroupedOption<Value>, 'options'> {
  options: React.ReactNode;
}


interface AutocompleteProps<Value> extends UseAutocompleteProps<Value, boolean, boolean, boolean> {
  icon: StaticImageData;
  placeholder?: string;
  showResultsOnEmpty?: boolean;
  renderOption: (params: UseAutocompleteRenderedOption<Value>) => React.ReactNode;
  renderGroup?: (params: _AutocompleteGroupedOption<Value>) => React.ReactNode;
  renderNoResults?: () => React.ReactNode;
  sx?: SxProps;
}

// Adapted from https://mui.com/base-ui/react-autocomplete/#using-a-portal
function Autocomplete<Value>(
  props: AutocompleteProps<Value>,
  ref: React.ForwardedRef<HTMLDivElement>,
) {
  const {
    getRootProps,
    getInputProps,
    getListboxProps,
    getOptionProps,
    groupedOptions,
    focused,
    popupOpen,
    anchorEl,
    setAnchorEl,
  } = useAutocomplete(props);

  const rootRef = useForkRef(ref, setAnchorEl);

  const renderGroupedOptions = () => {
    if (groupedOptions.length === 0) {
      return (
        <NoOptions>
          {props.renderNoResults !== undefined ?
            props.renderNoResults() :
            'No results'
          }
        </NoOptions>
      )
    }

    if (props.groupBy === undefined) {
      return renderOptions(
        groupedOptions as Value[],
        getOptionProps,
        props.renderOption,
      )
    } else {
      return (groupedOptions as AutocompleteGroupedOption<Value>[]).map((group) => {

        const options = renderOptions(
          group.options,
          getOptionProps,
          props.renderOption,
          group.index,
        )
        const renderGroupParams: _AutocompleteGroupedOption<Value> = {
          ...group,
          options,
        }

        return props.renderGroup === undefined ?
          defaultRenderGroup(renderGroupParams) :
          props.renderGroup(renderGroupParams)
      })
    }
  }

  const inputValue = getInputProps().value?.valueOf()
  const inputHasValue = inputValue !== undefined && inputValue.toString().length > 0
  const showResultsOnEmpty = props.showResultsOnEmpty === undefined ? true : props.showResultsOnEmpty
  const showResults = inputHasValue || showResultsOnEmpty

  const theme = useTheme()

  return (
    <Box
      sx={{
        display: 'flex',
        overflow: 'hidden',
        ...(props.sx || {})
      }}
    >
      <Root
        {...getRootProps()}
        ref={rootRef}
        sx={{
          display: 'flex',
          width: '100%',
          overflow: 'hidden',
        }}
      >
        <Box
          className={focused ? 'Mui-focused' : ''}
          sx={{
            width: '100%',
            height: '100%',
            overflowX: 'scroll',
            overflowY: 'clip',
            display: 'flex',
            gap: '10px',
            alignItems: 'center',
            bgcolor: 'utilityWhite.main',
            borderRadius: '30px',
            p: '14px 16px',
            border: '1px solid transparent',
            '&:focus-visible, &.Mui-focused': {
              border: `1px solid ${theme.palette.ground.grade75}`,
            },
            '&:hover, &:focus-visible, &.Mui-focused': {
              bgcolor: 'ground.grade100',
            },
            transition: theme.transitions.create(
              'all',
              {
                easing: theme.transitions.easing.ease,
                duration: theme.transitions.duration.ease,
              }
            ),
          }}
        >
          <Image
            src={props.icon}
            alt='Search icon'
          />

          <StyledInput
            {...getInputProps()}
            placeholder={props.placeholder}
            sx={{
              // display: 'flex',
              // flexGrow: 1,
              border: 'none',
              outline: 'none',
              p: '0',
              background: 'inherit',
              ...theme.typography.pBody,
              color: 'ground.grade10',
              '&::placeholder': {
                color: '#777777'
              },
            }}
          />
        </Box>
      </Root>
      {showResults && anchorEl && (
        <Popper
          open={popupOpen}
          anchorEl={anchorEl}
          slots={{
            root: StyledPopper,
          }}
          style={{
            width: anchorEl.offsetWidth,
          }}
        >
          <Box
            sx={{
              bgcolor: 'utilityWhite.main',
              borderRadius: '16px',
              m: '10px 0',
              p: '16px',
              border: `1px solid ${theme.palette.ground.grade75}`
            }}
          >
            <Listbox
              {...getListboxProps()}
              sx={{
                listStyle: 'none',
                maxHeight: '40rem',
                overflow: 'scroll',
                m: '0',
                p: '0',
              }}
            >
              {renderGroupedOptions()}
            </Listbox>
          </Box>
        </Popper>
      )}
    </Box >
  );
};

export default forwardRef(Autocomplete)


function renderOptions<Value>(
  options: Value[],
  getOptionProps: UseAutocompleteReturnValue<Value>['getOptionProps'],
  renderOption: AutocompleteProps<Value>['renderOption'],
  // required as per https://github.com/mui/material-ui/blob/next/packages/mui-material/src/Autocomplete/Autocomplete.js#L679C17-L679C33
  // if using groups
  indexOffset: number = 0,
): React.ReactNode {



  return options.map((option, index) => {
    const optionProps = getOptionProps({ option, index: indexOffset + index })
    // @ts-ignore
    const key = optionProps.key as string
    // @ts-ignore
    delete optionProps.key

    return (
      <Option key={key} {...optionProps}>
        {renderOption({ option, index })}
      </Option>
    )
  })
}

function defaultRenderGroup<Value>(params: _AutocompleteGroupedOption<Value>) {
  return (
    <li key={params.key}>
      <h1>{params.group}</h1>
      <ul>{params.options}</ul>
    </li>
  )
}

const Root = styled('div')(
  ({ theme }) => `

  &.Mui-focused {
  }

  &:hover {
  }

  &:focus-visible {
  }
`,
);

const StyledInput = styled('input')(
  ({ theme }) => `
`,
);

// ComponentPageTabs has z-index: 1000
const StyledPopper = styled('div')`
  position: relative;
  z-index: 1001;
`;

const Listbox = styled('ul')(
  ({ theme }) => `
  z-index: 1;
  
  `,
);

const Option = styled('li')(
  ({ theme }) => `
  list-style: none;
  cursor: default;
  outline: 1px solid transparent;
  padding: 8px;
  border-radius: 8px;
  transition: ${theme.transitions.create('all', { easing: theme.transitions.easing.ease, duration: theme.transitions.duration.ease })};

  &:last-of-type {
  }

  &:hover {
    cursor: pointer;
    background-color: ${theme.palette.utilityHighlight.main};
  }

  &[aria-selected=true] {
  }

  &.Mui-focused,
  &.Mui-focusVisible {
    background-color: ${theme.palette.utilityHighlight.main};
  }

  &.Mui-focusVisible {
    background-color: ${theme.palette.utilityHighlight.main};
  }

  &[aria-selected=true].Mui-focused,
  &[aria-selected=true].Mui-focusVisible {
  }
  `,
);

const NoOptions = styled('li')`
  list-style: none;
  cursor: default;
`;