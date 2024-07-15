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
import { styled } from '@mui/system';
import useForkRef from '@mui/utils/useForkRef';
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import { StaticImageData } from 'next/image';
import { forwardRef } from '@/lib/utils/types';


interface _AutocompleteGroupedOption<Value = string> extends Omit<AutocompleteGroupedOption<Value>, 'options'> {
  options: React.ReactNode;
}


interface AutocompleteProps<Value> extends UseAutocompleteProps<Value, boolean, boolean, boolean> {
  icon: StaticImageData;
  placeholder: string;
  renderOption: (params: UseAutocompleteRenderedOption<Value>) => React.ReactNode;
  renderGroup?: (params: _AutocompleteGroupedOption<Value>) => React.ReactNode;
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
      return <NoOptions>No results</NoOptions>
    }

    if (props.groupBy === undefined) {
      return renderOptions(
        groupedOptions as Value[],
        getOptionProps,
        props.renderOption
      )
    } else {
      return (groupedOptions as AutocompleteGroupedOption<Value>[]).map((group) => {

        const options = renderOptions(
          group.options,
          getOptionProps,
          props.renderOption
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

  return (
    <React.Fragment>
      <Root
        {...getRootProps()}
        ref={rootRef}
        className={focused ? 'Mui-focused' : ''}
      >
        <StyledInput {...getInputProps()} />
      </Root>
      {anchorEl && (
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
          <Listbox {...getListboxProps()}>
            {renderGroupedOptions()}
          </Listbox>
        </Popper>
      )}
    </React.Fragment >
  );
};

export default forwardRef(Autocomplete)


function renderOptions<Value>(
  options: Value[],
  getOptionProps: UseAutocompleteReturnValue<Value>['getOptionProps'],
  renderOption: AutocompleteProps<Value>['renderOption'],
): React.ReactNode {



  return options.map((option, index) => {
    const optionProps = getOptionProps({ option, index })
    // @ts-ignore
    const key = optionProps.key as string
    // @ts-ignore
    delete optionProps.key
    // }

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


// export function UseAutocompletePopper() {
//   const [value, setValue] = React.useState<(typeof top100Films)[number] | null>(
//     null,
//   );

//   const handleChange = (
//     event: React.SyntheticEvent,
//     newValue: (typeof top100Films)[number] | null,
//   ) => setValue(newValue);

//   return (
//     <Autocomplete options={top100Films} value={value} onChange={handleChange} />
//   );
// }

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
    };
  `,
);

const Option = styled('li')(
  ({ theme }) => `
  list-style: none;
  cursor: default;

  &:last-of-type {
  }

  &:hover {
    cursor: pointer;
  }

  &[aria-selected=true] {
  }

  &.Mui-focused,
  &.Mui-focusVisible {
  }

  &.Mui-focusVisible {
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