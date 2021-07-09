// Input component that takes a simple object and
// shows values based on a selected key

import React, {
  useState,
  useEffect,
  useCallback,
  useContext,
  useMemo
} from 'react';

import {VulcanContext} from '../../../../contexts/vulcan_context';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import CheckboxesInput from './checkboxes';
import {InputBackendComponent, InputSpecification, WithInputParams} from './input_types';
import {flattenOptions} from './multiple_multiselect_string_all';
import {Maybe, some} from "../../../../selectors/maybe";

// function DropdownCheckboxCombo({
//   dropdownValue,
//   handleSelect,
//   handleClickOption,
//   dropdownOptions,
//   checkboxInput
// }: {
//   dropdownValue: string;
//   dropdownOptions: string[];
//   handleSelect: (value: string) => void;
//   handleClickOption: (options: string[]) => void;
//   checkboxInput: InputSpecification;
// }) {
//   if (!checkboxInput) return null;
//
//   return (
//     <React.Fragment>
//       <div>
//         <DropdownAutocomplete
//           value={dropdownValue}
//           onSelect={handleSelect}
//           list={dropdownOptions}
//         />
//       </div>
//       <div className='checkbox-input-wrapper'>
//         <CheckboxesInput
//           input={checkboxInput}
//           onChange={(inputName: string, checkedOptions: string[]) =>
//             handleClickOption(checkedOptions)
//           }
//         />
//       </div>
//     </React.Fragment>
//   );
// }
//
// export default function SingleDropdownMulticheckbox({ data, onChange, value }: WithInputParams<{}>) {
//   // input.data for this component should be
//   // {'a': {experiment: ['1', '2'], tissue: ['a', 'b']}}
//   const [checkboxInput, setCheckboxInput] = useState({} as InputSpecification);
//   const [dropdownOptions, setDropdownOptions] = useState([] as string[]);
//   const allOptions: {[k: string]: string[]} = useMemo(
//     () => flattenOptions(data),
//     [data]
//   );
//   const allValues: string[] = useMemo(() =>
//     Object.values(allOptions).reduce((a, b) => [...a, ...b], []), [allOptions])
//
//   const [dropdownValue, setDropdownValue] = useState(null as Maybe<string>);
//   // default with all options selected.
//   useEffect(() => {
//     const optionKeys = Object.keys(allOptions);
//     setDropdownOptions(optionKeys);
//
//     if (null == value) {
//       onChange(some(allValues));
//     }
//   }, [allOptions, onChange, value, allValues]);
//
//   useEffect(() => {
//     // when options change for the dropdown menu,
//     //   set the checkboxes to the options
//     //   for the selected dropdownOption.
//     setCheckboxInput({
//       ...input,
//       name: currentDropdownValue,
//       data: allOptions ? allOptions[currentDropdownValue] : [],
//       value:
//         input.value && input.value[currentDropdownValue]
//           ? input.value[currentDropdownValue]
//           : null
//     });
//   }, [currentDropdownValue, input.value, allOptions, input]);
//
//   const handleSelect = useCallback(
//     (value: string) => {
//       onChange(internalDropdownParameter, value);
//     },
//     [internalDropdownParameter, onChange]
//   );
//
//   const handleClickOption = useCallback(
//     (options: string[]) => {
//       onChange(input.name, {
//         ...input.value,
//         [currentDropdownValue]: [...options]
//       });
//     },
//     [input, onChange, currentDropdownValue]
//   );
//
//   if (!input || !onChange) return null;
//
//   return (
//     <div>
//       <DropdownCheckboxCombo
//         checkboxInput={checkboxInput}
//         handleSelect={handleSelect}
//         handleClickOption={handleClickOption}
//         dropdownValue={currentDropdownValue}
//         dropdownOptions={dropdownOptions}
//       />
//     </div>
//   );
// };
