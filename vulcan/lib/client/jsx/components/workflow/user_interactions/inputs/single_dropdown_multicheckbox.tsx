// // Input component that takes a simple object and
// // shows values based on a selected key

// import React, {useState, useEffect, useCallback} from 'react';
// import * as _ from 'lodash';

// import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';

// import {useBufferedInputState} from './buffered_input_state';
// import CheckboxesInput from './checkboxes';
// import {InputBackendComponent, InputSpecification} from './input_types';

// function DropdownCheckboxCombo({
//   dropdownValue,
//   handleSelect,
//   handleClickOption,
//   cwlInputName,
//   dropdownOptions,
//   checkboxInput
// }: {
//   dropdownValue: string;
//   dropdownOptions: string[];
//   cwlInputName: string;
//   handleSelect: (cwlInputName: string, value: string) => void;
//   handleClickOption: (cwlInputName: string, options: string[]) => void;
//   checkboxInput: InputSpecification;
// }) {
//   if (!checkboxInput) return null;

//   return (
//     <React.Fragment>
//       <div>
//         <DropdownAutocomplete
//           defaultValue={dropdownValue}
//           onSelect={(value: string) => handleSelect(cwlInputName, value)}
//           list={dropdownOptions}
//         />
//       </div>
//       <div className='checkbox-input-wrapper'>
//         <CheckboxesInput
//           input={checkboxInput}
//           onChange={(inputName: string, checkedOptions: string[]) =>
//             handleClickOption(cwlInputName, checkedOptions)
//           }
//         />
//       </div>
//     </React.Fragment>
//   );
// }

// const SingleDropdownMulticheckbox: InputBackendComponent = ({
//   input,
//   onChange
// }) => {
//   const [checkboxInputs, setCheckboxInputs] = useState(
//     {} as {[key: string]: InputSpecification | null}
//   );

//   const [selected, setSelected] = useBufferedInputState<{
//     [key: string]: {[key: string]: string[]};
//   }>(input, {});

//   const [currentDropdownValues, setCurrentDropdownValues] = useState(
//     {} as {[key: string]: string}
//   );

//   const [dropdownOptions, setDropdownOptions] = useState(
//     {} as {[key: string]: string[]}
//   );
//   const {data} = input;

//   useEffect(() => {
//     if (!data) return;

//     setDropdownOptions(
//       Object.keys(data).reduce(
//         (
//           acc: {[key: string]: string[]},
//           nextCwlInput: string
//         ): {[key: string]: string[]} => {
//           acc[nextCwlInput] = Object.keys(data[nextCwlInput]);
//           return acc;
//         },
//         {}
//       )
//     );
//   }, [data]);

//   useEffect(() => {
//     // By default set the first option as selected for each CWL input.
//     setCurrentDropdownValues(
//       Object.entries(dropdownOptions).reduce(
//         (
//           acc: {[key: string]: string},
//           [nextCwlInput, options]: [string, string[]]
//         ) => {
//           acc[nextCwlInput] = options[0];
//           return acc;
//         },
//         {}
//       )
//     );
//   }, [dropdownOptions]);

//   useEffect(() => {
//     if (!initialized && !_.isEqual(dropdownOptions, {})) {
//       setInitialized(true);
//       if (input.value) {
//         // Set the boxes to any previous selection
//         setSelected(input.value);
//       }
//     }
//   }, [input.value, input.data, dropdownOptions]);

//   useEffect(() => {
//     // when options change for the dropdown menu,
//     //   set the checkboxes to the options
//     //   for the selected dropdownOption.
//     setCheckboxInputs(
//       Object.entries(currentDropdownValues).reduce(
//         (
//           acc: {[key: string]: InputSpecification},
//           [cwlInputName, dropdownValue]: [string, string]
//         ): {[key: string]: InputSpecification} => {
//           acc[cwlInputName] = {
//             ...input,
//             name: `${cwlInputName}-${dropdownValue}`,
//             data: data ? data[cwlInputName][dropdownValue] : [],
//             value:
//               selected && selected[cwlInputName]
//                 ? selected[cwlInputName][dropdownValue]
//                 : []
//           };
//           return acc;
//         },
//         {}
//       )
//     );
//   }, [currentDropdownValues, selected]);

//   useEffect(() => {
//     // Update state whenever a new selection comes in.
//     if (!_.isEqual(lastSelected, selected)) {
//       setLastSelected(selected);
//       if (
//         Object.values(
//           selected
//         ).every((nestedInput: {[key: string]: string[]}) =>
//           Object.values(nestedInput).every(
//             (selection: string[]) => selection.length > 0
//           )
//         )
//       ) {
//         onChange(input.name, selected);
//       } else {
//         onChange(input.name, null);
//       }
//     }
//   }, [selected]);

//   const handleSelect = useCallback(
//     (cwlInputName: string, value: string) => {
//       // Set the checkbox input to `null` first, to clear
//       //   out the component.
//       setCheckboxInputs({
//         ...checkboxInputs,
//         [cwlInputName]: null
//       });
//       setCurrentDropdownValues({
//         ...currentDropdownValues,
//         [cwlInputName]: value
//       });
//     },
//     [currentDropdownValues]
//   );

//   const handleClickOption = useCallback(
//     (cwlInputName: string, options: string[]) => {
//       setSelected({
//         ...selected,
//         [cwlInputName]: {
//           ...selected[cwlInputName],
//           [currentDropdownValues[cwlInputName]]: [...options]
//         }
//       });
//     },
//     [currentDropdownValues]
//   );

//   if (!input || !onChange) return null;

//   return (
//     <div>
//       {Object.keys(currentDropdownValues).map((cwlInputName: string) => {
//         let checkboxInput = checkboxInputs[cwlInputName];

//         if (!checkboxInput) return null;

//         return (
//           <DropdownCheckboxCombo
//             key={`${cwlInputName}-${checkboxInput.name}`}
//             cwlInputName={cwlInputName}
//             checkboxInput={checkboxInput}
//             handleSelect={handleSelect}
//             handleClickOption={handleClickOption}
//             dropdownValue={currentDropdownValues[cwlInputName]}
//             dropdownOptions={dropdownOptions[cwlInputName]}
//           />
//         );
//       })}
//     </div>
//   );
// };

// export default SingleDropdownMulticheckbox;
