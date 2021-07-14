// import React from 'react';
// import {mount, ReactWrapper} from 'enzyme';
// import MultipleMultiselectStringAllInput from '../multiple_multiselect_string_all';
// import {InputSpecification} from '../input_types';
//
// describe('MultipleMultiselectStringAllInput', () => {
//   let input: InputSpecification;
//   let onChange: jest.Mock;
//
//   function openDropdown(component: ReactWrapper, inputIndex: number) {
//     component
//       .find('.view_item')
//       .at(inputIndex)
//       .find('.add_item')
//       .first()
//       .simulate('click')
//       .update()
//       .find('.icon-wrapper')
//       .first()
//       .simulate('click')
//       .update();
//   }
//
//   function clickOption(
//     component: ReactWrapper,
//     inputIndex: number,
//     optionIndex: number
//   ) {
//     openDropdown(component, inputIndex);
//     component.find('li').at(optionIndex).simulate('click').update();
//   }
//
//   function renderedItemsList(component: ReactWrapper) {
//     return component.find('.delete_link').map((n) => n.text());
//   }
//
//   function renderedLabels(component: ReactWrapper) {
//     return component.find('.item_name').map((n) => n.text());
//   }
//
//   function renderedOptions(component: ReactWrapper) {
//     return component.find('li').map((n) => n.text());
//   }
//
//   beforeEach(() => {
//     input = {
//       type: 'doesnotmatter',
//       value: null,
//       label: 'Abcdef',
//       name: 'test-input',
//       data: {
//         'options-a': {
//           option1: ['1', '2', '3'],
//           option2: ['x', 'y', 'z']
//         },
//         'options-b': {
//           option3: ['9', '8', '7'],
//           option4: ['a', 'b', 'c']
//         }
//       }
//     };
//
//     onChange = jest.fn();
//   });
//
//   it('correctly manages child(ren) selects', () => {
//     const component = mount(
//       <MultipleMultiselectStringAllInput input={input} onChange={onChange} />
//     );
//     expect(component.find('.view_item').length).toEqual(4);
//     expect(renderedLabels(component)).toEqual([
//       'option1',
//       'option2',
//       'option3',
//       'option4'
//     ]);
//     clickOption(component, 0, 0);
//     expect(onChange).toHaveBeenLastCalledWith('test-input', {
//       option1: ['1']
//     });
//
//     clickOption(component, 1, 1);
//     expect(onChange).toHaveBeenLastCalledWith('test-input', {
//       option2: ['y']
//     });
//
//     clickOption(component, 2, 2);
//     expect(onChange).toHaveBeenLastCalledWith('test-input', {
//       option3: ['9']
//     });
//
//     clickOption(component, 3, 0);
//     expect(onChange).toHaveBeenLastCalledWith('test-input', {
//       option4: ['a']
//     });
//   });
//
//   it('correctly sets the lists when given a value', () => {
//     input.value = {
//       option1: ['1', '2'],
//       option2: ['y'],
//       option3: ['8', '9', '7'],
//       option4: ['c']
//     };
//     const component = mount(
//       <MultipleMultiselectStringAllInput input={input} onChange={onChange} />
//     );
//
//     expect(component.find('.view_item').length).toEqual(4);
//     expect(renderedLabels(component)).toEqual([
//       'option1',
//       'option2',
//       'option3',
//       'option4'
//     ]);
//     expect(renderedItemsList(component)).toEqual([
//       '1',
//       '2',
//       'y',
//       '8',
//       '9',
//       '7',
//       'c'
//     ]);
//   });
//
//   it('can remove a single entry with a value', () => {
//     input.value = {
//       option1: ['1', '2'],
//       option2: ['y'],
//       option3: ['8', '9', '7'],
//       option4: ['c']
//     };
//     const component = mount(
//       <MultipleMultiselectStringAllInput input={input} onChange={onChange} />
//     );
//
//     expect(component.find('.view_item').length).toEqual(4);
//     expect(renderedLabels(component)).toEqual([
//       'option1',
//       'option2',
//       'option3',
//       'option4'
//     ]);
//
//     expect(renderedItemsList(component)).toEqual([
//       '1',
//       '2',
//       'y',
//       '8',
//       '9',
//       '7',
//       'c'
//     ]);
//
//     component.find('.delete_link').last().simulate('click');
//     expect(onChange).toHaveBeenCalledWith('test-input', {
//       option1: ['1', '2'],
//       option2: ['y'],
//       option3: ['8', '9', '7'],
//       option4: []
//     });
//   });
// });
