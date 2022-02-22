import _ from 'lodash';
import { withDefault } from '../../../../../selectors/maybe';
import { inputValueNonEmpty } from '../../../../../selectors/workflow_selectors';
import {DataEnvelope, ValidationInputSpecification} from '../input_types';
import AllInnerValuesNotEmptyValidator, { AllInnerValuesNotEmptyValidatorStrong } from './all_inner_values_not_empty_validator';

function PlusSubsetValidator(
  subset_key:string,
  OtherValidator: (input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>) => string[],
  ) {
    return (input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>): string[] => {
      
      const out = OtherValidator({...input, data: input.value as DataEnvelope<any>});
      
      let subset_check: string[] = []
      if (input.value!=null) {
        
        let if_some_null: DataEnvelope<any> = {}
        if_some_null[subset_key] = {}
        const subset = withDefault(input.value,if_some_null)[subset_key] as DataEnvelope<(string|number|null)[][]>
        
        if (!_.isEmpty(subset)) {
          let bad_conditions: string[] = []
          subset['methods'].forEach( (def, index) => {
            if (def.some(x => x==null) || def.length<2) {
              bad_conditions.push(''+(index+1))
            }
          })
          subset_check = (bad_conditions.length==0) ? [] : ['Subset condition(s) '+ bad_conditions.toString() +' are incomplete.']
        }
      }
      
      out.push(...subset_check)
      return out
    }
  }
  
export default PlusSubsetValidator
