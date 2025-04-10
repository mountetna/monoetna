import _ from 'lodash';
import { withDefault } from '../../../../../selectors/maybe';
import {DataEnvelope, ValidationInputSpecification} from '../../input_types';
import { SelectionDefinition, SingleSelectionDefinition } from '../pieces/grouping_pieces';
import { NumericConstraint } from '../pieces/user_input_pieces';

function checkDef(def: SingleSelectionDefinition, index: number) {
  const out = [] as string[];
  if (!('col' in def) || (def['col']==null)) {
    out.push(`SelectionDefinition Constraint ${index+1} has no 'col' chosen`)
    return out;
  }
  if (!('def' in def) || def['def']==null || def['def'].length < 1) {
    out.push(`SelectionDefinition Constraint ${index+1} value selection is not filled`)
  }
  if ('def' in def && def['def']!=null && def['def'].length==4 && typeof def['def'][1]==='number') {
    const _def = def['def'] as NumericConstraint
    const left = _def[1];
    const right = _def[3];
    if (left==null || right==null) {
      out.push(`SelectionDefinition Constraint ${index+1} numeric value selection is not filled`)
    } else if (left>right || (left==right && (_def[0]!='exactly' || _def[2]!='exactly'))) {
      out.push(`SelectionDefinition Constraint ${index+1} numeric value selection leaves no valid numbers`)
    }
  }
  if (!('logic' in def) || (index > 0 && (!def['logic'] || !["AND", "OR"].includes(def['logic'])))) {
    out.push(`SelectionDefinition Constraint ${index+1} combination logic not established`)
  }
  return out;
}

function PlusSubsetValidator(
  subset_key: string,
  OtherValidator: (input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>) => string[],
  ) {
    return (input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>): string[] => {
      
      const out = OtherValidator({...input, data: input.value as DataEnvelope<any>});

      let subset_check: string[] = [];
      if (input.value!=null) {

        const if_blank: DataEnvelope<any> = {[subset_key]: false};
        const overall = withDefault(input.value,if_blank)
        if (subset_key in overall) {
          const subset = overall[subset_key] as SelectionDefinition | false;

          if (subset===false) return out;

          for (let ind = 0; ind < subset.length; ind++) {
            subset_check = subset_check.concat(checkDef(subset[ind], ind))
          }
        }
      }
      
      out.push(...subset_check);
      return out;
    };
  }
  
export default PlusSubsetValidator;
