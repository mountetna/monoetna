import {InputType} from 'etna-js/utils/input_types';
export {InputType} from 'etna-js/utils/input_types';
import BooleanInput from "./workspace/ui_definitions/inputs/components/boolean";
import CheckboxesInput from "./workspace/ui_definitions/inputs/components/checkboxes";
import { AnnotationEditorInput, DataTransformationInput } from "./workspace/ui_definitions/inputs/data_transformation";
import FloatInput from "./workspace/ui_definitions/inputs/components/float";
import { InputBackendComponent, InputValidator } from "./workspace/ui_definitions/input_types";
import IntegerInput from "./workspace/ui_definitions/inputs/components/integer";
import { MetisFileInput, MetisFolderInput } from "./workspace/ui_definitions/inputs/components/metis_items";
import MultipleInput from "./workspace/ui_definitions/inputs/multiple_input";
import MultiselectStringInput from "./workspace/ui_definitions/inputs/components/multiselect_string";
import NestedDropdownInput from "./workspace/ui_definitions/inputs/components/nested_dropdown";
// import NestedDropdownMultiChoicePieceRct from "./workspace/ui_definitions/inputs/pieces/nested_dropdown_multi_choice_piece";
import DiffExpSC from "./workspace/ui_definitions/inputs/scDGE";
import DropdownInput from "./workspace/ui_definitions/inputs/components/dropdown";
import DropdownMultiChoiceInput from "./workspace/ui_definitions/inputs/components/dropdown_multi_choice";
import SingleDropdownMulticheckbox from "./workspace/ui_definitions/inputs/single_dropdown_multicheckbox";
import StringInput from "./workspace/ui_definitions/inputs/components/string";
import { AnyDittoSeq, AnyPlotly, BarPlotly, DittoBarPlot, DittoDimPlot, DittoPlot, DittoScatterPlot, ScatterPlotly, YPlotly } from "./workspace/ui_definitions/inputs/components/visualizations";

import AllInnerKeysNotNullValidator from "./workspace/ui_definitions/inputs/validators/all_inner_keys_not_null_validator";
import AllInnerValuesNotEmptyValidator from "./workspace/ui_definitions/inputs/validators/all_inner_values_not_empty_validator";
import AllOutputValuesNotEmptyValidator, { AllOutputValuesNotEmptyAllowingEmptyArrayValidator } from "./workspace/ui_definitions/inputs/validators/all_output_values_not_empty_validator";
import AnnotationEditorValidator from "./workspace/ui_definitions/inputs/validators/annotation_editor_validator";
import { MetisFileValidator, MetisFolderValidator, MetisPathValidator } from "./workspace/ui_definitions/inputs/validators/metis_path_validators";
import { NotEmptyValidator, StronglyNotEmptyValidator } from "./workspace/ui_definitions/inputs/validators/not_empty_validator";
import PlusSubsetValidator from "./workspace/ui_definitions/inputs/validators/PlusSubsetValidator";

import LinkOutput from './workspace/ui_definitions/outputs/link';
import { PlotlyOutput, PlotOutput, PngOutput } from './workspace/ui_definitions/outputs/plot';
import ConsignmentOutput from './workspace/ui_definitions/outputs/consignment';
import RawOutput from './workspace/ui_definitions/outputs/raw';
import TwoGroupSelection from './workspace/ui_definitions/inputs/components/two_group_selection';
import { NestedDropdownMultiChoiceInput, NestedDropdownMultiChoiceAdvancedInput, NestedDropdownMultiChoiceBulkAddInput, NestedDropdownMultiChoiceReorderInput } from './workspace/ui_definitions/inputs/components/nested_dropdown_multi_choice';
import { MagmaRecordInput } from './workspace/ui_definitions/inputs/components/magma_record';

/*
InputComponents: 
  - For Input_Feed
  - I/O:
    - Some require inputs
    - Can produce outputs
  - Tracked as a step by snakemake
Defined here via:
- token - string utilized to request this input setup in a vulcan_config.yaml
- component - jsx defining the UI itself (more below)
- validator - jsx defining any when to block the user from 'Confirm'ing their selections, as well as the error messages to display (more below)
- outputKeys - string[] naming the 'value' keys managed by the 'component' jsx
- inputKeysRequired - string[] naming the 'data' keys expect by the 'component' jsx
- inputKeyOptional - string naming a single 'data' key only optionally expected by the 'component' jsx

'inputComponents' export holds all of these definitions.

components:
  - take as input:
    - current 'value', as Object where all values map to potential outputs
    - a standard 'onChange' function, not defined here, called whenever updating the current value
    - input 'dataElement'
    - (additional 'params', managed externally)
  - output:
    - A React Component

validators:
  - take as input:
    - current 'value'
    - input 'dataElement'
  - output:
    - string[] that is empty when there are no 'value' validation errors
*/

export const inputComponents = {} as {[k: string]: [
  InputBackendComponent<any, any, any>,
  InputValidator<any, any>,
  string[],
  string[],
  string | undefined
]};
function configureComponent<Value, DataElement>(
  type: InputType,
  component: InputBackendComponent<any, Value, DataElement>,
  validator: InputValidator<Value, DataElement>,
  outputKeys: string[],
  inputKeysRequired: string[],
  inputKeysOptional?: string
) {
  if (type in inputComponents) throw new Error(`Duplicate definition for ${type}`);
  inputComponents[type] = [component, validator, outputKeys, inputKeysRequired, inputKeysOptional];
}
configureComponent('default', StringInput, NotEmptyValidator, ['value'], []);
configureComponent('string', StringInput, NotEmptyValidator, ['value'], []);
configureComponent('float', FloatInput, NotEmptyValidator, ['value'], []);
configureComponent('int', IntegerInput, NotEmptyValidator, ['value'], []);
configureComponent('boolean', BooleanInput, NotEmptyValidator, ['value'], [], 'label');
configureComponent('checkbox', BooleanInput, NotEmptyValidator, ['value'], [], 'label');
configureComponent('checkboxes', CheckboxesInput, NotEmptyValidator, ['picked'], ['options']);
configureComponent('magma-record', MagmaRecordInput, NotEmptyValidator, ['record'], ['modelName', 'targetAtt', 'otherAttsShow'], 'hasAtts')
configureComponent('two-group-selection', TwoGroupSelection, AllOutputValuesNotEmptyValidator, ['g1', 'g2'], ['data_summary'], 'all_column_options')
configureComponent('multiselect-string', MultiselectStringInput, NotEmptyValidator, ['picked'], ['options']);
configureComponent('MetisFile', MetisFileInput, MetisFileValidator(), ['target'], []);
configureComponent('MetisCSVorTSV', MetisFileInput, MetisFileValidator('\\.[ct]sv$', 'csv or tsv file'), ['target'], []);
configureComponent('MetisFolder', MetisFolderInput, MetisFolderValidator(), ['target'], []);
configureComponent('MetisPath', MetisFileInput, MetisPathValidator(), ['target'], []);
configureComponent('dropdown', DropdownInput, StronglyNotEmptyValidator, ['picked'], ['options'], 'recommendation');
configureComponent('nested-dropdown', NestedDropdownInput, StronglyNotEmptyValidator, ['picked'], ['nestedOptions'], 'recommendation');
configureComponent('dropdown-multi-choice', DropdownMultiChoiceInput, StronglyNotEmptyValidator, ['picked'], ['options'], 'recommendation');
configureComponent('nested-dropdown-multi-choice', NestedDropdownMultiChoiceInput, StronglyNotEmptyValidator, ['picked'], ['nestedOptions'], 'recommendation');
configureComponent('nested-dropdown-multi-choice-bulk-add', NestedDropdownMultiChoiceBulkAddInput, StronglyNotEmptyValidator, ['picked'], ['nestedOptions'], 'recommendation');
configureComponent('nested-dropdown-multi-choice-reorder', NestedDropdownMultiChoiceReorderInput, StronglyNotEmptyValidator, ['picked'], ['nestedOptions'], 'recommendation');
configureComponent('nested-dropdown-multi-choice-advanced', NestedDropdownMultiChoiceAdvancedInput, StronglyNotEmptyValidator, ['picked'], ['nestedOptions'], 'recommendation');
// python plotly visualization
configureComponent('any-viz', AnyPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent('scatter-plotly', ScatterPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent('bar-plotly', BarPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent('y-plotly', YPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator), ['plot_setup'], ['data_frame'], 'optional');
// dittoSeq visualization
configureComponent('any-dittoseq', AnyDittoSeq, PlusSubsetValidator('cells_use',AllOutputValuesNotEmptyAllowingEmptyArrayValidator), ['plot_setup'], ['data_frame', 'continuous_opts', 'discrete_opts', 'all_opts', 'reduction_opts'], 'optional');
configureComponent('dittoseq-dim-plot', DittoDimPlot, PlusSubsetValidator('cells_use',AllOutputValuesNotEmptyAllowingEmptyArrayValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent('dittoseq-scatter-plot', DittoScatterPlot, PlusSubsetValidator('cells_use',AllOutputValuesNotEmptyAllowingEmptyArrayValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent('dittoseq-plot', DittoPlot, PlusSubsetValidator('cells_use',AllOutputValuesNotEmptyAllowingEmptyArrayValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent('dittoseq-bar-plot', DittoBarPlot, PlusSubsetValidator('cells_use',AllOutputValuesNotEmptyAllowingEmptyArrayValidator), ['plot_setup'], ['data_frame'], 'optional');
// NOT YET UPDATED
// configureComponent('multiple-string', MultipleInput(StringInput), AllInnerValuesNotEmptyValidator);
// configureComponent('data-transformation', DataTransformationInput, AllInnerKeysNotNullValidator);
// configureComponent('annotation-editor', AnnotationEditorInput, AnnotationEditorValidator);
// configureComponent('diff-exp-sc', DiffExpSC, PlusSubsetValidator('subset',AllOutputValuesNotEmptyAllowingEmptyArrayValidator));
// configureComponent('single-dropdown-multicheckbox', SingleDropdownMulticheckbox, NotEmptyValidator);
// configureComponent('multiple-multiselect-string-all', MultipleInput(MultiselectStringInput), AllInnerValuesNotEmptyValidator);


/*
Outputs: For Output_Feed, cannot produce an output, NOT tracked as a step by snakemake
How they work:
- Defined through a jsx definition of the component only.
- value in OUTPUT_TYPES object is the string that can be used for ui_component in a vulcan_config.
- Additionally, if the input data is very large, addinf the type to 'dontDownloadForOutputTypes' is a VERY good idea.
*/
export const OUTPUT_TYPES = {
  LINK: 'link',
  PLOT: 'plot',
  PLOTLY: 'plotly',
  PNG: 'png',
  CONSIGNMENT: 'consignment',
  RAW: 'raw'
};

export const dontDownloadForOutputTypes = ['default', OUTPUT_TYPES.LINK];
export const dataAndUrlForOutputTypes = [] as string[];

export const OUTPUTS = {
  default: LinkOutput,
  [OUTPUT_TYPES.LINK]: LinkOutput,
  [OUTPUT_TYPES.PLOTLY]: PlotlyOutput,
  [OUTPUT_TYPES.PLOT]: PlotOutput,
  [OUTPUT_TYPES.PNG]: PngOutput,
  [OUTPUT_TYPES.CONSIGNMENT]: ConsignmentOutput,
  [OUTPUT_TYPES.RAW]: RawOutput
};