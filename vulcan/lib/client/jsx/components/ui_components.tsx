import {InputType} from 'etna-js/utils/input_types';
export {InputType, DataEnvelope} from 'etna-js/utils/input_types';
import BooleanInput from "./workflow/user_interactions/inputs/components/boolean";
import CheckboxesInput from "./workflow/user_interactions/inputs/components/checkboxes";
import { AnnotationEditorInput, DataTransformationInput } from "./workflow/user_interactions/inputs/data_transformation";
import FloatInput from "./workflow/user_interactions/inputs/components/float";
import { InputBackendComponent, InputValidator } from "./workflow/user_interactions/inputs/input_types";
import IntegerInput from "./workflow/user_interactions/inputs/components/integer";
import { MetisFileInput, MetisFolderInput } from "./workflow/user_interactions/inputs/components/metis_items";
import MultipleInput from "./workflow/user_interactions/inputs/multiple_input";
import MultiselectStringInput from "./workflow/user_interactions/inputs/components/multiselect_string";
import NestedSelectAutocompleteInput from "./workflow/user_interactions/inputs/nested_select_autocomplete";
import NestedSelectAutocompleteMultiPickInput from "./workflow/user_interactions/inputs/nested_select_autocomplete_multi_choice";
import DiffExpSC from "./workflow/user_interactions/inputs/scDGE";
import SelectAutocompleteInput from "./workflow/user_interactions/inputs/components/select_autocomplete";
import SelectAutocompleteMultiPickInput from "./workflow/user_interactions/inputs/select_autocomplete_multi_choice";
import SingleDropdownMulticheckbox from "./workflow/user_interactions/inputs/single_dropdown_multicheckbox";
import StringInput from "./workflow/user_interactions/inputs/components/string";
import { AnyDittoSeq, AnyPlotly, BarPlotly, DittoBarPlot, DittoDimPlot, DittoPlot, DittoScatterPlot, ScatterPlotly, YPlotly } from "./workflow/user_interactions/inputs/components/visualizations";

import AllInnerKeysNotNullValidator from "./workflow/user_interactions/inputs/validators/all_inner_keys_not_null_validator";
import AllInnerValuesNotEmptyValidator from "./workflow/user_interactions/inputs/validators/all_inner_values_not_empty_validator";
import AllOutputValuesNotEmptyValidator, { AllOutputValuesNotEmptyAllowingEmptyArrayValidator } from "./workflow/user_interactions/inputs/validators/all_output_values_not_empty_validator";
import AnnotationEditorValidator from "./workflow/user_interactions/inputs/validators/annotation_editor_validator";
import { MetisFileValidator, MetisFolderValidator, MetisPathValidator } from "./workflow/user_interactions/inputs/validators/metis_path_validators";
import { NotEmptyValidator, StronglyNotEmptyValidator } from "./workflow/user_interactions/inputs/validators/not_empty_validator";
import PlusSubsetValidator from "./workflow/user_interactions/inputs/validators/PlusSubsetValidator";

import LinkOutput from './workflow/user_interactions/outputs/link';
import { PlotlyOutput, PlotOutput, PngOutput } from './workflow/user_interactions/outputs/plot';
import ConsignmentOutput from './workflow/user_interactions/outputs/consignment';
import RawOutput from './workflow/user_interactions/outputs/raw';

/*
Inputs: For Input_Feed, can produce an output, tracked as a step by snakemake
How they work:
- Defined through both a jsx definition of the component && as well as validator function via the `configureComponent` at the bottom of the section.
- value in INPUT_TYPES object is the string that can be used for ui_component in a vulcan_config.
*/
export const INPUT_TYPES = {
  INTEGER: 'int',
  FLOAT: 'float',
  BOOL: 'boolean',
  STRING: 'string',
  ARRAY: 'array',
  FILE: 'File',
  METIS_FILE: 'MetisFile',
  METIS_CSV_OR_TSV: 'MetisCSVorTSV',
  METIS_FOLDER: 'MetisFolder',
  METIS_FILE_OR_FOLDER: 'MetisPath',
  MULTISELECT_STRING: 'multiselect-string',
  SELECT_AUTOCOMPLETE: 'select-autocomplete',
  SELECT_AUTOCOMPLETE_MULTI_PICK: 'select-autocomplete-multi-pick',
  CHECKBOXES: 'checkboxes',
  NESTED_SELECT_AUTOCOMPLETE: 'nested-select-autocomplete',
  NESTED_SELECT_AUTOCOMPLETE_MULTI_PICK: 'nested-select-autocomplete-multi-pick',
  MULTIPLE_MULTISELECT_STRING_ALL: 'multiple-multiselect-string-all',
  MULTIPLE_STRING: 'multiple-string',
  SINGLE_DROPDOWN_MULTICHECKBOX: 'single-dropdown-multicheckbox',
  SCATTER_PLOTLY: 'scatter-plotly',
  BAR_PLOTLY: 'bar-plotly',
  Y_PLOTLY: 'y-plotly',
  DITTOSEQ_DIM_PLOT: 'dittoseq-dim-plot',
  DITTOSEQ_SCATTER_PLOT: 'dittoseq-scatter-plot',
  DITTOSEQ_BAR_PLOT: 'dittoseq-bar-plot',
  DITTOSEQ_PLOT: 'dittoseq-plot',
  ANY_DITTOSEQ: 'any-dittoseq',
  DIFF_EXP_SC: 'diff-exp-sc',
  DATA_TRANSFORMATION: 'data-transformation',
  ANNOTATION_EDITOR: 'annotation-editor',
  ANY_VIZ: 'any-viz'
};

export const components = {} as {[k: string]: [
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
  outputInternalNames: string[],
  inputNames: string[],
  optionalInputName?: string
) {
  if (type in components) throw new Error(`Duplicate definition for ${type}`);
  components[type] = [component, validator, outputInternalNames, inputNames, optionalInputName];
}
configureComponent(INPUT_TYPES.STRING, StringInput, NotEmptyValidator, ['value'], []);
configureComponent(INPUT_TYPES.FLOAT, FloatInput, NotEmptyValidator, ['value'], []);
configureComponent(INPUT_TYPES.INTEGER, IntegerInput, NotEmptyValidator, ['value'], []);
configureComponent(INPUT_TYPES.BOOL, BooleanInput, NotEmptyValidator, ['value'], []);
configureComponent(INPUT_TYPES.CHECKBOXES, CheckboxesInput, NotEmptyValidator, ['picked'], ['options']);
configureComponent(INPUT_TYPES.MULTISELECT_STRING, MultiselectStringInput, NotEmptyValidator, ['picked'], ['options']);
configureComponent(INPUT_TYPES.METIS_FILE, MetisFileInput, MetisFileValidator(), ['target'], []);
configureComponent(INPUT_TYPES.METIS_CSV_OR_TSV, MetisFileInput, MetisFileValidator('\\.[ct]sv$', 'csv or tsv file'), ['target'], []);
configureComponent(INPUT_TYPES.METIS_FOLDER, MetisFolderInput, MetisFolderValidator(), ['target'], []);
configureComponent(INPUT_TYPES.METIS_FILE_OR_FOLDER, MetisFileInput, MetisPathValidator(), ['target'], []);
configureComponent(INPUT_TYPES.SELECT_AUTOCOMPLETE, SelectAutocompleteInput, StronglyNotEmptyValidator, ['picked'], ['options'], 'recommendation');
configureComponent(INPUT_TYPES.SELECT_AUTOCOMPLETE_MULTI_PICK, SelectAutocompleteMultiPickInput, StronglyNotEmptyValidator);
configureComponent(INPUT_TYPES.NESTED_SELECT_AUTOCOMPLETE, NestedSelectAutocompleteInput, StronglyNotEmptyValidator);
configureComponent(INPUT_TYPES.NESTED_SELECT_AUTOCOMPLETE_MULTI_PICK, NestedSelectAutocompleteMultiPickInput, StronglyNotEmptyValidator);
configureComponent(INPUT_TYPES.MULTIPLE_STRING, MultipleInput(StringInput), AllInnerValuesNotEmptyValidator);
configureComponent(INPUT_TYPES.DATA_TRANSFORMATION, DataTransformationInput, AllInnerKeysNotNullValidator);
configureComponent(INPUT_TYPES.ANNOTATION_EDITOR, AnnotationEditorInput, AnnotationEditorValidator);
configureComponent(INPUT_TYPES.SCATTER_PLOTLY, ScatterPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent(INPUT_TYPES.BAR_PLOTLY, BarPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent(INPUT_TYPES.Y_PLOTLY, YPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent(INPUT_TYPES.ANY_VIZ, AnyPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent(INPUT_TYPES.DITTOSEQ_DIM_PLOT, DittoDimPlot, PlusSubsetValidator('cells_use',AllOutputValuesNotEmptyAllowingEmptyArrayValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent(INPUT_TYPES.DITTOSEQ_SCATTER_PLOT, DittoScatterPlot, PlusSubsetValidator('cells_use',AllOutputValuesNotEmptyAllowingEmptyArrayValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent(INPUT_TYPES.DITTOSEQ_PLOT, DittoPlot, PlusSubsetValidator('cells_use',AllOutputValuesNotEmptyAllowingEmptyArrayValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent(INPUT_TYPES.DITTOSEQ_BAR_PLOT, DittoBarPlot, PlusSubsetValidator('cells_use',AllOutputValuesNotEmptyAllowingEmptyArrayValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent(INPUT_TYPES.ANY_DITTOSEQ, AnyDittoSeq, PlusSubsetValidator('cells_use',AllOutputValuesNotEmptyAllowingEmptyArrayValidator), ['plot_setup'], ['data_frame'], 'optional');
configureComponent(INPUT_TYPES.DIFF_EXP_SC, DiffExpSC, PlusSubsetValidator('subset',AllOutputValuesNotEmptyAllowingEmptyArrayValidator));
configureComponent(INPUT_TYPES.SINGLE_DROPDOWN_MULTICHECKBOX, SingleDropdownMulticheckbox, NotEmptyValidator);
configureComponent(INPUT_TYPES.MULTIPLE_MULTISELECT_STRING_ALL, MultipleInput(MultiselectStringInput), AllInnerValuesNotEmptyValidator);

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

export const dontDownloadForOutputTypes = ['default', OUTPUT_TYPES.LINK]
export const dataAndUrlForOutputTypes = [OUTPUT_TYPES.PLOT, OUTPUT_TYPES.PNG]

export const OUTPUTS = {
  default: LinkOutput,
  [OUTPUT_TYPES.LINK]: LinkOutput,
  [OUTPUT_TYPES.PLOTLY]: PlotlyOutput,
  [OUTPUT_TYPES.PLOT]: PlotOutput,
  [OUTPUT_TYPES.PNG]: PngOutput,
  [OUTPUT_TYPES.CONSIGNMENT]: ConsignmentOutput,
  [OUTPUT_TYPES.RAW]: RawOutput
};