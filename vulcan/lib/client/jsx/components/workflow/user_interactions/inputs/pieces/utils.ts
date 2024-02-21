// Value container: individual ui elements of composite inputs are keyed, and with values stored per each key.
export type keyedValues = {[k: string]: any};

// Expected structure of multi-category options is a nested Object where the categories are keys containing an Object as their value, and literal options are the leaf keys which have 'null' as value.  Thus, all keys upstream of an leaf key reflect the category path the user would take to get to that option.
// **Important caveat: There are expected to be no duplicate leaf options!
export type nestedOptionSet = {[k: string]: null | nestedOptionSet};

// Some (many) inputs should be built to work for options provided as either 'nestedOptionSet' above OR a simple string array
export type OptionSet = string[] | nestedOptionSet;
