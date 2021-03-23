// Status
export const STATUS = {
  PENDING: 'pending',
  COMPLETE: 'complete',
  ERROR: 'error',
  RUNNING: 'running',
};

// Types
export const TYPE = {
  INTEGER: 'int',
  FLOAT: 'float',
  BOOL: 'boolean',
  STRING: 'string',
  ARRAY: 'array',
  FILE: 'File',
  MULTISELECT_STRING: 'multiselect-string',
  SELECT_AUTOCOMPLETE: 'select-autocomplete',
  CHECKBOXES: 'checkboxes'
};

// CWL Step RUN Sentinels
export const RUN = {
  UI_QUERY: 'ui-queries/',
  UI_OUTPUT: 'ui-outputs/'
};

// UI Output widgets
export const OUTPUT_COMPONENT = {
  LINK: 'link',
  PLOTLY: 'plotly',
  CONSIGNMENT: 'consignment',
  RAW: 'raw'
};
