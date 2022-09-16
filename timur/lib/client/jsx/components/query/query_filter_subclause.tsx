import React, {useMemo, useCallback, useState, useEffect} from 'react';
import Grid from '@material-ui/core/Grid';
import TextField from '@material-ui/core/TextField';
import FormControl from '@material-ui/core/FormControl';
import Autocomplete from '@material-ui/lab/Autocomplete';
import {makeStyles} from '@material-ui/core/styles';

import {Debouncer} from 'etna-js/utils/debouncer';
import {QuerySubclause} from '../../contexts/query/query_types';
import FilterOperator from './query_filter_operator';
import useQuerySubclause from './query_use_query_subclause';
import {QueryGraph} from '../../utils/query_graph';
import RemoveIcon from './query_remove_icon';
import Selector from './query_selector';
import {Attribute} from '../../models/model_types';

const useStyles = makeStyles((theme) => ({
  option: {
    width: 'max-content',
    whiteSpace: 'nowrap'
  },
  listbox: {
    width: 'max-content'
  },
  paper: {
    width: 'max-content'
  },
  operand: {
    paddingTop: '3px'
  }
}));

const QueryFilterSubclause = ({
  subclause,
  subclauseIndex,
  modelName,
  modelAttributes,
  graph,
  waitTime,
  eager,
  isColumnFilter,
  showRemoveIcon,
  patchSubclause,
  removeSubclause
}: {
  subclause: QuerySubclause;
  subclauseIndex: number;
  modelName: string;
  modelAttributes: Attribute[];
  graph: QueryGraph;
  waitTime?: number;
  eager?: boolean;
  isColumnFilter: boolean;
  showRemoveIcon: boolean;
  patchSubclause: (subclause: QuerySubclause) => void;
  removeSubclause: () => void;
}) => {
  const [operandValue, setOperandValue] = useState('' as string | number);
  const [previousOperandValue, setPreviousOperandValue] = useState(
    '' as string | number
  );
  const [previousAttributeName, setPreviousAttributeName] = useState(
    subclause.attributeName
  );
  const [debouncer, setDebouncer] = useState(
    () => new Debouncer({windowMs: waitTime, eager})
  );
  const classes = useStyles();

  // Clear the existing debouncer and accept any new changes to the settings
  useEffect(() => {
    const debouncer = new Debouncer({windowMs: waitTime, eager});
    setDebouncer(debouncer);
    return () => debouncer.reset();
  }, [waitTime, eager]);

  const {attributeType, distinctAttributeValues, fetchDistinctAttributeValues} =
    useQuerySubclause({
      subclause,
      graph,
      modelName
    });

  const filterOperator = useMemo(() => {
    return new FilterOperator({
      subclause,
      isColumnFilter
    });
  }, [subclause, isColumnFilter]);

  useEffect(() => {
    // When user selects a different attribute, update the type
    if (attributeType !== subclause.attributeType) {
      patchSubclause({
        ...subclause,
        attributeType
      });
    }
  }, [attributeType, subclause, patchSubclause]);

  useEffect(() => {
    // When component loads, if subclause already populated then
    //    fetch the pre-selected attribute values.
    if (
      '' !== subclause.attributeName &&
      filterOperator.hasPrepopulatedOperandOptions()
    ) {
      fetchDistinctAttributeValues();
    }
  }, []);

  useEffect(() => {
    // When user selects a different attribute, update the pre-populated options
    if (
      previousAttributeName !== subclause.attributeName &&
      '' !== subclause.attributeName &&
      filterOperator.hasPrepopulatedOperandOptions()
    ) {
      setPreviousAttributeName(subclause.attributeName);
      fetchDistinctAttributeValues();
    }
  }, [
    subclause.attributeName,
    filterOperator,
    previousAttributeName,
    fetchDistinctAttributeValues
  ]);

  const handleAttributeSelect = useCallback(
    (attributeName: string) => {
      patchSubclause({
        ...subclause,
        attributeName,
        operator: '',
        operand: ''
      });
    },
    [subclause, patchSubclause]
  );

  const handleOperatorSelect = useCallback(
    (operator: string) =>
      patchSubclause({
        ...subclause,
        operator: filterOperator.magmify(operator)
      }),
    [subclause, patchSubclause, filterOperator]
  );

  const handleOperandChange = useCallback(
    (operand: string) => {
      patchSubclause({
        ...subclause,
        operand: filterOperator.formatOperand(operand)
      });
    },
    [patchSubclause, subclause, filterOperator]
  );

  const handleOperandChangeWithDebounce = useCallback(
    (value: string) => {
      debouncer.ready(() => handleOperandChange(value));
      setOperandValue(value);
    },
    [handleOperandChange, debouncer]
  );

  // When the operand value changes, follow it
  useEffect(() => {
    if (subclause.operand !== previousOperandValue) {
      debouncer.reset();
      setOperandValue(subclause.operand);
      setPreviousOperandValue(subclause.operand);
    }
  }, [subclause.operand, debouncer, previousOperandValue]);

  let uniqId = (idType: string): string =>
    `${idType}-Select-${Math.random().toString()}`;

  return (
    <Grid container>
      <Grid item xs={showRemoveIcon ? 4 : 5}>
        {modelAttributes.length > 0 ? (
          <Selector
            canEdit={true}
            label='attribute'
            name={subclause.attributeName}
            onSelect={handleAttributeSelect}
            choiceSet={modelAttributes.map((a) => a.attribute_name)}
          />
        ) : null}
      </Grid>
      <Grid item xs={3}>
        <Selector
          label={`operator-${subclauseIndex}`}
          canEdit={true}
          name={filterOperator.prettify() || ''}
          choiceSet={Object.keys(filterOperator.options())}
          onSelect={handleOperatorSelect}
        />
      </Grid>
      <Grid item xs={4}>
        {filterOperator.hasOperand() ? (
          <FormControl fullWidth>
            {filterOperator.hasPrepopulatedOperandOptions() &&
            distinctAttributeValues.length > 0 ? (
              <Autocomplete
                classes={{
                  option: classes.option,
                  listbox: classes.listbox,
                  paper: classes.paper,
                  root: classes.operand
                }}
                id={uniqId(`operand-${subclauseIndex}`)}
                freeSolo
                fullWidth
                options={distinctAttributeValues}
                renderInput={(params) => <TextField {...params} />}
                onInputChange={(e, v, r) => {
                  // Only send event if user manually clears the value
                  //   or selects a non-empty-string option.
                  if ('' !== v || 'reset' !== r)
                    handleOperandChangeWithDebounce(v || '');
                }}
                inputValue={operandValue.toString()}
                data-testid='operand-autocomplete'
              />
            ) : (
              <TextField
                id={uniqId(`operand-${subclauseIndex}`)}
                value={operandValue}
                onChange={(e) =>
                  handleOperandChangeWithDebounce(e.target.value as string)
                }
              />
            )}
          </FormControl>
        ) : null}
      </Grid>
      {showRemoveIcon ? (
        <Grid item xs={1} container justify='flex-end'>
          <RemoveIcon
            showRemoveIcon={showRemoveIcon}
            onClick={removeSubclause}
            label='subclause'
          />
        </Grid>
      ) : null}
    </Grid>
  );
};
export default QueryFilterSubclause;
