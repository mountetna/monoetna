import React, {useMemo, useCallback, useState, useEffect, useContext} from 'react';
import Grid from '@material-ui/core/Grid';
import TextField from '@material-ui/core/TextField';
import FormControl from '@material-ui/core/FormControl';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import Autocomplete from '@material-ui/lab/Autocomplete';

import {Debouncer} from 'etna-js/utils/debouncer';
import {QuerySubclause} from '../../contexts/query/query_types';
import FilterOperator from './query_filter_operator';
import useQuerySubclause from './query_use_query_subclause';
import useQueryClause from './query_use_query_clause';
import {QueryGraph} from '../../utils/query_graph';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import QueryNumber from './query_number';
import RemoveIcon from './query_remove_icon';
import Selector from './query_selector';
import {Attribute} from 'etna-js/models/magma-model';

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
  filter_operand: {
    flex: '1'
  },
  operand: {
    paddingTop: '3px'
  }
}));

const QueryFilterSubclause = ({
  subclause,
  subclauseIndex,
  modelName,
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
  waitTime?: number;
  eager?: boolean;
  isColumnFilter: boolean;
  showRemoveIcon: boolean;
  patchSubclause: (subclause: QuerySubclause) => void;
  removeSubclause: () => void;
}) => {
  const { state: {graph} } = useContext(QueryGraphContext);

  const {modelAttributes} = useQueryClause({
    modelName, graph, isColumnFilter
  });

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
        attributeType: '',
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

  const [ removeHint, setRemoveHint ] = useState(false);

  return <Grid item container alignItems='center' style={{ textDecoration: removeHint ? 'line-through' : 'none' }} >
    {
      !isColumnFilter && <>
        <QueryNumber
          setRemoveHint={ showRemoveIcon ? setRemoveHint : undefined }
          onClick={ showRemoveIcon ? removeSubclause : undefined}
          number={subclauseIndex}
          level={2}/>
        <Selector
          canEdit={true}
          label='attribute'
          placeholder='attribute'
          color='secondary'
          name={subclause.attributeName}
          onSelect={handleAttributeSelect}
          choiceSet={modelAttributes.map((a) => a.attribute_name)}
        />
      </>
    }
    <Selector
      label={`operator-${subclauseIndex}`}
      canEdit={true}
      name={filterOperator.prettify() || ''}
      placeholder='satisfies'
      color='purple'
      choiceSet={Object.keys(filterOperator.options())}
      onSelect={handleOperatorSelect}
    />
    {
      filterOperator.hasOperand() &&
      <FormControl className={ classes.filter_operand } variant="standard">
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
            placeholder='value'
            options={distinctAttributeValues}
            renderInput={(params) => <TextField variant="standard" {...params} />}
            onInputChange={(e, v, r) => {
              // Only send event if user manually clears the value
              //   or selects a non-empty-string option.
              if ('' !== v || 'reset' !== r) {
                handleOperandChangeWithDebounce(v || '');
              }
            }}
            inputValue={operandValue.toString()}
            data-testid='operand-autocomplete'
          />
        ) : (
          <TextField
            variant="standard"
            id={uniqId(`operand-${subclauseIndex}`)}
            value={operandValue}
            fullWidth
            placeholder='value'
            onChange={(e) =>
              handleOperandChangeWithDebounce(e.target.value as string)
            } />
        )}
      </FormControl>
    }
  </Grid>;
};
export default QueryFilterSubclause;
