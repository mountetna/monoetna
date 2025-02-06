import React, {useState, useEffect, useMemo} from 'react';
import {DataEnvelope} from '../../input_types';
import Paper from '@material-ui/core/Paper';
import FormControl from '@material-ui/core/FormControl';
import InputLabel from '@material-ui/core/InputLabel';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';
import ArrowDropDownIcon from '@material-ui/icons/ArrowDropDown';
import ArrowDropUpIcon from '@material-ui/icons/ArrowDropUp';
import LowPriorityIcon from '@material-ui/icons/LowPriority';
import DragHandleIcon from '@material-ui/icons/DragHandle';
import Typography from '@material-ui/core/Typography';
import {some} from '../../../../../selectors/maybe';
import SelectAutocompleteInput from '../components/select_autocomplete';
import {DragDropContext, Droppable, Draggable} from 'react-beautiful-dnd';
import {arrayLevels, checkboxPiece} from './user_input_pieces';

const LevelComponent = (props: any) => {
  return (
    <Draggable
      draggableId={`level-${props.levelIndex}`}
      index={props.levelIndex}
    >
      {(provided: any) => (
        <Grid container
          ref={provided.innerRef}
          {...provided.draggableProps}
          {...provided.dragHandleProps}
        >
          <Grid item>
            <DragHandleIcon fontSize='small' color='secondary'/>
          </Grid>
          <Grid item>
            <Typography component='text' fontSize='8'>{props.level}</Typography>
          </Grid>
        </Grid>
      )}
    </Draggable>
  );
};

const onDragEnd = (fullValue: string | string[], changeFxn: Function) => {
  return (result: any)  => {
    if (!result.destination) {
      return;
    }

    if (result.destination.index === result.source.index) {
      return;
    }

    const newValues = Array.from(fullValue);
    const [removed] = newValues.splice(result.source.index, 1);
    newValues.splice(result.destination.index, 0, removed);

    changeFxn(newValues);
  };
};

const DragDrop = (values: string[], HandleOnDragEnd: ReturnType<typeof onDragEnd>, styling: Object = {}) => {
return <DragDropContext onDragEnd={HandleOnDragEnd}>
  <Droppable droppableId='columns'>
    {(provided: any) => (
      <Paper ref={provided.innerRef} {...provided.droppableProps}
        style={styling}>
        {(values).map((level: string, index: number) => {
          return (
            <LevelComponent key={index} level={level} levelIndex={index} />
          );
        })}
        {provided.placeholder}
      </Paper>
    )}
  </Droppable>
</DragDropContext>
};

export function ReorderCollapsiblePiece(
  key: string,
  changeFxn: Function,
  value: string[] | null = [] as string[],
  label: string = 'Reorder'
) {
  const [open, setOpen] = useState(false);
  const value_use = value==null ? [] : value;
  const disabled=value_use.length<2;

  const openToggle = <Grid item>
    <IconButton aria-label='open-close' size='small' color='secondary' disabled={disabled} onClick={()=>setOpen(!open)}>
      <Grid item container>
        <Grid item>
          <LowPriorityIcon fontSize='small'/>
        </Grid>
        <Grid item>
          {open ?
            <ArrowDropUpIcon fontSize='small'/> :
            <ArrowDropDownIcon fontSize='small'/>
          }
        </Grid>
      </Grid>
    </IconButton>
  </Grid>;

  return (
    <Grid key={key} item container direction='row'>
      {openToggle}
      <Grid item style={{paddingTop: '8px'}}>
        <FormControl>
          <InputLabel shrink disabled={disabled}
            style={{width: 'max-content', userSelect: 'none'}}
            onClick={()=>{if(!disabled)setOpen(!open)}}
          >
            {label}
          </InputLabel>
          {open && value_use.length>=2 ? DragDrop(value_use, onDragEnd(value_use, changeFxn), {paddingTop: '15px'}) : null}
        </FormControl>
      </Grid>
    </Grid>
  );
};

export function ReorderVizPiece(
  key: string,
  changeFxn: Function,
  value: string | string[] = 'unordered',
  label: string,
  full_data: DataEnvelope<any[]>,
  data_target: string | string[] | null,
  discrete_data: string[],
  can_randomize: boolean = false
) {
  const data_targ = Array.isArray(data_target) ? data_target[data_target.length - 1] : data_target
  const onMethodSelect = (e: any, picked: string | null) => {
    const newValue =
      picked == 'custom'
        ? arrayLevels(Object.values(full_data[data_targ as string]))
        : picked;
    changeFxn(newValue, key);
  };

  const chosen_case = Array.isArray(value) ? 'custom' : value;
  const order_options = useMemo(() => {
    let options = ['increasing', 'decreasing', 'unordered']
    if (can_randomize) {
      options.push('randomize')
    }
    if (data_targ != null && discrete_data.includes(data_targ)) {
      options.push('custom');
    }
    return options;
  }, [data_target, discrete_data, can_randomize]);

  // Reset to 'increasing' if data_target changes while 'custom', (also hit on page refresh!)
  useEffect(() => {
    if (data_targ != null && chosen_case == 'custom') {
      // Needed if new data_targ is discrete or if new data levels aren't captured
      let needs_reset = !discrete_data.includes(data_targ);
      if (!needs_reset) {
        needs_reset =
          arrayLevels(Object.values(full_data[data_targ as string])).filter(
            (val) => !value.includes(val)
          ).length > 0;
      }
      if (needs_reset) onMethodSelect(null, 'increasing');
    }
  }, [full_data, data_target, discrete_data]);

  const case_dropdown = (
    <SelectAutocompleteInput
      key={key}
      label={label}
      value={some(chosen_case)}
      data={{a: order_options}}
      onChange={() => {}}
      onChangeOverride={onMethodSelect}
    />
  );

  return (
    <div key={key} style={{paddingTop: 8}}>
      {case_dropdown}
      {!Array.isArray(value) ? null : DragDrop(value, onDragEnd(value, (r:any) => changeFxn(r, key)))}
    </div>
  );
}

export function ReorderCustomOnlyPiece(
  key: string,
  changeFxn: Function,
  value: string | string[] = 'make',
  label: string,
  full_data: DataEnvelope<any[]>,
  data_target: string | string[] | null,
  discrete_data: string[]
) {

  const data_targ = useMemo( () => {
    if (data_target == null) return null
    return Array.isArray(data_target) ? data_target[data_target.length - 1] : data_target
  }, [data_target])
  const levels = useMemo( () => {
    if (data_targ == null || !discrete_data.includes(data_targ)) return null
    return arrayLevels(Object.values(full_data[data_targ as string]))
  }, [full_data, data_targ, discrete_data])
  const canReorder = useMemo( () => {
    return data_targ != null && discrete_data.includes(data_targ)
  }, [data_targ, discrete_data])

  const startOrClear = (doReorder: boolean, x?: any) => {
    const new_full = doReorder ? levels : 'make'
    changeFxn(new_full, key);
  };

  // Reset to 'off'-mode if data_target changes, (also hit on page refresh, so ultimate check should be levels validity!)
  useEffect(() => {
    if (data_targ != null && value != 'make') {
      // Needed if new data_targ is discrete or if new data levels aren't captured
      let needs_reset = !discrete_data.includes(data_targ);
      if (!needs_reset) {
        needs_reset =
          arrayLevels(Object.values(full_data[data_targ as string])).filter(
            (val) => !value.includes(val)
          ).length > 0;
      }
      if (needs_reset) startOrClear(false);
    }
  }, [full_data, data_target, discrete_data]);

  return (
    <div key={key} style={{paddingTop: 8}}>
      {checkboxPiece(
        key,
        startOrClear,
        value != 'make',
        label,
        !canReorder
      )}
      {!Array.isArray(value) ? null : DragDrop(value, onDragEnd(value, (r:any) => changeFxn(r, key)))}
    </div>
  );
}

