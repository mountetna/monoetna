import React, {useEffect, useMemo} from 'react';
import {DataEnvelope} from './input_types';
import {Paper} from '@material-ui/core';
import {some} from '../../../../selectors/maybe';
import SelectAutocompleteInput from './select_autocomplete';
import {DragDropContext, Droppable, Draggable} from 'react-beautiful-dnd';
import {arrayLevels} from './user_input_pieces';

const LevelComponent = (props: any) => {
  return (
    <Draggable
      draggableId={`level-${props.levelIndex}`}
      index={props.levelIndex}
    >
      {(provided: any) => (
        <div
          ref={provided.innerRef}
          {...provided.draggableProps}
          {...provided.dragHandleProps}
        >
          {props.level}
        </div>
      )}
    </Draggable>
  );
};

/*
reorderPiece is configured around the VizualizationUI.
It could be generalized fairly easily, if needed.
*/
export function ReorderPiece(
  key: string,
  changeFxn: Function,
  value: string | string[] = 'unordered',
  label: string,
  full_data: DataEnvelope<any[]>,
  data_target: string | null,
  discrete_data: string[]
) {
  const handleOnDragEnd = (result: any) => {
    if (!result.destination) {
      return;
    }

    if (result.destination.index === result.source.index) {
      return;
    }

    const newValues = Array.from(value);
    const [removed] = newValues.splice(result.source.index, 1);
    newValues.splice(result.destination.index, 0, removed);

    changeFxn(newValues, key);
  };
  const onMethodSelect = (e: any, picked: string | null) => {
    const newValue =
      picked == 'custom'
        ? arrayLevels(Object.values(full_data[data_target as string]))
        : picked;
    changeFxn(newValue, key);
  };

  const chosen_case = Array.isArray(value) ? 'custom' : value;
  const order_options = useMemo(() => {
    if (data_target != null && discrete_data.includes(data_target)) {
      return ['increasing', 'decreasing', 'unordered', 'custom'];
    }
    return ['increasing', 'decreasing', 'unordered'];
  }, [data_target, discrete_data]);

  // Reset to 'increasing' if data_target changes while 'custom', (also hit on page refresh!)
  useEffect(() => {
    if (data_target != null && chosen_case == 'custom') {
      // Needed if new data_target is discrete or if new data levels aren't captured
      let needs_reset = !discrete_data.includes(data_target);
      if (!needs_reset) {
        needs_reset =
          arrayLevels(Object.values(full_data[data_target as string])).filter(
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

  const reorder_custom = !Array.isArray(value) ? null : (
    <DragDropContext onDragEnd={handleOnDragEnd}>
      <Droppable droppableId='columns'>
        {(provided: any) => (
          <Paper ref={provided.innerRef} {...provided.droppableProps}>
            {(value as string[]).map((level: string, index: number) => {
              return (
                <LevelComponent key={index} level={level} levelIndex={index} />
              );
            })}
            {provided.placeholder}
          </Paper>
        )}
      </Droppable>
    </DragDropContext>
  );

  return (
    <div key={key} style={{paddingTop: 8}}>
      {case_dropdown}
      {reorder_custom}
    </div>
  );
}
