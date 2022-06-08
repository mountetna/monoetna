import React from 'react';
import { DataEnvelope } from './input_types';
import { Paper } from '@material-ui/core';
import { some } from "../../../../selectors/maybe";
import SelectAutocompleteInput from './select_autocomplete';
import { DragDropContext, Droppable, Draggable } from 'react-beautiful-dnd';
import { arrayLevels } from './subsetDataFrame_piece';

const LevelComponent = (props: any) => {
  return(
    <Draggable
      draggableId={`level-${props.levelIndex}`}
      index={props.levelIndex}
    >
      {(provided: any) => (
        <div
          ref={provided.innerRef}
          {...provided.draggableProps}
          {...provided.dragHandleProps}
          style={{padding:2}}
        >
          {props.level}
        </div>
      )}
    </Draggable>
  )
}

export function reorderPiece(
  key: string, changeFxn: Function, value: string | string[] = "unordered",
  label: string, full_data: DataEnvelope<any[]>, data_target: string | null) {

  const handleOnDragEnd = (result: any) => {
    if (!result.destination) {
      return;
    }

    if (result.destination.index === result.source.index) {
      return;
    }

    const newValues = Array.from(value)
    const [removed] = newValues.splice(result.source.index, 1);
    newValues.splice(result.destination.index, 0, removed);

    changeFxn(newValues, key)
  }

  const chosen_case = Array.isArray(value) ? "custom" : value

  const case_dropdown = <SelectAutocompleteInput
      key={key}
      label={label}
      disabled={data_target == null}
      value={some(chosen_case)}
      data={{a: ['unordered', 'increasing', 'decreasing', 'custom']}}
      onChange={() => {}}
      onChangeOverride={ (e, picked) => {
        console.log({picked})
        const newValue = picked=='custom' ? arrayLevels(Object.values(full_data[data_target as string])) : picked
        changeFxn(newValue, key)
      }}
    />
  
  const reorder_custom = !Array.isArray(value) ? null :
    <DragDropContext onDragEnd={handleOnDragEnd} style={{paddingLeft:10, padding:2}}>
      <Droppable droppableId='columns'>
        {(provided: any) => (
          <Paper
            ref={provided.innerRef}
            {...provided.droppableProps}
            // style={{display: 'inline-flex'}}
          >
            {(value as string[]).map((level: string, index: number) => {
              return (
                <LevelComponent
                  key={index}
                  level={level}
                  levelIndex={index}
                />
              );
            })}
            {provided.placeholder}
          </Paper>
        )}
      </Droppable>
    </DragDropContext>
  
  return(
    <div key={key} style={{paddingTop: 8}}>
      {case_dropdown}
      {reorder_custom}
    </div>
  )
}