import React from 'react';

import {Draggable} from 'react-beautiful-dnd';

import QueryColumnSelector from './query_column_selector';

const DraggableQueryColumnSelector = (props: any) => {
  return (
    <Draggable
      draggableId={`column-${props.columnIndex}`}
      index={props.columnIndex}
    >
      {(provided: any) => (
        <div
          ref={provided.innerRef}
          {...provided.draggableProps}
          {...provided.dragHandleProps}
        >
          <QueryColumnSelector {...props} />
        </div>
      )}
    </Draggable>
  );
};

export default DraggableQueryColumnSelector;
