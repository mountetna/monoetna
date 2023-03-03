import React from 'react';

const ModelNode = ({model_name, center, size, selected, handler, disabled, numSiblings}) =>

  center ? (
    <div
      className={`model_node ${selected ? 'selected' : ''} ${
        disabled ? 'disabled' : ''
      }`}
      style={{top: center.y, left: center.x, '--model-width': `${600 / (numSiblings + 1) - 15}px` }}
      onClick={() => {
        if (!disabled) handler(model_name);
      }}
    >
      {model_name}
    </div>
  ) : null;

export default ModelNode;
