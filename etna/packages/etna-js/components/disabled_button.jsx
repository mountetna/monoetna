import React from 'react';

export default function DisabledButton({
  id,
  disabled,
  label,
  onClick,
  className
}) {
  let classes = disabled ? 'button disabled' : 'button';

  if (className) {
    classes = `${classes} ${className}`;
  }
  return (
    <input
      id={id}
      type='button'
      className={classes}
      value={label}
      disabled={disabled}
      onClick={onClick}
    />
  );
}
