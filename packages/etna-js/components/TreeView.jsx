import React, {useState, useEffect} from 'react';

export default function TreeView({ options = [], selected = [], onChange = (path) => null, ItemView = TreeViewDefaultItemView } = {}) {
  const [selectedState, setSelected] = useState([]);
  // selected = selected == null ? selectedState : selected;

  return <div>
    Hello
  </div>;
}

function TreeViewInner(options, selected, onChange, ItemView) {

}

export function TreeViewDefaultItemView({ children, path }) {
  return <div>
    { children }
  </div>;
}