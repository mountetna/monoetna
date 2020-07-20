import React, {useState, useEffect} from 'react';
import {CheckBox} from "./CheckBox";
import CollapsibleArrow from "./CollapsibleArrow";
import {trace} from "../utils/asserts";

require('./TreeView.css')

const parentSymbol = Symbol('parent');

export default function TreeView({
  options = [],
  selected: selectedProp,
  collapsible = false,
  onSelectionsChange = (newState) => null,
  ItemLabelView = TreeViewDefaultItemLabelView
} = {}) {
  const [selectedState, setSelected] = useState(() => initializeState(options));

  useEffect(() => {
    if (selectedProp == null) setSelected(initializeState(options));
  }, [options, selectedProp])

  const selected = selectedProp == null ? selectedState : selectedProp;

  return <div className='etna-tree-view'>
    {Nodes(options, selected)}
  </div>;

  function Nodes(options, selected, path = []) {
    const inner = (node, options, selected, path) => WrapNode(
      <ItemLabelView node={node} path={path} onChange={() => select(node, path, options, selected)}
                     checked={isSelected(options, selected)}/>,
      Nodes(options, selected, path.concat([node])));

    return options.map(([node, childOptions = []]) => <div className='child' key={node}>
      {inner(node, childOptions, selected[node], path)}
    </div>)
  }

  function WrapNode(label, children) {
    return collapsible ? <CollapsibleArrow label={label}>{children}</CollapsibleArrow> : <div>
      <div>{label}</div>
      {children}
    </div>
  }

  function isSelected(options, selected) {
    if (options.length === 0) return selected === true;
    return !options.find(([optionNode, childOptions = []]) => !isSelected(childOptions, selected[optionNode]))
  }

  function select(node, path, nodeOptions, nodeSelected) {
    const newState = { ...selected };
    const newSelected = path.reduce((s, p) => (s[p] = { ...s[p] }), newState);
    newSelected[node] = initializeState(nodeOptions, nodeSelected, !isSelected(nodeOptions, nodeSelected));


    setSelected(newState);
    onSelectionsChange(newState);
  }
}

export function TreeViewDefaultItemLabelView({ node, onChange, checked, path }) {
  return <span>
    <CheckBox data-node={node} data-path={path} onChange={onChange} checked={checked}/>
    {node}
  </span>
}

function initializeState(options, parent = undefined, leafValue = true) {
  if (options.length === 0 && parent !== undefined) return leafValue;
  return options.reduce((o, [node, options = []]) => (o[node] = initializeState(options, o, leafValue), o), {});
}

export function getSelectedLeaves(selectedState) {
  return Object.keys(selectedState).reduce((r, k) => r.concat(
    selectedState[k] === true ? [k] : getSelectedLeaves(selectedState[k])), []);
}